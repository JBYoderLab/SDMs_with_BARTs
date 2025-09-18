# Bayesian regression, a worked example
# Jeremy B. Yoder, 18 Sept 2025

# In the last few decades, biologists have increasingly adopted Bayesian 
# approaches, which can accommodate more complex "prior" understandings of the 
# world --- but which are also more computationally complex, and often require 
# you to make explicit choices about those priors. Bayesian methods are 
# popular in phylogenetics and population genetics, but we can get a taste
# of this philosophical choice in "simple" linear regression.

# Here, we'll use anexample data set describing plants' floral symmetry, floral
# visitor diversity, and other ecological factors (from Yoder et al 2020,
# Biology Letters). The visitor species count data are highly non-normal --- 
# and because the analysis involves a comparison among species, we really should
# account for the counfounding effect of phylogenetic relationships among those
# species. The solution Yoder et al chose was Bayesian regression, using 
# the `brms` package in R.

# This script assumes your working directory is the project directory; you may want to change this
# setwd("~/Documents/Active_projects/SDMs_with_BARTs")
rm(list=ls()) 

# Load necessary packages 
library("tidyverse")
library("fitdistrplus") # Useful for describing the distribution of data
# install.packages("fitdistrplus")
library("brms") # For running Bayesian regression

set.seed(20250925)

#-------------------------------------------------------------------------
# Load the flower visitor data

visitors <- read.csv("data/flower_visitors_all.csv")
glimpse(visitors)

# This is a somewhat complex data set --- it includes a bunch of different 
# measures of pollinator diversity and pollinator sharing between co-occurring 
# plant species, and loadings on four PC axes derived from phylogenetic
# relationships among the plant species in the data set, AND the latitude of the 
# community in which the plant species was observed. The latitude is relevant
# because species diversity is higher in the tropics, so we may want to include
# it as a covariate in any model involving species diversity, like the visitor
# species counts. Column names with the `.sc` suffix have been rescaled for 
# convenience; we'll use the scaled latitude and PC axes.

hist(visitors$n.poll)

#-------------------------------------------------------------------------
# First, let's visualize our predictors, just to get a sense of what's going on:

# effect of floral symmetry?
ggplot(visitors, aes(x=symmetry, y=log10(n.poll))) + 
  geom_boxplot() + 
  theme_bw()

# effect of latitude?
ggplot(visitors, aes(x=lat.sc, y=log10(n.poll))) + 
  geom_point() + 
  geom_smooth(method="lm") + 
  theme_bw()

# effects of phylogeny?
ggplot(visitors, aes(x=pc1.sc, y=log10(n.poll))) + 
  geom_point() + 
  geom_smooth(method="lm") +
  theme_bw()

ggplot(visitors, aes(x=pc2.sc, y=log10(n.poll))) + 
  geom_point() + 
  geom_smooth(method="lm") +
  theme_bw()

# Looks like many of these have some potential contribution, but we'll need to
# test them formally, ideally all together.

#-------------------------------------------------------------------------
# Choose an error distribution

# As we've seen, visitor species number is highly right-skewed:
hist(visitors$n.poll)

# The `descdist()` function gives us a visual assessment of distributions that
# might be a decent fit for this data, by bootstrapping the data:
descdist(visitors$n.poll, boot=100) 
# Looks like beta or gamma will be pretty good.

#-------------------------------------------------------------------------
# Perform Bayesian regression

# Now we'll use `brms` to fit some candidate models. For all models, we want to
# use a random effect of community ID (in this version of the data, the `web` 
# column) to control for the fact that each community has different baseline 
# diversity of plants and flower-visiting animals.

# We'll fit three candidate models. For each of these, the `brm()` function 
# essentially creates a little C++ program to simulate data as part of fitting
# the model --- this is where priors come in. We mostly won't mess with the 
# default priors, but we do need to specify an error distribution, which 
# we do in a very similar way to frequentist generalized modeling. 

# Based on our check above, we'll use a gamma distribution. The `iter` option 
# specifies how long to run the simulation, and `warmup` specifies how many 
# iterations to treat as "burn-in" before retaining iterated parameters to 
# estimate posterior values.

# The `cores` option tells `brm` to spread the work of simulating over four 
# independent processor cores, which speeds things up. 

# For a first attempt, we'll fit our most complete model, accounting for 
# possible effects of floral symmetry, latitude, phylogeny, and an interaction
# between latitude and symmetry. This may take substantial processing time, but 
# you'll get progress readouts to track the simulations.
npoll.MSxLP <- brm(n.poll ~ symmetry * lat.sc + pc1.sc + pc2.sc + (1|web), data=visitors, family="gamma", iter=5000, warmup=3000, cores=4)

# This is a point at which Bayesian inference clearly differs from frequentist
# methods. Rather than assessing whether estimated effects differ from zero (are
# "significant") by conducting a t-test of the estimate, or assessing whether
# a predictor explains significant variation using an F-test, we look at the 
# posterior distribution of estimates from the simulation underlying the model.

summary(npoll.MSxLP)

# Each parameter estimate is accompanied by the 95% CI derived from the 
# simulation procedure --- an effect is "significant" in this framework if the
# 95% CI doesn't cross zero. You can see that there's a signficant negative
# effect of zygomorphic floral symmetry, and a significant positive effect of
# latitude, for instance.

# If we plot the model, we can get visualizations of the posterior parameter 
# distributions overall, and in time-series
plot(npoll.MSxLP)

# The time-series figures let us confirm that the simulation has sampled from a
# "stationary" distribution, that is, from a range of values close to an optimum.
# If the time-series plots don't show clear trends, or look like "fuzzy 
# caterpillars", they're probably moving mostly at random around an optimal 
# value, and the variation in the posterior values they sample is only due to 
# the random element of the simulation. This is also reflected in the ESS values 
# in the summary readout --- a bigger ESS, and Rhat close to 1, means more 
# statistically independent samples were taken after warmup, and the posterior
# estimates are reliable. (The brm function will throw a warning if ESS is too 
# low to be reliable for estimation.)

# If we did see trends in the time series, or poor ESS values in the summary,
# we'd know the estimates were poor representations of the posterior, and we'd
# revise the simulation to run longer, and set the `warmup` parameter longer,
# to try to give the MCMC more time to find a stable optimum.

# We can also fit simpler models, for comparison. We'll try the original model
# without the interaction between symmetry and latitude, without phylogeny, and
# with just the random effect of community/web. 
npoll.MSLP <- brm(n.poll ~ symmetry + lat.sc + pc1.sc + pc2.sc + (1|web), data=visitors, family="gamma", iter=5000, warmup=3000, cores=4)
npoll.MSL <- brm(n.poll ~ symmetry + lat.sc + (1|web), data=visitors, family="gamma", iter=5000, warmup=3000, cores=4)
npoll.M <- brm(n.poll ~ (1|web), data=visitors, family="gamma", iter=5000, warmup=3000, cores=4)

# To compare models fitted in this way, we use a procedure called "leave one 
# out" (LOO) cross-validation. This is (very roughly) analogous to comparing AIC
# scores for alternative models ---  
npoll.comp <- LOO(npoll.M, npoll.MSL, npoll.MSLP, npoll.MSxLP) 

npoll.comp 
# Our best-fit model (with the lowest ELPD, or expected log pointwise predictive 
# density) is MSxLP, the most complex model.

# We can also get an estimate of the variation explained using this function:
bayes_R2(npoll.MSxLP) # rsq = 0.28




