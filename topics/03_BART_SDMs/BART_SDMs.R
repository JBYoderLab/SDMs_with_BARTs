# Species distribution modeling with BARTs in R
# Jeremy B. Yoder, 22 Sept 2025

# Clears the environment and load key packages
rm(list=ls())

# This script assumes your working directory is the project directory; you may want to change this
# setwd("~/Documents/Active_projects/SDMs_with_BARTs")

library("tidyverse") 
library("sf")
library("terra")

library("rnaturalearth")
library("rnaturalearthdata")
library("geodata")

library("embarcadero")
# devtools::install_github("cjcarlson/embarcadero", build_vignettes = TRUE, force = TRUE)


set.seed(20250915)

#-----------------------------------------------------------
# Load data

# Species presence/pseudo-absence points generated previously
PA <- read.csv("output/JT_presence-pseudoabsence-envs.txt")
glimpse(PA)

# Environmental data and accessory stuff
MojExt <- c(-119.5, -112, 33, 38.5) # Extent of the Mojave desert

# climate (Bioclim, previously downloaded at max res)
BClim <- crop(worldclim_tile("bio", -120, 33, path="../data/Yucca/"), MojExt)

plot(BClim[[1]]) # one layer from the SpatRaster stack

# NB: the BioClim layers are (still) incorrectly sorted, per https://github.com/issues/created?issue=rspatial%7Cgeodata%7C38
# But that issue thread lets us know how to fix this ...
# "It seems that the order of bioclim variables is like the following: 1, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 2, 3 etc"

names(BClim) <- c("AMT", "MTWaQ", "MTCQ", "AP", "PWeM", "PDM", "PS", "PWeQ", "PDQ", "PWaQ", "PCQ", "MDR", "ITH", "TS", "MTWaM", "MTCM", "TAR", "MTWeQ", "MTDQ")

# As a reminder, the bioclim variables are averages calculated from 1970-2000;
# the 19 variables are

# BIO1 = Annual Mean Temperature (AMT)
# BIO2 = Mean Diurnal Range (Mean of monthly (max temp - min temp)) (MDR)
# BIO3 = Isothermality (BIO2/BIO7) (×100) (ITH)
# BIO4 = Temperature Seasonality (standard deviation ×100) (TS)
# BIO5 = Max Temperature of Warmest Month (MWaMT)
# BIO6 = Min Temperature of Coldest Month (MCMT)
# BIO7 = Temperature Annual Range (BIO5-BIO6) (TAR)
# BIO8 = Mean Temperature of Wettest Quarter (MWeQT)
# BIO9 = Mean Temperature of Driest Quarter (MDQT)
# BIO10 = Mean Temperature of Warmest Quarter (MWaQT)
# BIO11 = Mean Temperature of Coldest Quarter (MCQT)
# BIO12 = Annual Precipitation (AP)
# BIO13 = Precipitation of Wettest Month (PWeM)
# BIO14 = Precipitation of Driest Month (PDM)
# BIO15 = Precipitation Seasonality (Coefficient of Variation) (PS)
# BIO16 = Precipitation of Wettest Quarter (PWeQ)
# BIO17 = Precipitation of Driest Quarter (PDQ)
# BIO18 = Precipitation of Warmest Quarter (PWaQ)
# BIO19 = Precipitation of Coldest Quarter (PCQ)

#-------------------------------------------------------------------------
# Fit a species distribution model with BARTs

# vector of the Bioclim variables
xnames <- colnames(PA)[5:23]

# Fit a simple BART model with all the predictors
jtBART <- bart(y.train=PA[,"JT"], x.train=PA[,xnames], keeptrees=TRUE)

{png("topics/03_BART_SDMs/jtBART_summary.png", width=500, height=500)

summary(jtBART)

}
dev.off()

# As before, we want to simplify the model from the full set of 19 bioclim variables
# One option is to do this with stepwise BART model training. This will be slow!
jtBART.step <- bart.step(y.data=as.numeric(PA[,"JT"]), x.data=PA[,xnames], full=FALSE, quiet=TRUE)

# It may be useful to save a model object for later work. 
# First, we need to "touch" the model state to make sure the trees are saved with the rest of the model.
# This is a `dbarts` feature that is supposed to save storage.
invisible(jtBART.step$fit$state)
# Then write out the model to an Rdata file
write_rds(jtBART.step, file="output/models/jtBART.step.rds") 
jtBART.step <- read_rds(file="output/models/jtBART.step.rds")

summary(jtBART.step)

stepX <- attr(jtBART.step$fit$data@x, "term.labels")
stepX # these are the predictors kept in the stepwise model


#-------------------------------------------------------------------------
# More formal predictor selection with varimp.diag()

# Alternatively, use varimp.diag() for predictor selection. This method fits
# models of varying complexity (numbers of trees), and compares the contributions 
# of each predictor in more complex models (with more trees) versus simpler
# models (with fewer trees). In simpler models, each branch point in a tree
# is under more constraint to be "right", or informative --- a more complex model
# with more trees can "afford" to waste branches on relatively uninformative
# features. This means predictors that have larger contributions specifically in 
# simpler models are more likely to be highly informative.

# In practice, this means training a bunch of models, so this may be very slow!
jtVarimp <- varimp.diag(y.data=as.numeric(PA[,"JT"]), x.data=PA[,xnames])

# a useful fix for a minor issue
jtVarimp$data <- jtVarimp$data %>% mutate(trees = factor(trees, c(10,20,50,100,150,200)))

jtVarimp$labels$group <- "Trees"
jtVarimp$labels$colour <- "Trees"

write_rds(jtVarimp, file="output/models/jt_varimp.rds") # save the work
jtVarimp <- read_rds(file="output/models/jt_varimp.rds") # read it back in

# Write out the varimp diagnostic plot
{png(file="topics/03_BART_SDMs/jt_varimp_plot.png", width=750, height=500)
  
  jtVarimp + 
    theme_bw(base_size=18) +
    theme(legend.position="inside", legend.position.inside=c(0.8, 0.7), axis.text.x=element_text(angle=90))  # okay nice
  
}
dev.off()

# In that diagnostic plot, we're looking for the predictors that have visibly
# higher "relative contribution" in simpler models compared, as opposed to 
# predictors that have relatively similar relative contribution regardless of
# model complexity. Making the dividing line is more art than science, but 
# we can defensibly draw it at TS, temperature seasonality, and all predictors
# to its left.

# Define the top predictor set
topX <- c("PS", "PDQ", "PWaQ", "MTCM", "PDM", "MTWeQ", "TS")

# Train a new model with those top predictors
jtBARTtop <- bart(y.train=PA[,"JT"], x.train=PA[,topX], keeptrees=TRUE)

# Make sure trees are stored with the model object
invisible(jtBARTtop$fit$state)
# Then write out the model to an Rdata file
write_rds(jtBARTtop, file="output/models/jtBARTtop.rds") 
jtBARTtop <- read_rds(file="output/models/jtBARTtop.rds")


{png("topics/03_BART_SDMs/jtBARTtop_summary.png", width=500, height=500)
  
summary(jtBARTtop)
  
}
dev.off()

# Did our predictor set from varimp.diag() differ from stepwise fitting?
setdiff(stepX, topX) # stepwise fitting selected MTCQ, ITH, and TAR
setdiff(topX, stepX) # all the topX predictors were found in stepwise fitting

# Check the varimp plot again and see where the three extra predictors found
# by stepwise fitting lie in the ranking.

# Compare summary outputs for "topX" and stepwise models
summary(jtBARTtop)
summary(jtBART.step)

# Show the relative importance of predictors in a model:
varimp(jtBARTtop)

# varimp() also generates a plot, if you want:
{png("topics/03_BART_SDMs/jtBARTtop_varimp.png", width=500, height=500)
varimp(jtBARTtop, plot=TRUE)
}
dev.off()


#-------------------------------------------------------------------------
# And, now, predict to the full range of the Mojave 

# Unfortunately this will take some extra work, because BART model objects 
# don't quite work with SpatRaster objects, the data format used in terra

# Make a NaN vector with a value for every cell in a layer of BClim
pred_bart <- predict(jtBARTtop, stack(BClim), splitby=20, quantiles=c(0.025,0.975))
pred_bart # note that we have three layers in this 
# These are the mean and specified quantiles, the upper and lower bounds of the 95% posterior density

# Save (and reload, if later)
writeRaster(pred_bart, "output/jtBARTtop_SDM_pred.tiff", overwrite=TRUE) 
pred_bart <- rast("output/jtBARTtop_SDM_pred.tiff")

# Mask to the same "joshua tree range" we created earlier
# Create a spatial polygon defining the range from which pseudo-absences are drawn
jt_range <- read.csv("data/JT_obs.txt", sep="\t") %>% # original presence records
  st_as_sf(coords=c("lon", "lat"), crs=4326) %>% # coverted to sf, scaled in degrees
  st_transform(crs=3857) %>% # transformed to scaling in meters
  st_buffer(50000) %>% st_union() %>% # buffer by ... 10km?
  st_convex_hull() %>% # Convex hull around the resulting polygon
  st_simplify(preserveTopology=TRUE, dTolerance=5000) %>% st_buffer(10000) %>% 
  st_transform(crs=4326) %>% st_as_sf() # back to lat-lon

pred_bart.masked <- mask(pred_bart, jt_range)

# reformat as a dataframe, for figure generation
jtBARTtop.df <- cbind(crds(pred_bart.masked), as.data.frame(pred_bart.masked)) %>% 
  rename(prJT = jtBARTtop_SDM_pred_1, lo95=jtBARTtop_SDM_pred_2, up95=jtBARTtop_SDM_pred_3, lon=x, lat=y)
glimpse(jtBARTtop.df)

#-------------------------------------------------------------------------
# Generate map figures
states <- ne_states(country="united states of america", returnclass="sf")
coast <- ne_coastline(scale=10, returnclass="sf")

{png("topics/03_BART_SDMs/jtBARTtop_predicted.png", width=750, height=750)
  
  ggplot() + 
    geom_sf(data=coast, color="slategray2", linewidth=3) + 
    geom_sf(data=states, fill="cornsilk2", color="antiquewhite4") + 
    geom_tile(data=jtBARTtop.df, aes(x=lon, y=lat, fill=prJT)) +
    geom_sf(data=states, fill=NA, color="antiquewhite4") + 
    
    scale_fill_gradient(low="#efedf5", high="#756bb1", name="Pr(present), BART") +
    
    coord_sf(xlim = c(-119.5, -112), ylim = c(33, 38.5), expand = FALSE) + 
    theme_bw(base_size=18) + 
    theme(axis.title=element_blank(), panel.background=element_rect(fill="slategray3"), panel.grid.major=element_blank(), legend.position="top", legend.direction="horizontal", legend.key.width=unit(50, "points"))
  
}
dev.off()

# illustrate the mean predictions with original data overlaid
{png("topics/03_BART_SDMs/jtBARTtop_predicted_data.png", width=750, height=750)
  
ggplot() + 
    geom_sf(data=coast, color="slategray2", linewidth=3) + 
    geom_sf(data=states, fill="cornsilk2", color="antiquewhite4") + 
    geom_tile(data=jtBARTtop.df, aes(x=lon, y=lat, fill=prJT)) +
    geom_sf(data=states, fill=NA, color="antiquewhite4") + 
    geom_point(data=filter(PA, JT==1), aes(x=lon, y=lat), shape=20, color="forestgreen", size=1) + 
    geom_point(data=filter(PA, JT==0), aes(x=lon, y=lat), shape=20, color="darkorange", size=1) +
    
    scale_fill_gradient(low="#efedf5", high="#756bb1", name="Pr(present), BART") +
    
    coord_sf(xlim = c(-119.5, -112), ylim = c(33, 38.5), expand = FALSE) + 
    theme_bw(base_size=18) + 
    theme(axis.title=element_blank(), panel.background=element_rect(fill="slategray3"), panel.grid.major=element_blank(), legend.position="top", legend.direction="horizontal", legend.key.width=unit(50, "points"))
  
}
dev.off()

# illustrate the width of the 95% posterior density interval
{png("topics/03_BART_SDMs/jtBARTtop_predicted_95pctWidth.png", width=750, height=750)
  
  ggplot() + 
    geom_sf(data=coast, color="slategray2", linewidth=3) + 
    geom_sf(data=states, fill="cornsilk2", color="antiquewhite4") + 
    geom_tile(data=jtBARTtop.df, aes(x=lon, y=lat, fill=up95-lo95)) +
    geom_sf(data=states, fill=NA, color="antiquewhite4") + 
    
    scale_fill_gradient(low="#feedde", high="#a63603", name="Confidence interval width") +
    
    coord_sf(xlim = c(-119.5, -112), ylim = c(33, 38.5), expand = FALSE) + 
    theme_bw(base_size=18) + 
    theme(axis.title=element_blank(), panel.background=element_rect(fill="slategray3"), panel.grid.major=element_blank(), legend.position="top", legend.direction="horizontal", legend.key.width=unit(50, "points"))
  
}
dev.off()


#-------------------------------------------------------------------------
# Compare the BRT and BART results

# reload the BRT model
library("gbm")

# reload the BRT predictions
jtBRT.pred <- rast("output/jt_BRT_SDM_pred.tiff")

BRT.pred.masked <- mask(jtBRT.pred, jt_range) # mask to the range polygon

# compare classification accuracy for BART and BRT
predictedBART <- terra::extract(pred_bart, PA[,c("lon", "lat")])[,2]
predictedBRT <- terra::extract(jtBRT.pred, PA[,c("lon", "lat")])[,2]

# calculate AUC (higher is better)
auc(PA$JT, predictedBART)
auc(PA$JT, predictedBRT)

# determine a classification cutoff for the BRT
predBRT <- prediction(predictedBRT, PA$JT)
tssBRT <- performance(predBRT, "sens", "spec")
tss.list <- (tssBRT@x.values[[1]] + tssBRT@y.values[[1]] - 1)
tss.df <- data.frame(alpha = tssBRT@alpha.values[[1]], tss = tss.list)
threshBRT <- min(tss.df$alpha[which(tss.df$tss == max(tss.df$tss))])

# get cutoff for the BART
summary(jtBARTtop)
threshBART <- 0.56

# Confusion matrix for BRT: 0,1 is observed; FALSE,TRUE is model-predicted
table(PA$JT, predictedBRT>threshBRT) 

# Confusion matrix for BART: 0,1 is observed; FALSE,TRUE is model-predicted
table(PA$JT, predictedBART>threshBART) 

# illustrate differences
BARTvBRT <- pred_bart.masked - BRT.pred.masked # difference in Pr(present)

BARTvBRT.df <- cbind(crds(BARTvBRT), as.data.frame(BARTvBRT)) %>% rename(prDiff = jtBARTtop_SDM_pred, lon = x, lat = y)
glimpse(BARTvBRT.df)


{png("topics/03_BART_SDMs/BARTvBRT_predictions.png", width=1000, height=1000)

ggplot() + 
    geom_sf(data=coast, color="slategray2", linewidth=3) + 
    geom_sf(data=states, fill="cornsilk2", color="antiquewhite4") + 
    geom_tile(data=BARTvBRT.df, aes(x=lon, y=lat, fill=prDiff)) +
    geom_sf(data=states, fill=NA, color="antiquewhite4") + 
    
    scale_fill_gradient2(low="#e66101", mid="#efedf5", high="#5e3c99", name="Pr(present)\nBART - BRT") +
    
    coord_sf(xlim = c(-119.5, -112), ylim = c(33, 38.5), expand = FALSE) + 
    theme_bw(base_size=18) + 
    theme(axis.title=element_blank(), panel.background=element_rect(fill="slategray3"), panel.grid.major=element_blank(), legend.position="top", legend.direction="horizontal", legend.key.width=unit(50, "points"))

}
dev.off()




