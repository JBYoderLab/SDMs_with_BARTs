# Species distribution modeling with BARTs in R
# Predictor partial effects, in multiple formats
# Jeremy B. Yoder, 1 Oct 2025

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
BClim <- crop(worldclim_tile("bio", -120, 33, path="data"), MojExt)

plot(BClim[[5]]) # one layer from the SpatRaster stack

names(BClim) <- c("AMT", "MDR", "ITH", "TS", "MWaMT", "MCMT", "TAR", "MWeQT", "MDQT", "MWaQT", "MCQT", "AP", "PWeM", "PDM", "PS", "PWeQ", "PDQ", "PWaQ", "PCQ")

plot(BClim[["AMT"]]) # Annual mean temperature
plot(BClim[["AP"]]) # Annual precipitation

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
# using previously-identified best predictors

# vector of the Bioclim variables
topX <- c("PS", "PDQ", "PWaQ", "MTCM", "PDM", "MTWeQ", "TS")

# Train a new model with those top predictors
jtBARTtop <- bart(y.train=PA[,"JT"], x.train=PA[,topX], keeptrees=TRUE)

summary(jtBARTtop)


# Estimate and visualize predictor partial effects; this will be quite slow!
jtBART.partials <- partial(jtBARTtop, topX, trace=FALSE, smooth=5) 

write_rds(jtBART.partials, file="output/models/jtBARTtop.partials.rds") # save the results
jtBART.partials <- read_rds("output/models/jtBARTtop.partials.rds") # reload later

jtBART.partials

# Build a multi-panel partials figure
partvals <- data.frame(predictor=rep(topX, each=nrow(jtBART.partials[[1]]$data)), do.call("rbind", lapply(jtBART.partials, function(x) x$data))) %>% 
  mutate(predictor=factor(predictor, topX))

{png("topics/04_BART_interpretation/jtBART_partials.png", width=1000, height=500)

ggplot(partvals) + 
  geom_ribbon(aes(x=x, ymin=q05, ymax=q95), fill="#99d8c9") + 
  geom_line(aes(x=x, y=med), color="white", linewidth=1) + 
  
  # consider adjusting row number and figure dimensions
  facet_wrap("predictor", nrow=2, scale="free") + 
  
  labs(y="Marginal Pr(Present)") + 
  theme_bw(base_size=18) + theme(axis.title.x=element_blank(), panel.spacing=unit(0.2,"in"))

}
dev.off()

# It can be instructive to compare the partials to the distributions of each
# predictor in the original training data. Let's set that up ...
training_predictors <- PA %>% dplyr::select(ID, JT, all_of(topX)) %>% 
  pivot_longer(all_of(topX), names_to="predictor", values_to="value") %>%
  mutate(predictor=factor(predictor, topX), JT=JT==1)

glimpse(training_predictors)

{png("topics/04_BART_interpretation/jtBART_training.png", width=1000, height=500)
  
ggplot(training_predictors, aes(x=value, fill=JT, group=JT)) + 
  geom_histogram(position="dodge") +
  
  # consider adjusting row number and figure dimensions
  facet_wrap("predictor", nrow=2, scale="free") + 
  
  scale_fill_manual(values=c("#dfc27d","#018571"), "JT present?") + 
  
  labs(y="Grid cells") + 
  theme_bw(base_size=18) + 
  theme(axis.title.x=element_blank(), panel.spacing=unit(0.2,"in"), legend.position="inside", legend.position.inside=c(0.9,0.3))
  
}
dev.off()


#-------------------------------------------------------------------------
# Spatial partial effects

# Another useful illustration is to project partial effects in space
# We call these "spatial partials" or "spartials"

# This will also be quite slow!
# You can also break up the task by setting x.vars=[one predictor]
jtSpartials <- spartial(jtBARTtop, stack(BClim), x.vars=topX)
names(jtSpartials) <- topX

writeRaster(jtSpartials, "output/jtBARTtop_spartials.tiff", overwrite=TRUE)
jtSpartials <- rast("output/jtBARTtop_spartials.tiff")
names(jtSpartials) <- topX

# reformat for visualization
jtSpartials.df <- cbind(crds(jtSpartials), as.data.frame(jtSpartials)) %>% 
                  rename(lon = x, lat = y) %>%
                  pivot_longer(all_of(topX), names_to="Predictor", values_to="Marginal_pr")

glimpse(jtSpartials.df)

{png("topics/04_BART_interpretation/jtBART_spartials.png", width=1000, height=500)

ggplot() + 
  geom_tile(data=jtSpartials.df, aes(x=lon, y=lat, fill=Marginal_pr)) +
  
  # consider adjusting row number and figure dimensions
  facet_wrap("Predictor", nrow=2, scale="free") + 
    
  scale_fill_gradient(low="#f6eff7", high="#016c59", name="Marginal\nPr(present)") +
  
  theme_bw(base_size=18) + 
  theme(axis.title=element_blank(), panel.spacing=unit(0.2,"in"), legend.position="inside", legend.position.inside=c(0.9,0.3))
  
}
dev.off()




