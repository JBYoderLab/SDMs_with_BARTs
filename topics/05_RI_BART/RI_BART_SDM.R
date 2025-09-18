# Species distribution modeling with RI BARTs in R
# Jeremy B. Yoder, 18 Sept 2024

# Clears the environment and load key packages
rm(list=ls())

# This script assumes your working directory is the project directory; you may want to change this
# setwd("~/Documents/Active_projects/SDMs_with_BARTs")

library("tidyverse") 
library("sf")
# library("terra")
library("raster")

library("rnaturalearth")
library("rnaturalearthdata")
library("geodata")

library("embarcadero")
# devtools::install_github("cjcarlson/embarcadero", build_vignettes = TRUE, force = TRUE)


set.seed(20250915)

#-----------------------------------------------------------
# Load and prepare data

# Joshua tree occurrence records, and pseudoabsences, with environmental data
# (This was generated in the code for topic 01, `BART_SDM.R`)
PA <- read.csv("output/JT_presence-pseudoabsence-envs.txt")
glimpse(PA)

# Environmental data and accessory stuff
MojExt <- c(-119.5, -112, 33, 38.5) # Extent of the Mojave desert

# climate (Bioclim, previously downloaded at max res)
BClim <- crop(worldclim_tile("bio", -120, 33, path="../data/Yucca/"), MojExt)
BClim

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


# Visualize the occurrence points on a map of the US Southwest
# As you will recall, there's quite a bit of heterogeneity in sampling density.
states <- ne_states(country="united states of america", returnclass="sf")
coast <- ne_coastline(scale=10, returnclass="sf")

{png("topics/05_RI_BART/JT_occurrences_map.png", height=750, width=750)
  ggplot() + 
    geom_sf(data=coast, color="slategray2", linewidth=3) + 
    geom_sf(data=states, fill="cornsilk2", color="antiquewhite4") + 
    geom_point(data=filter(PA, JT==1), aes(x=lon, y=lat), shape=21, color="forestgreen", size=2) + 
    coord_sf(xlim = c(-119.5, -112), ylim = c(33, 38.5), expand = FALSE) + 
    theme_bw(base_size=18) + theme(axis.title=element_blank(), panel.background=element_rect(fill="slategray3"), panel.grid.major=element_blank(), legend.position="none")
}
dev.off()

# The records are a LOT denser in Nevada than they are the other states where
# Joshua trees grow. That reflects different data sources --- Nevada is covered
# by more detailed surveys. 

# Does that variation in sampling density mean a SDM trained on this data could
# be confounding Nevada-like environments with environments that actually predict
# Joshua tree presence?

# A random-intercept model will let us answer this.

#-------------------------------------------------------------------------
# First, generate a new environmental layer for state IDs

MojStates <- raster::rasterize(filter(states, name_en %in% c("California", "Nevada", "Utah", "Arizona")), BClim[[1]], "name_en") 
plot(MojStates) # confirm we've got a raster layer with state IDs

# but we need numeric values to provide as predictors for BARTs, so
MojStates.numeric <- MojStates %>% as.numeric() 
plot(MojStates.numeric) # there we go

# add state as a layer to our RasterStack of predictors
envs <- raster::stack(c(BClim, MojStates.numeric))
names(envs)[[nlayers(envs)]] <- "State"
plot(envs[["State"]])

# add state as a column to our data frame of presence and pseudoabsence locations
PA$State <- raster::extract(envs[["State"]], PA[,c("lon", "lat")])

glimpse(PA) # should now have a new numeric-value column


#-------------------------------------------------------------------------
# Fit a species distribution model with random-intercept BARTs

# Adding a random-intercept effect to a BART model lets us account for 
# heterogeneity among data sources. Recall that our occurrence records
# are much denser in Nevada --- we can treat state as an RI predictor to
# account for this.

# We can simply add the RI predictor into a model with predictors found by
# previous predictor selection ...
topX <- c("PS", "PDQ", "PWaQ", "MTCM", "PDM", "MTWeQ", "TS")

# Then fit the model, specifying State as an RI predictor ("group.by")
jtRIBART <- rbart_vi(as.formula(paste(paste('JT', paste(topX, collapse=' + '), sep = ' ~ '), 'State', sep=' - ')), data = PA, group.by = PA[,'State'], n.chains = 1, k = 2, power = 2, base = 0.95, keepTrees = TRUE)

# as before
invisible(jtRIBART$fit[[1]]$state)
write_rds(jtRIBART, file="output/models/jtRIBART.rds") 
jtRIBART <- read_rds(file="output/models/jtRIBART.rds")

{png("topics/05_RI_BART/JT_RI_BART_summary.png", height=750, width=750)
  
summary(jtRIBART)

}
dev.off()

# We can inspect the estimated RI effect to see whether it makes sense to include:

{png("topics/05_RI_BART/BART_RI_estimates.png", width=500, height=350)

plot.ri(jtRIBART, temporal=FALSE) + 
  scale_x_discrete(labels=c("Arizona", "California", "Nevada", "Utah")) + # relabel numeric coding
  labs(title="Random intercept effects of state", x="RI term") + 
  theme_bw(base_size=18)

}
dev.off()

# The RI effect is essentially a different baseline probability of Joshua tree
# presence for each state, and this shows that NV does have a higher baseline!

#-------------------------------------------------------------------------
# Make predictions while controlling for the RI term

# We can do this directly because the inputs are formatted as a RasterStack
pred_RI <- predict(jtRIBART, envs[[topX]], splitby=20, ri.data=envs[["State"]], ri.name='State', ri.pred=FALSE)

# save (and reload, if later)
writeRaster(pred_RI, "output/jtRIBART_SDM_pred.tiff", overwrite=TRUE) 
pred_RI <- rast("output/jtRIBART_SDM_pred.tiff")

# Mask to the same "joshua tree range" we created earlier
# Create a spatial polygon defining the range from which pseudo-absences are drawn
jt_range <- read.csv("data/JT_obs.txt", sep="\t") %>% # original presence records
  st_as_sf(coords=c("lon", "lat"), crs=4326) %>% # coverted to sf, scaled in degrees
  st_transform(crs=3857) %>% # transformed to scaling in meters
  st_buffer(50000) %>% st_union() %>% # buffer by ... 10km?
  st_convex_hull() %>% # Convex hull around the resulting polygon
  st_simplify(preserveTopology=TRUE, dTolerance=5000) %>% st_buffer(10000) %>% 
  st_transform(crs=4326) %>% st_as_sf() # back to lat-lon

pred_RI.masked <- mask(pred_RI, jt_range)

# reformat as a dataframe, for figure generation
jtRIBART.df <- cbind(crds(pred_RI.masked), as.data.frame(pred_rast.masked)) %>% rename(prJT = jtRIBART_SDM_pred, lon=x, lat=y)
glimpse(jtRIBART.df)

#-------------------------------------------------------------------------
# Generate a map figure
states <- ne_states(country="united states of america", returnclass="sf")
coast <- ne_coastline(scale=10, returnclass="sf")

{png("topics/05_RI_BART/jtRIBART_predicted.png", width=750, height=750)
  
  ggplot() + 
    geom_sf(data=coast, color="slategray2", linewidth=3) + 
    geom_sf(data=states, fill="cornsilk2", color="antiquewhite4") + 
    geom_tile(data=jtRIBART.df, aes(x=lon, y=lat, fill=prJT)) +
    geom_sf(data=states, fill=NA, color="antiquewhite4") + 
    
    scale_fill_gradient(low="#efedf5", high="#756bb1", name="Pr(present)") +
    
    coord_sf(xlim = c(-119.5, -112), ylim = c(33, 38.5), expand = FALSE) + 
    theme_bw(base_size=18) + 
    theme(axis.title=element_blank(), panel.background=element_rect(fill="slategray3"), panel.grid.major=element_blank(), legend.position="top", legend.direction="horizontal", legend.key.width=unit(50, "points"))
  
}
dev.off()


#-------------------------------------------------------------------------
# Compare the RI BART and BART results

# Reload our previously generated BART model predictions
pred_basic <- rast("output/jtBARTtop_SDM_pred.tiff")

# Calculate the difference in Pr(Present)
# between the RI model and the previous non-RI model
RIvBART <- pred_RI - pred_basic 

RIvBART.df <- cbind(crds(RIvBART), as.data.frame(RIvBART)) %>% rename(prDiff = jtRIBART_SDM_pred, lon = x, lat = y)
glimpse(RIvBART.df)

{png("topics/05_RI_BART/jtRIBART_vs_BART.png", width=750, height=750)
  
  ggplot() + 
    geom_sf(data=coast, color="slategray2", linewidth=3) + 
    geom_sf(data=states, fill="cornsilk2", color="antiquewhite4") + 
    geom_tile(data=RIvBART.df, aes(x=lon, y=lat, fill=prDiff)) +
    geom_sf(data=states, fill=NA, color="antiquewhite4") + 
    
    scale_fill_gradient2(low="#e66101", mid="#efedf5", high="#5e3c99", name="Difference in Pr(present),\nRI BART vs BART") +
    
    coord_sf(xlim = c(-119.5, -112), ylim = c(33, 38.5), expand = FALSE) + 
    theme_bw(base_size=18) + 
    theme(axis.title=element_blank(), panel.background=element_rect(fill="slategray3"), panel.grid.major=element_blank(), legend.position="top", legend.direction="horizontal", legend.key.width=unit(50, "points"))
  
}
dev.off()




