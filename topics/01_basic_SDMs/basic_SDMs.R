# Species distribution modeling with boosted regression trees in R
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

library("dismo")
library("gbm")


set.seed(20250915)

#-----------------------------------------------------------
# Species occurrence data

# Joshua tree occurrence records, from published data: https://doi.org/10.5061/dryad.6s67t
# ... supplemented with validated iNaturalist records and field notes, and thinned a bit.

jtOcc <- read.csv("data/JT_obs.txt", sep="\t")

glimpse(jtOcc)

# Visualize the points on a map of the US Southwest
# As you can see, there's quite a bit of heterogeneity in sampling density.
states <- ne_states(country="united states of america", returnclass="sf")
coast <- ne_coastline(scale=10, returnclass="sf")

{png("topics/01_basic_SDMs/JT_occurrences_map.png", height=750, width=750)
ggplot() + 
  geom_sf(data=coast, color="slategray2", linewidth=3) + 
  geom_sf(data=states, fill="cornsilk2", color="antiquewhite4") + 
  geom_point(data=jtOcc, aes(x=lon, y=lat), shape=21, color="forestgreen", size=2) + 
  coord_sf(xlim = c(-119.5, -112), ylim = c(33, 38.5), expand = FALSE) + 
  theme_bw(base_size=18) + theme(axis.title=element_blank(), panel.background=element_rect(fill="slategray3"), panel.grid.major=element_blank(), legend.position="none")
}
dev.off()


#-------------------------------------------------------------------------
# Climate data

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
# Create pseudo-absence records

# Create a spatial polygon defining the range from which pseudo-absences are drawn
jt_range <- jtOcc %>% # original presence records
  st_as_sf(coords=c("lon", "lat"), crs=4326) %>% # coverted to sf, scaled in degrees
  st_transform(crs=3857) %>% # transformed to scaling in meters
  st_buffer(50000) %>% st_union() %>% # buffer by ... 10km?
  st_convex_hull() %>% # Convex hull around the resulting polygon
  st_simplify(preserveTopology=TRUE, dTolerance=5000) %>% st_buffer(10000) %>% 
  st_transform(crs=4326) %>% st_as_sf() # back to lat-lon

# visualize this range polygon with the original data
{png("topics/01_basic_SDMs/JT_occurrences_range_map.png", height=750, width=750)
ggplot() + 
  geom_sf(data=coast, color="slategray2", linewidth=3) + 
  geom_sf(data=states, fill="cornsilk2", color="antiquewhite4") + 
  geom_sf(data=jt_range, fill="cornsilk3", color=NA, alpha=0.5) +
  geom_point(data=jtOcc, aes(x=lon, y=lat), shape=21, color="forestgreen", size=2) + 
  coord_sf(xlim = c(-119.5, -112), ylim = c(33, 38.5), expand = FALSE) + 
  theme_bw(base_size=18) + theme(axis.title=element_blank(), panel.background=element_rect(fill="slategray3"), panel.grid.major=element_blank(), legend.position="none")
}
dev.off()

# mask a BIOCLIM layer to the range polygon, take the cell coordinates
pseudabs <- terra::mask(BClim[[1]], vect(st_transform(jt_range, crs=crs(BClim[[1]])))) %>% 
  crds(df=TRUE) %>% rename(lon=x, lat=y) %>% 
  filter(!(paste(lon,lat) %in% paste(jtOcc$lon, jtOcc$lat))) %>% # remove cells represented in presence
  slice_sample(n=nrow(jtOcc)) # pull a sample the same size as presence

glimpse(pseudabs) # should see the same number of rows as jtOcc

# How do the pseudoabsence points look?
{png("topics/01_basic_SDMs/JT_occurrences_pseudoabsences_range_map.png", height=750, width=750)
  ggplot() + 
    geom_sf(data=coast, color="slategray2", linewidth=3) + 
    geom_sf(data=states, fill="cornsilk2", color="antiquewhite4") + 
    geom_sf(data=jt_range, fill="cornsilk3", color=NA, alpha=0.5) +
    geom_point(data=jtOcc, aes(x=lon, y=lat), shape=21, color="forestgreen", size=2) + 
    geom_point(data=pseudabs, aes(x=lon, y=lat), shape=20, color="darkorange", size=1) +
    coord_sf(xlim = c(-119.5, -112), ylim = c(33, 38.5), expand = FALSE) + 
    theme_bw(base_size=18) + theme(axis.title=element_blank(), panel.background=element_rect(fill="slategray3"), panel.grid.major=element_blank(), legend.position="none")
}
dev.off()


#-------------------------------------------------------------------------
# Assemble environmental values for our presences and pseudoabsences 

pres_envs <- data.frame(terra::extract(BClim, jtOcc[,c("lon","lat")]))
glimpse(pres_envs) # BioClim values for the presence locations

abs_envs <- data.frame(terra::extract(BClim, pseudabs[,c("lon","lat")]))
glimpse(abs_envs) # And ditto for the pseudoabsences

# Assemble a single dataframe with presence/absence coordinates and the
# associated BioClim values:
PA <- rbind(data.frame(lon=jtOcc$lon, lat=jtOcc$lat, JT=1, pres_envs), 
          data.frame(lon=pseudabs$lon, lat=pseudabs$lat, JT=0, abs_envs)) %>% 
          mutate(ID = row_number()) %>% filter(!is.na(AMT))

glimpse(PA) # should be ~twice the size of the presence data set

# You may want to write out this data frame for later /read back in 
if(!dir.exists("output")) dir.create("output")

write.table(PA, "output/JT_presence-pseudoabsence-envs.txt", sep=",", col.names=TRUE, row.names=FALSE, quote=FALSE)
PA <- read.csv("output/JT_presence-pseudoabsence-envs.txt")
glimpse(PA)

#-------------------------------------------------------------------------
# Fit a species distribution model with Boosted Regression Trees

# vector of the Bioclim variables
xnames <- colnames(PA)[5:23]

# Fit a BRT model with defaults provided by `dismo`
jtBRT <- gbm.step(data=PA, gbm.x=xnames, gbm.y="JT", family = "bernoulli", tree.complexity = 5, learning.rate = 0.05, bag.fraction = 0.5)

jtBRT

# Examine predictor contributions and write out a figure
{png("topics/01_basic_SDMs/jtBRT_predictor_contributions.png", width=500, height=1500)
summary(jtBRT)
}
dev.off()

# Try to simplify the model from the full set of 19 Bioclim variables
jtBRT.simp <- gbm.simplify(jtBRT)

{png("topics/01_basic_SDMs/jtBRT.simp_predictor_contributions.png", width=500, height=750)
summary(jtBRT.simp)
}
dev.off()

# Fit a new stepwise model using the optimal predictors identified (in my run, this is 9, your mileage may vary)
jtBRT2 <- gbm.step(data=PA, gbm.x=jtBRT.simp$pred.list[[9]], gbm.y="JT", family = "bernoulli", tree.complexity = 5, learning.rate = 0.05, bag.fraction = 0.5)

write_rds(jtBRT2, file="output/models/jtBRT.rds") # save the results
jtBRT2 <- read_rds("output/models/jtBRT.rds") # reload later


{png("topics/01_basic_SDMs/jtBRT2_predictor_contributions.png", width=500, height=750)
summary(jtBRT2)
}
dev.off()


# Examine predictor partial effects
{png("topics/01_basic_SDMs/jtBRT_partials.png", width=1000, height=750)
gbm.plot(jtBRT2, n.plots=10, plot.layout=c(3,4), write.title=FALSE)
}
dev.off()


# And, finally, predict to the full range of the Mojave:
jtBRT.pred <- predict(BClim, jtBRT2, n.trees=jtBRT2$gbm.call$best.trees, type="response")

writeRaster(jtBRT.pred, "output/jt_BRT_SDM_pred.tiff", overwrite=TRUE) 
jtBRT.pred <- rast("output/jt_BRT_SDM_pred.tiff")

jtBRT.pred.masked <- mask(jtBRT.pred, jt_range)

# reformat as a dataframe, for figure generation
jtBRT.pred.df <- cbind(crds(jtBRT.pred.masked), as.data.frame(jtBRT.pred.masked)) %>% rename(prJT = lyr1, lon=x, lat=y)
glimpse(jtBRT.pred.df)

{png("topics/01_basic_SDMs/jtBRT_predicted.png", width=750, height=750)

ggplot() + 
  geom_sf(data=coast, color="slategray2", linewidth=3) + 
  geom_sf(data=states, fill="cornsilk2", color="antiquewhite4") + 
  geom_tile(data=jtBRT.pred.df, aes(x=lon, y=lat, fill=prJT)) +
  geom_sf(data=states, fill=NA, color="antiquewhite4") + 

  scale_fill_gradient(low="#efedf5", high="#756bb1", name="Pr(present), BRT") +
  
  coord_sf(xlim = c(-119.5, -112), ylim = c(33, 38.5), expand = FALSE) + 
  theme_bw(base_size=18) + 
  theme(axis.title=element_blank(), panel.background=element_rect(fill="slategray3"), panel.grid.major=element_blank(), legend.position="top", legend.direction="horizontal", legend.key.width=unit(50, "points"))

}
dev.off()

{png("topics/01_basic_SDMs/jtBRT_predicted_data.png", width=750, height=750)
  
  ggplot() + 
    geom_sf(data=coast, color="slategray2", linewidth=3) + 
    geom_sf(data=states, fill="cornsilk2", color="antiquewhite4") + 
    geom_tile(data=jtBRT.pred.df, aes(x=lon, y=lat, fill=prJT)) +
    geom_sf(data=states, fill=NA, color="antiquewhite4") + 
    geom_point(data=filter(PA, JT==1), aes(x=lon, y=lat), shape=20, color="forestgreen", size=1) + 
    geom_point(data=filter(PA, JT==0), aes(x=lon, y=lat), shape=20, color="darkorange", size=1) +
    
    scale_fill_gradient(low="#efedf5", high="#756bb1", name="Pr(present), BRT") +
    
    coord_sf(xlim = c(-119.5, -112), ylim = c(33, 38.5), expand = FALSE) + 
    theme_bw(base_size=18) + 
    theme(axis.title=element_blank(), panel.background=element_rect(fill="slategray3"), panel.grid.major=element_blank(), legend.position="top", legend.direction="horizontal", legend.key.width=unit(50, "points"))
  
}
dev.off()







