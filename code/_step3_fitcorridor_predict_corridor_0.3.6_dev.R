# fitcorridor demo code step 3: predict corridor between a new set of endpoints
# 29 March 2023
# Tristan Nuñez
# tnunez@uw.edu
# fitcorridor version 0.3.5
# built with R 4.2.3

library(ggplot2)
library(parallel)
library(gdistance)
library(rgdal)
library(pROC)
source("./code/fitcorridor_functions_0.3.5.R")

# read in the dataframe with the fitted coefficients 
allbetas <- read.csv("./test1/iso(HF)/beta_by_groups.csv")

# let's take the average of the estimates for each covariate. 
avgbeta <- aggregate(est ~ vars, data=allbetas, FUN = "mean")
# then arrange so that the estimate for loglambda comes first, followed by the variable estimates in the order they appear in the model formula

fittedcoefs <- c(avgbeta$est[avgbeta$vars=="loglambda"],
                 avgbeta$est[avgbeta$vars=="iso(HF)"])
fittedcoefs

# read in the covariates
covars <- brick("./data/covariates/HF.tif")
# specify the Coordinate Reference System
mycrs <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0"
crs(covars) <- mycrs 
# assign names to layers
names(covars) <- c("HF") # if more than one layer, use this format: c("HF", "DEM")


# then we make a pairs of points between which we will model a corridor
# we'll do two separate pairs of points.
#first pair of points
long <- c(-115, -116, -117, -115)
lat <- c(42, 43.5, 42.5, 44)
id <- c("track1", "track1", "track2", "track2")
lonlat <-data.frame(cbind(long,lat))
newpoints <- SpatialPointsDataFrame(lonlat, 
                                    data = data.frame(id),
                                    proj4string = CRS("EPSG:4326"))

# generate a list where each element is a SpatialPointsDataFrame belonging to the same id
newpoints.list <- split(newpoints, newpoints$id)
newpoints.list

predicted.surface <- PredictCorr(points.list = newpoints.list,
                        fittedbetas = fittedcoefs,
                             model.string = "iso(HF)",
                             covars = covars,
                             directions = 8)
#This code may throw the following warning; ok to ignore:
# as(<dsCMatrix>, "dgTMatrix") is deprecated since Matrix 1.5-0; do as(as(.,  "generalMatrix"), "TsparseMatrix") instead
# Appears connected to gdistance package reliance on the matrix package and lack of updating of the gdistance package. 


# output of PredictCorr is a RasterStack of the predicted corridors
predicted.surface

plot(predicted.surface[[1]]) # plot predicted surface
points(newpoints.list[[1]], col="red") # plot endpoints
