#' Vignette Illustrating Statistical Corridor Modeling using Maximum Likelihood
#' To accompany "A statistical framework for predicting migration corridors"
#' Nuñez, T., Hurley, M., Graves, T. Ortega, A. Sawyer, H., Fattebert, J., Merkle, J., and Kauffman, M.   
#' July 12 2022
#' fitcorridor v1.0.0
#' Tristan Nuñez

#' Overview
#' Inputs required for modeling:
#' 1. tracks = Migration tracks to which we will fit the models, in shapefile format (.shp)
#' 2. covars = Spatial covariates in raster or raster stack format (developed with GeoTiffs, see raster package for other formats). Must be in same projection as the tracks shapefile. 
#' 3. names(covars) = Names of spatial covariates in order of raster stack, e.g., "c("SHRUB", "HF")"
#' 4. model.to.fit = Model to fit, e.g. "iso(SHRUB)+iso(HF)", where SHRUB and HF are covariates
#' 5. Miscellaneous parameters:
#'  a. mycrs = Coordinate reference system of tracks and covariates
#'  b. covar.aggregate = Factor by which to aggregate covariate raster
#'  c. startval.string = Start values for optimization
#'  d. n.inside.cores = Number of cores to use for parallel processing
#'  e. neighbors = Number of neighbors for cost-distance calculation (8 or 16)
#'  f. param.factor = Coefficient scaling value to aid in model fitting (for anisotropic [denoted using aniso() in the model formula] variables)
#' See manuscript for explanation of variable names, e.g. Lh, Lts, beta_h, beta_ts

# 1. Set your working directory to the location of the fitcorridor folder
setwd("C:/fitcorridor")

# 2. Install needed packages
# If the packages in "packages" below are already installed, you can skip this step
packages <- c("gdistance", "parallel", "rgdal",
               "ggplot2", "ggspatial", "rgeos")

# this function checks if each package is installed and if not, it installs it
# (from https://gist.github.com/stevenworthington/3178163)
ipak <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg))
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}
#call the ipak function to install and / or load packages
ipak(packages)

# 3. Source the code used for model fitting. sourcing this code will load the needed packages
source("./fitcorridor_functions.R")

# 4. Read in tracks
# HF_sim_tracks.shp consists of a set of tracks simulated from a human footprint surface with 
# a beta of -1
tracks<-readOGR("./data/tracks/HF_sim_tracks.shp", stringsAsFactors = F)
plot(tracks)
# set the crs (must be the same for tracks and covariates);
# note that setting the crs is not the same as changing the crs
mycrs <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0"
crs(tracks) <- mycrs
# NOTE: tracks must have a column called "trackid" with a unique identifier for each track, e.g."deer101_fa19"
names(tracks) # check that "trackid" is a column in tracks
head(tracks)
# 5. Read in spatial covariates
# read in covariates from a tif or other format (see the help for the brick function for more information)
# in this case, we're reading in the Human Footprint layer from Leu et al. 2008, cropped to the Boise, ID area.
# Higher values represent more human influence
covars <- brick("./data/covariates/HF.tif")
# specify the Coordinate Reference System
crs(covars) <- mycrs 
# assign names to layers
names(covars) <- c("HF") # if more than one layer, use this format: c("HF", "DEM")
# plot to check 
plot(covars)
plot(tracks, add=T, pch=".")

#6. Fit the model using the FitSingleModel function
FitSingleModel(
  tracks = tracks, # must have column called trackid, same crs as covars
  covars = covars, # same crs as tracks, names(covars) corresponds to variable names in model
  model.to.fit = "iso(HF)",  # no spaces; use + to separate terms, e.g. "iso(HF)+aniso(DEM)"
  neighbors = 8,  # 8 or 16; 8 recommended at coarser resolutions
  covar.aggregate = 1,  # factor by which to aggregate covariate rasters. run ?raster::aggregate for details
  covar.standardize = c(1),  # 1 to scale variable (dividing by standard deviation), 0 for no scaling. same length as number of terms in model. Generally 1 for iso() variables and 0 for aniso() variables. 
  param.factor = c(1),  # 1 for no change. Generally 1 for iso() variables. For aniso() variables, try 10 or 100 for elevation (in meters), and 100 or 1000 for variables with Julian day units, such as snow-off date
  startval.string = c(-5, 0),  # the first term is the start value for log(lambda), typically a value between -5 and -10 works well; the others correspond to the model terms (in order), and 0 (no effect) generally works well
  dataset = "SimHerd1",  # name of track dataset; can be altered as needed
  out.folder.path = "./fittedmods/",  # place where outputs are written. 
  n.inside.cores = 4, # number of cores for parallel processing
  Lt.tiff = T # T/F if determines if Lt and Lh and a plot of Lh will be written out 
)

#' 7. Check the "./fittedmods/SimHerd1/iso(HF)" folder to explore the results! See fitcorridor_rmarkdown_vignette.html for a description of the outputs. 

