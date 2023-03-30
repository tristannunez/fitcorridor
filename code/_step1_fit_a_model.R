# fitcorridor demo code step 1: fit a model
# 30 March 2023
# Tristan Nuñez
# tnunez@uw.edu
# fitcorridor version 0.3.6
# built with R 4.2.3
# key package versions: gdistance_1.6, raster_3.6-20, terra_1.7-18

# General Notes: This workflow relies on the gdistance package, which relies on the raster package, which has recently been replaced by the terra package. Until the gdistance package is updated to use terra instead of raster, there are likely to be bugs that arise from the reliance on an outdated package and the ongoing development of the terra package. If there are errors in running the demo code, try updating the packages (especially gdistance, raster, and terra) to their latest versions (and R as well, if your version is more than several months old). If you are running the latest versions of these packages and encounter an error, please let the package author know.

# install following packages if not already installed
# reinstalling may be useful in resolving errors if the code does not run successfully
# install.packages("gdistance")
# install.packages("gridExtra")
# install.packages("pROC")
# install.packages("raster")
# install.packages("terra")
# install.packages("sf")
# install.packages("snow")
# install.packages("reshape2")

# 1. source required packages (install if not already installed)
library(gdistance)
library(rgdal)
library(parallel)
library(ggplot2)
library(gridExtra)
library(ggspatial)
library(rgeos)
library(pROC)

# 2. Source the code used for model fitting

source("./code/fitcorridor_functions.R")

# 3. Read in tracks
# HF_sim_tracks.shp consists of a set of tracks simulated from a human footprint surface
# The AID column is the animal ID
# The trackid is unique to each movement track segment; the corridor will be modeled between the first and the last location in the movement track segment. 

tracks<-readOGR("./data/tracks/HF_sim_tracks.shp", stringsAsFactors = F)
plot(tracks)
head(tracks)

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
covars <- stack("./data/covariates/HF.tif")

# If you wish to speed up run time on this demo code, you can aggregate the covariates to a coarser spatial resolution. Do this by uncommenting and then running the following line.
# covars <- stack(aggregate(covars, fact=4))

# specify the Coordinate Reference System
crs(covars) <- mycrs 
# assign names to layers
names(covars) <- c("HF") # if more than one layer, use this format: c("HF", "DEM")
# plot to check 
plot(covars)
class(covars)
plot(tracks, add=T, pch=".")

# a "groupid" column must be specified. The cost distance model will be optimized using the likelihoods from all the tracks with the same groupid. In this case, the unique AID (animal ID) is used as the group. 
tracks$groupid <- tracks$AID

FitModels(tracks.spdf = tracks,
          covars = covars,
          model.string = "iso(HF)", # to add terms, use "+" with no spaces (e.g., "iso(HF)+iso(DEM)"). The iso() specifies that cost distances are calculated isotropically. 
          directions = 8, # 8 or 16
          startvals = c(-15, 0),# the first term is the start value for log(lambda); the remaining terms are start values for the coefficients of the model terms (in order). Generally, use values between -8 and -15 for loglambda (experimentation required), and 0 for beta coefficients.
          dataset = "./test1", # can be anything, but needs to start with "./", this folder is where outputs are written in the working directory
          n.inside.cores = 2,
          n.outside.cores = 2,
          Lt.tiff = TRUE, # faster if False
          plotindivtracks=TRUE, # faster if False
          cropfactor = 3) # cropfactor determines how large an area to calculate the cost surfaces around each individual track, as a multiple of the geographic extent of the track

#6.  Check the "./test1/iso(HF)" folder to explore the results!

# 7. development environment for troubleshooting
# this code was developed with the following package versions; if issues arise, use sessionInfo() to compare the versions here with the versions on your machine. If your environment uses a newer package version and throws an error, please let Tristan know at tnunez@uw.edu

# sessionInfo() results for development environment:
# other attached packages:
#   [1] pROC_1.18.0     rgeos_0.6-2     ggspatial_1.1.7 gridExtra_2.3  
# [5] ggplot2_3.3.6   rgdal_1.6-5     gdistance_1.6   Matrix_1.5-3   
# [9] igraph_1.3.2    raster_3.6-20   sp_1.6-0       
# 
# loaded via a namespace (and not attached):
#   [1] Rcpp_1.0.10        plyr_1.8.7         pillar_1.8.0      
# [4] compiler_4.2.3     class_7.3-21       tools_4.2.3       
# [7] lifecycle_1.0.1    tibble_3.1.8       gtable_0.3.0      
# [10] lattice_0.20-45    pkgconfig_2.0.3    rlang_1.0.4       
# [13] cli_3.3.0          DBI_1.1.3          rstudioapi_0.13   
# [16] e1071_1.7-11       terra_1.7-18       withr_2.5.0       
# [19] dplyr_1.0.9        generics_0.1.3     vctrs_0.4.1       
# [22] classInt_0.4-7     grid_4.2.3         tidyselect_1.1.2  
# [25] glue_1.6.2         sf_1.0-8           R6_2.5.1          
# [28] fansi_1.0.3        purrr_0.3.4        magrittr_2.0.3    
# [31] units_0.8-0        scales_1.2.0       codetools_0.2-19  
# [34] assertthat_0.2.1   colorspace_2.0-3   KernSmooth_2.23-20
# [37] utf8_1.2.2         proxy_0.4-27       munsell_0.5.0    

# end 7. development environment section