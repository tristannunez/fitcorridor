# fitcorridor demo code step 2: crossvalidate results
# 30 March 2023
# Tristan Nuñez
# tnunez@uw.edu
# fitcorridor version 0.3.6
# built with R 4.2.3

# NOTE: This script only works with relatively small landscapes due to way the calc_ROC() function (see fitcorridor_functions.R) is currently written, because it reads in the entire landscape into memory. Will require updating for landscapes with large numbers of cells. 
# ALSO NOTE: As currently written, this code calculates the ROC and AUC metrics on the training data. The code will need to be modified for cases where separate test data are used to calculate ROC and AUC metrics by changing the object "tracks.spdf" to consist of test data, and specifying the desired beta coefficients for the test by changing the object "focalbetadf". 

# this script calculates the ROC and AUC for null and fitted models, and writes out a plot comparing the fitted and null models, as well as a .csv that can be later used to compare among different models
library(ggplot2)
library(parallel)
library(gdistance)
library(rgdal)
library(pROC)
source("./code/fitcorridor_functions.R")

# read in the dataframe with the fitted coefficients 
#allbetas <- read.csv("./socal_300m/iso(slope)/beta_ts.csv")
allbetas <- read.csv("./test1/iso(HF)/beta_by_groups.csv")
# need to set nvars, the number of variables in the model
nvars <- 1 # for iso(HF), nvars = 1; for iso(slope)+iso(shrubs), nvars=2, and so on

# read in the tracks for validation (in this case it's the same as the ones for training)
tracks<-readOGR("./data/tracks/HF_sim_tracks.shp", stringsAsFactors = F)
mycrs <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0"
crs(tracks) <- mycrs
plot(tracks)
head(tracks)

tracks.spdf <- tracks
tracks.spdf$groupid <- tracks$AID
track.list.by.groupid <- split(tracks.spdf, tracks.spdf$groupid) # make a list of tracks

covars <- brick("./data/covariates/HF.tif")
# specify the Coordinate Reference System
crs(covars) <- mycrs 
# assign names to layers
names(covars) <- c("HF") # if more than one layer, use this format: c("HF", "DEM")

# this for loop goes through each animal to calculate the proportion of landscape needed to capture each proportion of points for the fitted model as well as for a null model; generates a dataframe that is then written out.
focalgroupid <-unique(allbetas$groupid)[1]
for(focalgroupid in unique(allbetas$groupid)){
  print(focalgroupid)
  focaltracks <- tracks.spdf[tracks.spdf$groupid==focalgroupid,]
  focalbetadf <- allbetas[allbetas$groupid==focalgroupid,]
  
  # first calculate the ROC object for fitted models
  roc_fitted <- XVal(tracks = focaltracks,
                     fittedbetas = focalbetadf$est,
                     model.string = unique(focalbetadf$model)[1],
                     covars = covars,
                     cropfactor = 3,
                     ncores = 4, # uses just one set of cores; one core for each track 
                     directions = 8)
  roc_fitted$groupid <- focalgroupid
  roc_fitted$model <- unique(focalbetadf$model)[1]
  
  roc.df.fitted <- data.frame(sensitivity = roc_fitted$sensitivities,
                              specificity = roc_fitted$specificities,
                              model = unique(focalbetadf$model)[1],
                              groupid = focalgroupid,
                              type = "fitted")
  
  # then calculate the ROC object for a null model using the same lambda as the fitted model 
  roc_null <- XVal(tracks = focaltracks,
                   fittedbetas = c(focalbetadf$est[1], rep(0,nvars)),
                   model.string = unique(focalbetadf$model)[1],
                   covars = covars,
                   cropfactor = 3,
                   ncores = 4, 
                   directions = 8)
  roc_null$groupid <- focalgroupid
  roc_null$model <- "NULL"
  
  roc.df.null <- data.frame(sensitivity = roc_null$sensitivities,
                            specificity = roc_null$specificities,
                            model = "NULL",
                            groupid = focalgroupid,
                            type = "NULL")
  
  auc.df <- data.frame(auc = c(roc_fitted$auc, roc_null$auc),
                       type = c("fitted", "null"),
                       model = c(unique(focalbetadf$model)[1], "null"),
                       groupid = c(focalgroupid, focalgroupid))
  
  if(focalgroupid==unique(allbetas$groupid)[1]){
    out.df <- rbind(roc.df.fitted, roc.df.null)
    auc.df.out <- auc.df
    
  } else {
    out.df <- rbind(out.df, rbind(roc.df.fitted, roc.df.null))
    auc.df.out <- rbind(auc.df.out, auc.df)
  } 
}  

#This code may throw the following warning; ok to ignore:
# as(<dsCMatrix>, "dgTMatrix") is deprecated since Matrix 1.5-0; do as(as(.,  "generalMatrix"), "TsparseMatrix") instead
# Appears connected to gdistance package reliance on the matrix package and lack of updating of the gdistance package. 

write.csv(out.df, "./test1/iso(HF)/ROC.csv")
write.csv(auc.df.out, "./test1/iso(HF)/AUC.csv")

ROC.plot <-ggplot(out.df)+
  geom_line(aes(x=specificity, y=sensitivity, color=model))+
  facet_wrap(.~groupid)+
  scale_x_reverse()+
  theme_bw()

AUC.boxplot <- ggplot(auc.df.out)+
  geom_boxplot(aes(y=auc, x=model, color=model))+
  theme_bw()

pdf("./test1/iso(HF)/rocplots.pdf",
    width=11, height=8.5, onefile=T)
ROC.plot
AUC.boxplot
dev.off()

# TODO: Calculate ROC and AUC for all groupids combined



