# fitcorridor functions
# 29 March 2023
# Tristan Nuñez
# tnunez@uw.edu
# fitcorridor version 0.3.6
# built with R 4.2.3
# key package versions: gdistance_1.6, raster_3.6-20, terra_1.7-18


#' TODO: 
#' 1. write function that translates output from optim into summary statistics.
#' 2. have objective function calculate NLL from a null model
#' 3. write null NLL to output dataframe
#' 

FitModels <- function(tracks.spdf,
                      covars,
                      model.string = "iso(HF)", #"iso(SHRUB)+iso(HF)" # no spaces
                      directions = 8, # 8 or 16
                      startvals = c(-15, 0),# NULL #c(-15 0 0) # the first term is the start value for log(lambda); the others correspond to the model terms (in order)
                      dataset = "./defaultherd", # name of output folder 
                      n.inside.cores = NULL,
                      n.outside.cores = NULL,
                      Lt.tiff = T ,
                      plotindivtracks = FALSE,
                      cropfactor = 3 #,
                      #codefile = "./code/code/corridor_fitting_functions_fitmultiple_1.0.R"
){

# Some simple error checking
  stopifnot("tracks.spdf is not a SpatialPointsDataFrame" = class(tracks.spdf) == "SpatialPointsDataFrame")
  stopifnot("no column labeled 'trackid' in tracks.spdf" = "trackid" %in% names(tracks.spdf))
  stopifnot("no column labeled 'groupid' in tracks.spdf" = "groupid" %in% names(tracks.spdf))
  stopifnot("CRS of 'tracks.spdf' and CRS of 'covars' not identical" = compareCRS(crs(tracks.spdf), crs(covars)))
  stopifnot("covars is not a RasterStack" = class(covars)=="RasterStack" | class(covars)=="RasterBrick" | class(covars)=="Raster")
  stopifnot("'startvals' needs to have one more element than the number of variables in 'model.string'" = length(startvals) == length(strsplit(model.string, "\\+")[[1]]) + 1)
  stopifnot("directions must be either 4, 8, or 16" = directions == 4 | directions == 8 | directions == 16)

  
  # next 4 lines pull out info from model.string 
  model.split1 <- unlist(strsplit(as.character(model.string), split = "\\+"))
  ms2 <- sapply(strsplit(model.split1, "\\(") ,"[[", 2)
  variables <- sapply(strsplit(ms2, "\\)"),"[[", 1)

# to test    
#  stopifnot("variable in model not one of named layers in covars" = variables %in% names(covars))
  
  nvars <- length(variables)
  npars <- nvars+1 # add 1 for lambda
  # set up out files
  dataset <- as.character(dataset)
  if(!dir.exists(paste(#out.folder.path, 
    dataset, sep=""))){
    dir.create(paste(#out.folder.path, 
      dataset, sep=""))  
  }
  foldername <- model.string
  folderpath <- paste(#out.folder.path, 
    dataset, "/", foldername, sep="") # might need to tweak
  dir.create(folderpath)
  params_nll_out_file_optim <- paste(folderpath, "/params_nll_out_optim.csv", sep="")
  
#  covars <- scale(covars, center=F)
  
  #covars <- aggregate(covars, fact=covar.aggregate)
  
  model.list.out <- MakeModelList(model.string = as.character(model.string),
                                  covars = covars
                                  #covar.standardize = covar.standardize
                                  )
  focal.model.list <- model.list.out$focal.model.list

  track.list.by.groupid <- split(tracks.spdf, tracks.spdf$groupid)
  
  model.list.by.track.by.groupid <- list()
  for(a in 1:length(track.list.by.groupid)){
    print(a)
    tracks <- track.list.by.groupid[[a]]
    track.list <- split(tracks, tracks$trackid)
    model.list.by.track <- lapply(X=track.list, FUN=MakeModelList2,
                                model.string = as.character(model.string),
                                covars = covars, 
                                cropfactor = cropfactor
                                #covar.standardize = covar.standardize
                                )
    model.list.by.track.by.groupid[[a]]<-model.list.by.track
  }

gc()  

  cl.wrap <- parallel::makeCluster(n.outside.cores, type = "SOCK") 
  Packages.wrap <- c("gdistance", "rgdal", "ggplot2", "reshape2", "gridExtra", "ggspatial")
  clusterExport(cl.wrap, "Packages.wrap", envir = environment())
  clusterEvalQ(cl.wrap, lapply(Packages.wrap, library, character.only = TRUE)) # export needed library

  clusterExport(
    cl.wrap,
    c(
      "MakeL_tStack2",
      "NormcorrListToLikelihoodSurf",
      "PlotTracks2",
      "ModelListToTransMatList",
      "ModelListToTransMatListToCombinedTransMat2",
      "NormcorrListToProbList",
      "CalcNormCorr2",
      "MakeModelList2",
      "TransMatToNLL2",
      "ObjFunc2",
      "LLProfileOptimize2",
      "TrackSPDFtoTrackSPLFLines" 
    ),
    envir = environment())  # this specifies that objects exported will come first within the function
  clusterExport(cl.wrap, c("startvals", 
                           #"param.factor", 
                           "directions"), envir=environment()) 
  clusterExport(cl.wrap, c("model.list.by.track.by.groupid", "covars"), envir=environment()) 
  
  df.list.par  <- parallel::parLapplyLB(cl=cl.wrap,
                                                X=model.list.by.track.by.groupid,
                                                fun=LLOptimize3,
                                                optim.start.par=startvals,
                                              #  param.factor=param.factor,
                                                folderpath=folderpath, 
                                                model.split1=model.split1,  # TODO: replace model.split1 with model.string, adjust dependencies
                                                Lt.tiff = Lt.tiff,
                                                covars= covars,
                                                n.inside.cores = n.inside.cores, 
                                                plotindivtracks = plotindivtracks)
  
  stopCluster(cl.wrap)

gc()

  short.out <- do.call(rbind,sapply(df.list.par, "[", 1))
  row.names(short.out) <- NULL
  short.out$groupid <- names(track.list.by.groupid)
  long.out <- do.call(rbind,sapply(df.list.par, "[", 2))
  row.names(long.out)<-NULL
  long.out$groupid <- rep(names(track.list.by.groupid), each=npars)
  
  
  combinedtracks <- do.call(rbind,track.list.by.groupid)
  ntracks<- length(unique(paste0(combinedtracks$groupid,combinedtracks$trackid)))
  ngroupids <-length(track.list.by.groupid)
  
  
  # next line edit 3/5
  parnames <- c("loglambda", model.split1)
  #  parnames <- c("loglambda", variables)
  sds <- as.data.frame(t(sapply(short.out[,parnames], sd, na.rm = TRUE) ))
  means <- as.data.frame(t(sapply(short.out[,parnames], mean, na.rm = TRUE) ))
  
  
  combinednll <- sum(short.out$nll)
  allconverge <- sum(short.out$convergence)==0
  outvals <- c(combinednll, means$loglambda, sds$loglambda, allconverge)
  paramnames <-c("nll_all_groups", "loglambda.mean", "loglambda.sd", "alltrackconverge")
  df.summary<-as.data.frame(t(as.data.frame(outvals, row.names = paramnames)))
  tracklevelall <- long.out
  tracklevelall$model <- model.string
  tracklevelall$dataset <- dataset # dataset needs to be changed to groupid
  tracklevelall.out.name <- paste(folderpath, "beta_by_groups.csv", sep="/")
  write.csv(tracklevelall, tracklevelall.out.name, row.names = F)
  reff.out <-rbind(means, sds)
  reff.out$descr <- c("beta_group", "stddev.of.beta_groups")
  for(pn in parnames){
    cis<-long.out[long.out$vars==pn,1:4]
    track.cis.avgs <- apply(cis, 2, FUN=mean)      
    cisonly <- (c(track.cis.avgs["lCI"],track.cis.avgs["uCI"]))
    df2 <- data.frame(z=cisonly)
    colnames(df2)<-pn
    if(pn == parnames[1]){
      dfout<-df2
    } else{
      dfout <-cbind(dfout, df2)
    }
  }
  dfout$descr <- c("lower95CI_mean_of_hessian_method", "upper95CI_mean_of_hessian_method")
  reff.out <- rbind(reff.out, dfout)
  #### insert CI calculation from the standard error of the mean
  inline <- reff.out
  sigma <- inline[2,]
  sigma <- sigma[-length(sigma)]
  mu <- inline[1,]
  mu <- mu[-length(mu)]
  sem <- sigma/(sqrt(ngroupids))
  uCI <- mu+(1.96*sem)
  lCI <- mu-(1.96*sem)
  sem2 <- sem 
  sem2$descr <- "std.err.mean"
  uCI2 <- uCI
  uCI2$descr <- "upper95CI_SEMmethod"
  lCI2 <- lCI
  lCI2$descr <- "lower95CI_SEMmethod"
  signif <- uCI 
  signif[1,] <- 0
  signif[uCI>0 & lCI>0] <- 1 
  signif[uCI<0 & lCI<0] <- 1
  signif$descr <- "signif95CI_SEMmethod"
  coefs.with.cis <- rbind(uCI2, inline[1,], lCI2, signif, inline[-1,], sem2)
  write.csv(coefs.with.cis, 
            paste0(folderpath, "/", "betas_all_groups_w_95CIs.csv"),
            row.names = F)
  #### end CI calculation
  # now calculate AIC
  k <- npars # nvars plus one for the shape parameter 
  # 2k - 2 log(L) AIC formula
  # 2k + 2*nll AIC formula for nll, where nll is the negative log likelihood
  AIC <- 2*k + 2*combinednll
  #  AICc
  # AIC + (2k^2 + 2k / (n-k-1))
  AICc <- AIC + ((2*(k^2) + 2*k)/(ngroupids-k-1))
  #### end AIC calculation
  # add AIC and AICc to df.summary
  df.summary$AIC <- AIC
  df.summary$AICc <- AICc
  df.summary$k <- k
  df.summary$ntracks <- ntracks
  df.summary$ngroupids <- length(track.list.by.groupid) # number of animal ids
  df.summary$model <- model.string
  df.summary$dataset <- dataset
  df.summary$resolution <- res(covars)[1] # assumes square cells
  df.summary$directions <- directions
  # write df.summary dataframe out
  out.allinfo.table.name <- paste(folderpath, "AIC_etc.csv", sep="/")
  write.csv(df.summary, out.allinfo.table.name, row.names = F)
  # BOXPLOT of beta_hat coefficient estimates 
  tracklevel.plot <- tracklevelall[tracklevelall$vars!="loglambda",]
  plot.coefs <- tracklevel.plot %>% 
    ggplot(aes(x=factor(vars),y=est)) +
    geom_boxplot(outlier.shape=NA) + 
    geom_jitter(width = 0.1, alpha=0.3)+  #  labs(fill = "herd") + 
    xlab("Variable")+
    ylab("Coefficient Estimate")+
    theme_bw(base_size = 16)+
    theme(axis.text.x = element_text(angle = 45, hjust = 1))+
    geom_hline(yintercept = 0, linetype='dotted', col='black')+
    ggtitle(paste0(dataset, "\n", model.string), subtitle = paste0("AIC:", round(AIC, digits=1)))+
    theme(plot.title = element_text(size =10, face = "bold"),
          plot.subtitle = element_text(size =8, face = "bold"))
  png(paste0(folderpath, "/betas_boxplot.png"),
      width=4*1.62, height=4+1,
      units="in",
      bg="transparent", res=300)
  print(plot.coefs)
  dev.off()  


  if(Lt.tiff){
    Lt.out <- sapply(df.list.par, "[", 3)
    Lt.stack <-raster::stack(Lt.out)
    
    Lt.stack.terra <- terra::rast(Lt.stack)
    terra::writeRaster(x = Lt.stack.terra, 
                        filename = paste0(folderpath, "/Lt_multiband.tif"), 
                        overwrite=TRUE)
    
    Lh.raster <-sum(Lt.stack, na.rm=T)
    Lh.raster.terra <- terra::rast(Lh.raster)
    
    terra::writeRaster(Lh.raster.terra, paste0(folderpath, "/Lh_singleband.tif"), overwrite=TRUE)
    
    tracks<-do.call(rbind, track.list.by.groupid) # this puts all tracks from all animals together for plotting
    plottracks <- tracks[,"trackid"]
    names(plottracks) <- "id"
    # tracks.spdf to lines
    l <- TrackSPDFtoTrackSPLFLines(plottracks) 
    # convert to dataframes
    Lcorr.ras <- ggspatial::df_spatial(Lh.raster)
    trackpoint.gg <- ggspatial::df_spatial(plottracks)
    trackline.gg <- ggspatial::df_spatial(l, feature_id=id)
    # generate plot
    L_t.plot.simp <-ggplot()+
      geom_raster(data=Lcorr.ras, aes(x=x, y=y, fill=band1))+
      coord_fixed()+
      scale_fill_viridis_c(option = "viridis", direction=1, begin = .3)+
      geom_path(data=trackline.gg, aes(x,y, group=id), 
                size = 0.5, colour="black",  
                linetype = "dashed", alpha=0.3)+
      #      geom_point(data=trackpoint.gg, aes(x,y,group=id), colour="red", size=0.1)+
      theme_bw()+
      labs(fill="Likelihood")+
      ggtitle(label=dataset, subtitle = paste0(model.string , "\n", paste0("AIC:", round(AIC, digits=1))))
    png(paste0(folderpath, "/L_h.png"),
        width=180/25.4, height=180/25.4,
        units="in",
        bg="transparent", res=300)
    print(L_t.plot.simp)
    dev.off()  
    
    rm(Lt.out)
  }

rm(df.list.par)  
  }

# Objective function for calculating track likelihood for use in optim
# X are the parameters
ObjFunc2 <- function(X, model.list.by.track, #param.factor, 
                     cl.inside){
  model.list.working <- model.list.by.track
  pars0 <- as.numeric(X)  # pars will be a vector of:c(rate or sigma, coef1, coef2, ...) 
  print(pars0)
  pars <- pars0 #*param.factor
  lambda <- exp(pars[1])  #this changes loglambda to lambda
  
  for(t in 1:length(model.list.working)){
   # print(t)
  for(i in 2:length(pars)){
  #  print(pars[i])
    model.list.working[[t]]$focal.model.list[[i-1]][[4]]<-pars[i]
  }
  }
  
  trans.mat.by.track <- ModelListToTransMatListToCombinedTransMat2(model.list.by.track = model.list.working, directions=directions)  
  tracknll <- parallel::parLapplyLB(cl=cl.inside,
                                    X=trans.mat.by.track, 
                                    TransMatToNLL2, rate=lambda)
  groupnll <- sum(unlist(tracknll), na.rm=T)
  print(c(pars, groupnll))
  return(groupnll)
} # end ObjFunc2

## 2. TransMatToNLL2
# calculates a normalized corridor from a transition matrix, then calculates likelihood of points based on
# an exponential distribution.(combination of old CalcNormCorr and NormcorrListToProbList.lap)
TransMatToNLL2 <- function(X, rate){ #X is list of tracks with individual trans.mat
  track <- X[[2]]
  trans.mat <- X[[1]]
  startpoint <- track[1,] 
  endpoint <- track[length(track),]
  corr <- gdistance::accCost(trans.mat, startpoint) + gdistance::accCost(trans.mat, endpoint)
  norm.corr <- corr - min(values(corr), na.rm=T) 
  corr.ys <- values(norm.corr)
  exp.prob <- dexp(corr.ys, rate=rate)
  prob.corr <- norm.corr
  prob.corr[] <- exp.prob
  droppoints <- c(1, length(track)) # drop endpoints
  point.probs <- raster::extract(prob.corr, track[-droppoints,])
  probsAll <- values(prob.corr)
  ll.multinom <- log(point.probs)-log(sum(probsAll)) 
  nll.multinom <- -sum(ll.multinom)
  return(nll.multinom)
}

## 3. ModelListToTransMatListToCombinedTransMat2
# Takes the model list and transforms it into a transition matrix for calculating cost distances
ModelListToTransMatListToCombinedTransMat2 <- function(model.list.by.track, directions){

  for(t in 1:length(model.list.by.track)){
  model.list <- model.list.by.track[[t]][[1]]
  trans.mat.list <- lapply(model.list, ModelListToTransMatList, directions=directions) # takes a while...
  nulls <- unlist(lapply(trans.mat.list,function(X){length(X@transitionMatrix@x)==0}))
  trans.mat.list.nonulls <- trans.mat.list[!nulls]
  if(length(trans.mat.list.nonulls) > 0){
    for(i in 1:length(trans.mat.list.nonulls)){
      foc <-trans.mat.list.nonulls[[i]]
      if(max(foc@transitionMatrix@x[foc@transitionMatrix@x != Inf]) != -Inf){
        foc@transitionMatrix@x[foc@transitionMatrix@x==Inf] <- max(foc@transitionMatrix@x[foc@transitionMatrix@x!=Inf]) }# set any Inf cells to max of non-Inf cells
      if(i==1){
        tm.mult <- foc
      } else {
        tm.mult <- tm.mult * foc # multiply transition matrices together ... testing 12/10
      }
    }
    tm.mult2 <- tm.mult
  } else {  # this is what happens if the transition matrix is uniform ...
    ras <-model.list[[1]][[1]]
    ras[]<-1
    tm.mult2 <-gdistance::transition(ras, mean, directions = directions)
  }
  tm.corr <- gdistance::geoCorrection(tm.mult2, type="c")
  tm.corr@transitionMatrix@x[tm.corr@transitionMatrix@x==0]<-.Machine$double.xmin
  
  track <- model.list.by.track[[t]][[2]]
  
  sublist <- list(tm.corr, track)
  if(t==1){
    outlist <- list(sublist)
  } else {
    outlist[[t]]<-sublist
  }
  }
  # test 3/5/2023
  rm(tm.corr, tm.mult2, tm.mult, track, sublist)
  # end test
  return(outlist)
} # end ModelListToTransMatListToCombinedTransMat2

## 4. ModelListToTransMatList
ModelListToTransMatList <- function(X, directions){
  covar <- X[[1]]
  form <- X[[2]]
  formtype <- X[[3]]
  if(formtype == "isotropic"){
    coef <- X[[4]]
    cond <- exp(eval(parse(text=form)))
    tr1 <- gdistance::transition(cond, transitionFunction=mean, directions=directions)
    # geoCorrection for isotropic happens after transition matrices are combined
    # in first part of the TrackNLL.ObjFunc function
  } else if (formtype == "anisotropic"){
    coef <- X[[4]]
    transition.matrix <- gdistance::transition(covar, transitionFunction=function(x){x[2] - x[1]}, directions=directions, symm=FALSE) # calculate a transition matrix
    delt <- gdistance::geoCorrection(transition.matrix, type="c") 
    adj <- raster::adjacent(covar, cells=1:ncell(covar), pairs=TRUE,  directions=directions)
    speed <- delt
    speed[adj] <-  exp(eval(parse(text=form))) # note that "form" is "(abs(delt[adj])*coef)"
    tr1 <- speed
  } 
}

## XX. MakeModelList
MakeModelList <- function(model.string, covars#, covar.standardize
                          ){  
  model.split1 <- unlist(strsplit(model.string, split = "\\+"))
  formulas <- sapply(strsplit(model.split1, "\\("),"[[", 1)
  ms2 <- sapply(strsplit(model.split1, "\\(") ,"[[", 2)
  variables <- sapply(strsplit(ms2, "\\)"),"[[", 1)
  nvars <- length(variables)
  tropy.string <- formulas
  tropy.string[formulas == "iso"] <- "isotropic"
  tropy.string[formulas == "aniso"] <- "anisotropic"
  formula.string <- formulas
  formula.string[formulas == "iso"] <- "coef*covar"
  formula.string[formulas == "aniso"] <- "(abs(delt[adj])*coef)"
  for(n in 1:nvars){
    foc.covar <- covars[[as.character(variables[n])]]
    # if(covar.standardize[n]==1){
    #   foc.covar <- scale(foc.covar, center=F)
    # }
    l.out <- list(foc.covar,
                  formula.string[n],
                  tropy.string[n])  
    if(n==1){
      focal.model.list<- list(l.out)
    } else{
      focal.model.list[[n]]<-l.out
    }
  } # end generate model list
  return(list(focal.model.list = focal.model.list))
}  # end MakeModelList



## 5.CalcNormCorr2
# calculates a normalized corridor from a transition matrix, returns a list of tracks and the normalized corridor
CalcNormCorr2 <- function(X){ #X is trans.mat.by.track, the combined list of tracks and trans.mat, use with Lapply
  track <- X[[2]]
  trans.mat <- X[[1]]
  startpoint <- track[1,] 
  endpoint <- track[length(track),]
  AC <- gdistance::accCost(trans.mat, startpoint)
  BC <- gdistance::accCost(trans.mat, endpoint)
  corr <- AC + BC #(A.1 + B.1)/2 # calculates least cost corridor
  D <- corr -  min(values(corr), na.rm=T) 
  out.list <- list(track, D) # returns a list of tracks with the normalized corridor...
  return(out.list)
}

## 6. NormcorrListToProbList was NormcorrListToProbList.lap
NormcorrListToProbList <- function(X, rate){ 
  norm.corr <- X[[2]] 
  onetrack <- X[[1]]
  startpoint<- onetrack[1,]
  endpoint <- onetrack[nrow(onetrack),]
  corr.ys <- values(norm.corr)
  exp.prob <- dexp(corr.ys, rate=rate)
  prob.corr <- norm.corr
  prob.corr[] <- exp.prob
  point.probs <- raster::extract(prob.corr, onetrack)
  probsAll <- values(prob.corr)
  ll.multinom <- log(point.probs)-log(sum(probsAll)) 
  nll.multinom <- -sum(ll.multinom)
  return(nll.multinom)
}

PlotTracks2 <- function(pars, model.list.by.track, #param.factor, 
                        cl.inside=NULL,
                       writeindivtracksdir=NULL){
  model.list.working <- model.list.by.track
  pars <- pars#*param.factor
  loglambda <- pars[1]#exp(pars[1])  #this changes loglambda to lambda
  
  for(t in 1:length(model.list.working)){
    for(i in 2:length(pars)){
      model.list.working[[t]]$focal.model.list[[i-1]][[4]]<-pars[i]
    }
  }
  
  trans.mat.by.track <- ModelListToTransMatListToCombinedTransMat2(model.list.by.track = model.list.working, directions=directions)  
  if(is.null(cl.inside)){
  norm.corr.by.track <- lapply(X=trans.mat.by.track, FUN=CalcNormCorr2)
      } else {
    norm.corr.by.track <- parallel::parLapplyLB(cl=cl.inside,
                                       X=trans.mat.by.track, 
                                       fun=CalcNormCorr2)
  }
  
  L_t <- lapply(X=norm.corr.by.track, FUN=NormcorrListToLikelihoodSurf, loglambda=loglambda)
  
  for(i in 1:length(model.list.by.track)){
    
    onetrack <- model.list.by.track[[i]][[2]]
    onetrack <- onetrack[,"trackid"]
    names(onetrack) <- "id"
    l <- TrackSPDFtoTrackSPLFLines(onetrack) 
    Lcorr.ras <- ggspatial::df_spatial(L_t[[i]])
    trackpoint.gg <- ggspatial::df_spatial(onetrack)
    trackline.gg <- ggspatial::df_spatial(l, feature_id=id)
    # generate plot
    L_t.plot.simp <-ggplot()+
      geom_raster(data=Lcorr.ras, aes(x=x, y=y, fill=band1))+
      coord_fixed()+
      scale_fill_viridis_c(option = "viridis", direction=1, begin = .3)+
      geom_path(data=trackline.gg, aes(x,y, group=id), 
                size = 0.05, colour="black", 
                linetype = "dashed")+
      geom_point(data=trackpoint.gg, aes(x,y,group=id), colour="red", size=0.1)+
      theme_bw()+
      labs(fill="Likelihood")
    # write plot
    png(paste(writeindivtracksdir,"/trackplot_", unique(model.list.by.track[[i]][[2]]$groupid),"_",unique(model.list.by.track[[i]][[2]]$trackid), ".png", sep=""),
        width=7, height=7, 
        units="in", pointsize =8,
        bg="transparent", res=300)
    print(L_t.plot.simp)
    dev.off()
  }
}

NormcorrListToLikelihoodSurf <- function(X, loglambda){   #note here rate is loglambda ... need to fix this
  D_t <- X[[2]]
  onetrack <- X[[1]]
  startpoint<- onetrack[1,]
  endpoint <- onetrack[nrow(onetrack),]
  corr.ys <- values(D_t)
  exp.prob <- dexp(corr.ys, rate=exp(loglambda))
  P_t <- D_t
  P_t[] <- exp.prob
  L_t <- P_t/(sum(P_t[]))
  return(L_t)
}

MakeL_tStack2 <- function(pars, model.list.by.track, #param.factor, 
                          cl.inside=NULL, originalcovars,
                          directions){
  model.list.working <- model.list.by.track
  pars <- pars#*param.factor
  loglambda <- pars[1]
  for(t in 1:length(model.list.working)){
    for(i in 2:length(pars)){
      model.list.working[[t]]$focal.model.list[[i-1]][[4]]<-pars[i]
    }
  }
  trans.mat.by.track <- ModelListToTransMatListToCombinedTransMat2(model.list.by.track = model.list.working, directions=directions)  
  if(is.null(cl.inside)){
    norm.corr.by.track <- lapply(X=trans.mat.by.track, FUN=CalcNormCorr2)
  } else {
    norm.corr.by.track <- parallel::parLapplyLB(cl=cl.inside,
                                                X=trans.mat.by.track, 
                                                fun=CalcNormCorr2)
  }
  
  L_t.list <- lapply(X=norm.corr.by.track, FUN=NormcorrListToLikelihoodSurf, loglambda=loglambda)
  
extended <-lapply(L_t.list, FUN=raster::extend, y=originalcovars)
L_t.stack <-raster::stack(extended)
return(L_t.stack)
}

TrackSPDFtoTrackSPLFLines <- function(tracks.spdf){
  x <- lapply(split(tracks.spdf, tracks.spdf$id), function(x) Lines(list(Line(coordinates(x))), x$id[1L]))
  lines <- SpatialLines(x)
  data <- data.frame(id = unique(tracks.spdf$id))
  rownames(data) <- data$id
  l <- SpatialLinesDataFrame(lines, data)
}

## XX. MakeModelList2 - this one crops by each track, requires lapply to return list of model lists...
# cropfactor sets the proportion increase in the extent of the track that is used to crop the env covariates
MakeModelList2 <- function(X, cropfactor=3, model.string, covars
                           #, covar.standardize
                           ){  
  track <- X
  if(!is.null(cropfactor)){
  cropextent <- raster::extent(track)*cropfactor
  covars <- raster::crop(covars, cropextent)
  }
  model.split1 <- unlist(strsplit(model.string, split = "\\+"))
  formulas <- sapply(strsplit(model.split1, "\\("),"[[", 1)
  ms2 <- sapply(strsplit(model.split1, "\\(") ,"[[", 2)
  variables <- sapply(strsplit(ms2, "\\)"),"[[", 1)
  nvars <- length(variables)
  tropy.string <- formulas
  tropy.string[formulas == "iso"] <- "isotropic"
  tropy.string[formulas == "aniso"] <- "anisotropic"
  formula.string <- formulas
  formula.string[formulas == "iso"] <- "coef*covar"
  formula.string[formulas == "aniso"] <- "(abs(delt[adj])*coef)"
  for(n in 1:nvars){
    foc.covar <- covars[[as.character(variables[n])]]
    # if(covar.standardize[n]==1){
    #   foc.covar <- scale(foc.covar, center=F)
    # }
    l.out <- list(foc.covar,
                  formula.string[n],
                  tropy.string[n])  
    if(n==1){
      focal.model.list<- list(l.out)
    } else{
      focal.model.list[[n]]<-l.out
    }
  } # end generate model list
  return(list(focal.model.list = focal.model.list, track=track))
}  # end MakeModelList2




LLOptimize3 <- function (X,
                                optim.start.par,
                                # param.factor,
                                n.inside.cores, 
                                folderpath, 
                                model.split1,  # TODO: replace model.split1 with model.string, adjust dependencies
                                Lt.tiff,
                                covars,
                         plotindivtracks) {
  model.list.by.track <- X
  cl.inside <- makeCluster(n.inside.cores) #
  Packages.inside <- c("gdistance", "ggplot2", "ggspatial")
  clusterExport(cl.inside, "Packages.inside",envir = environment())
  clusterEvalQ(cl.inside, lapply(Packages.inside, library, character.only = TRUE)) # export needed library
  
  clusterExport(
    cl.inside,
    c(
      "MakeL_tStack2",
      "NormcorrListToLikelihoodSurf",
      "PlotTracks2",
      "ModelListToTransMatList",
      "ModelListToTransMatListToCombinedTransMat2",
      "NormcorrListToProbList",
      "CalcNormCorr2",
      "MakeModelList2",
      "TransMatToNLL2",
      "ObjFunc2",
      #"LLProfile2",
      "LLProfileOptimize2",
      "TrackSPDFtoTrackSPLFLines" 
    ),
    envir = environment())  # this specifies that objects exported will come first within the function
  clusterExport(cl.inside, c("model.list.by.track", "startvals",
                             #"param.factor", 
                             "directions"), envir=environment()) 
  
  # repeatcounter <- 0
  # repeat{  
  #   global.optim <- NULL
  #   repeatcounter <- repeatcounter + 1
  #   
  #   #     test <- ObjFunc2(X=optim.start.par,
  #   # #                    track.list=track.list,
  #   #                     model.list.by.track=model.list.by.track,
  #   #                     param.factor = param.factor,
  #   #                     cl.inside = cl.inside)
  #   #     
  #   
  global.optim <- optim(optim.start.par, ObjFunc2,
                          method="Nelder-Mead",
                          model.list.by.track=model.list.by.track,
                          #param.factor=param.factor,
                          cl.inside=cl.inside,
                          hessian = T,
                          control = list(maxit=200))
  
  gc()
  #   
  #   # calculate likelihood profile, see if can find a better nll than optim
  #   LL_list <- LLProfile2(fittedpars=global.optim$par,
  #                         model.list.by.track = model.list.by.track, 
  #                         #track.list=track.list, 
  #                         formulas = model.split1,
  #                         cl.inside=cl.inside, 
  #                         # param.factor=param.factor,
  #                         folderpath = folderpath, 
  #                         searchgrid = exp(seq(-6,3,.5))) # third element of seq() adjusts grain, TODO make this a variable
  #   LLprofile.df <- do.call(rbind, LL_list[[2]])
  #   
  #   # uncomment next line to write out LL profile data
  #   write.csv(LLprofile.df, paste(folderpath, "/LLprofiledata_", unique(model.list.by.track[[1]]$track$groupid), "_", repeatcounter,".csv", sep=""))
  #   
  #   LL.CIs.df <-LL_list[[1]]
  #   LL.profile.nll <- min(LL.CIs.df$minnll, na.rm=T)
  #   
  #   globminind<-which.min(LLprofile.df$nlls.avec1)
  #   optim.start.par <- as.numeric(LLprofile.df[globminind,1:length(optim.start.par)])
  #   
  #   if(LL.profile.nll > (global.optim$value - 0.1) | repeatcounter > 10){
  #     break
  #   }  
  # }# end repeat loop  
  
  fit <- global.optim  
  vc <-  tryCatch({solve(fit$hessian)}, error = function(e){vc <- NA})
  SE <-   tryCatch({SE <- sqrt(diag(vc))}, error = function(e){SE <- NA} )
  
  lower.raw <- fit$par-1.96*SE
  upper.raw <- fit$par+1.96*SE
  
  parnames <- c("loglambda", model.split1)
  
  df.pars <- t(as.data.frame(fit$par, row.names = parnames))
  df.pars <- as.data.frame(df.pars)
  rownames(df.pars)<-NULL
  
  df.out<-data.frame(nll=fit$value, convergence=fit$convergence)
  df.out <- cbind(df.out, df.pars)
  
  parnames <- c("loglambda", model.split1)
  
  df.e <- as.data.frame(fit$par, row.names = parnames)
  df.l <- as.data.frame(lower.raw, row.names = parnames)
  df.u <- as.data.frame(upper.raw, row.names = parnames)
  df.SE <- as.data.frame(SE, row.names = parnames)
  df.c <-cbind(df.l, df.e, df.u, df.SE)
  names(df.c)<- c("lCI", "est", "uCI", "SE")
  df.c$vars <- parnames
  df.c$nll <- fit$value
  
  if(plotindivtracks){
  writeindivtracksdir=paste(folderpath, "indivtrackfits", sep="/")
  if(!dir.exists(writeindivtracksdir)){
    dir.create(writeindivtracksdir)}
  
  PlotTracks2(pars=as.numeric(df.pars), 
              model.list.by.track = model.list.by.track,
              #param.factor=param.factor,
              cl.inside=cl.inside,
              writeindivtracksdir=writeindivtracksdir)
  }
  
  if(Lt.tiff){ #Lt.tiff is a boolean
    
    L_t.stack <- MakeL_tStack2(pars=as.numeric(df.pars), 
                               model.list.by.track = model.list.by.track,
                               #param.factor=param.factor,
                               cl.inside=cl.inside,
                               originalcovars = covars,
                               directions = directions)
    
  } else {
    L_t.stack <- NULL
  }
  
  stopCluster(cl.inside)
  
  return(list(df.out, df.c, L_t.stack))
} 










####### DEPRECATED FUNCTIONS
# LLProfileOptimize2 is old function used for calculating profiles. 
LLProfileOptimize2 <- function (X,
                                optim.start.par,
                                # param.factor,
                                n.inside.cores, 
                                folderpath, 
                                model.split1,  # TODO: replace model.split1 with model.string, adjust dependencies
                                Lt.tiff = F,
                                covars) {
  model.list.by.track <- X
  cl.inside <- makeCluster(n.inside.cores) #
  Packages.inside <- c("gdistance", "ggplot2", "ggspatial")
  clusterExport(cl.inside, "Packages.inside",envir = environment())
  clusterEvalQ(cl.inside, lapply(Packages.inside, library, character.only = TRUE)) # export needed library
  
  clusterExport(
    cl.inside,
    c(
      "MakeL_tStack2",
      "NormcorrListToLikelihoodSurf",
      "PlotTracks2",
      "ModelListToTransMatList",
      "ModelListToTransMatListToCombinedTransMat2",
      "NormcorrListToProbList",
      "CalcNormCorr2",
      "MakeModelList2",
      "TransMatToNLL2",
      "ObjFunc2",
      
      #"LLProfile2",
      "LLProfileOptimize2",
      "TrackSPDFtoTrackSPLFLines" 
    ),
    envir = environment())  # this specifies that objects exported will come first within the function
  clusterExport(cl.inside, c("model.list.by.track", "startvals",
                             #"param.factor", 
                             "directions"), envir=environment()) 
  
  # repeatcounter <- 0
  # repeat{  
  #   global.optim <- NULL
  #   repeatcounter <- repeatcounter + 1
  
  #     test <- ObjFunc2(X=optim.start.par,
  # #                    track.list=track.list,
  #                     model.list.by.track=model.list.by.track,
  #                     param.factor = param.factor,
  #                     cl.inside = cl.inside)
  #     
  
  global.optim <- optim(optim.start.par, ObjFunc2, 
                        method="Nelder-Mead",
                        model.list.by.track=model.list.by.track,
                        #param.factor=param.factor,
                        cl.inside=cl.inside,
                        hessian = T,
                        control = list(maxit=200))
  
  # calculate likelihood profile, see if can find a better nll than optim
  # LL_list <- LLProfile2(fittedpars=global.optim$par,
  #                           model.list.by.track = model.list.by.track, 
  #                           #track.list=track.list, 
  #                           formulas = model.split1,
  #                           cl.inside=cl.inside, 
  #                          # param.factor=param.factor,
  #                           folderpath = folderpath, 
  #                           searchgrid = exp(seq(-6,3,.5))) # third element of seq() adjusts grain, TODO make this a variable
  # LLprofile.df <- do.call(rbind, LL_list[[2]])
  # 
  # uncomment next line to write out LL profile data
  #   write.csv(LLprofile.df, paste(folderpath, "/LLprofiledata_", unique(model.list.by.track[[1]]$track$groupid), "_", repeatcounter,".csv", sep=""))
  #   
  #   LL.CIs.df <-LL_list[[1]]
  #   LL.profile.nll <- min(LL.CIs.df$minnll, na.rm=T)
  #   
  #   globminind<-which.min(LLprofile.df$nlls.avec1)
  #   optim.start.par <- as.numeric(LLprofile.df[globminind,1:length(optim.start.par)])
  #   
  #   if(LL.profile.nll > (global.optim$value - 0.1) | repeatcounter > 10){
  #     break
  #   }  
  # }# end repeat loop  
  
  fit <- global.optim  
  vc <-  tryCatch({solve(fit$hessian)}, error = function(e){vc <- NA})
  SE <-   tryCatch({SE <- sqrt(diag(vc))}, error = function(e){SE <- NA} )
  
  lower.raw <- fit$par-1.96*SE
  upper.raw <- fit$par+1.96*SE
  
  parnames <- c("loglambda", model.split1)
  
  df.pars <- t(as.data.frame(fit$par, row.names = parnames))
  df.pars <- as.data.frame(df.pars)
  rownames(df.pars)<-NULL
  
  df.out<-data.frame(nll=fit$value, convergence=fit$convergence)
  df.out <- cbind(df.out, df.pars)
  
  parnames <- c("loglambda", model.split1)
  
  df.e <- as.data.frame(fit$par, row.names = parnames)
  df.l <- as.data.frame(lower.raw, row.names = parnames)
  df.u <- as.data.frame(upper.raw, row.names = parnames)
  df.SE <- as.data.frame(SE, row.names = parnames)
  df.c <-cbind(df.l, df.e, df.u, df.SE)
  names(df.c)<- c("lCI", "est", "uCI", "SE")
  df.c$vars <- parnames
  df.c$nll <- fit$value
  
  if(plotindivtracks){
    writeindivtracksdir=paste(folderpath, "indivtrackfits", sep="/")
    if(!dir.exists(writeindivtracksdir)){
      dir.create(writeindivtracksdir)}
    
    PlotTracks2(pars=as.numeric(df.pars), 
                model.list.by.track = model.list.by.track,
                #param.factor=param.factor,
                cl.inside=cl.inside,
                writeindivtracksdir=writeindivtracksdir)
  }
  
  if(Lt.tiff){ #Lt.tiff is a boolean
    L_t.stack <- MakeL_tStack2(pars=as.numeric(df.pars), 
                               model.list.by.track = model.list.by.track,
                               #param.factor=param.factor,
                               cl.inside=cl.inside,
                               originalcovars = covars,
                               directions = directions)
    
  } else {
    L_t.stack <- NULL
  }
  
  stopCluster(cl.inside)
  
  return(list(df.out, df.c, L_t.stack))
} 


#########################################################
# CROSSVALIDATION FUNCTIONS
#########################################################


#' XVal takes a dataframe of tracks and a model and returns an ROC object
XVal <- function(tracks,
                 fittedbetas,
                 model.string,
                 covars,
                 cropfactor,
                 ncores, 
                 directions){
  
  track.list <- split(tracks, tracks$trackid)
  model.list.by.track <- lapply(X=track.list, FUN=MakeModelList2,
                                model.string = as.character(model.string),
                                covars = covars,
                                cropfactor = cropfactor)
  cl.inside <- makeCluster(ncores) #
  Packages.inside <- c("gdistance")
  clusterExport(cl.inside, "Packages.inside", envir = environment())
#  clusterExport(cl.inside, "Packages.inside")
  clusterEvalQ(cl.inside, lapply(Packages.inside, library, character.only = TRUE)) # export needed library -- takes a long time; possible to 
  L_t.stack <- MakeL_tStack2(pars=as.numeric(fittedbetas), 
                             model.list.by.track = model.list.by.track,
                             cl.inside=cl.inside,
                             originalcovars = covars, 
                             directions = directions)
  roc.obj <- calc_ROC(probabilitycorridorstack = L_t.stack,
           tracks = tracks)
  stopCluster(cl.inside)
  return(roc.obj)
} # end XVal

#Deprectated
#PctPointsPerPercentileByTrack <- function(probabilitycorridorstack, tracks){
#   ecdf_fun <- function(x,perc) ecdf(x)(perc)
#   track.list <- split(tracks, tracks$trackid)
#   for(l in 1:length(track.list)){
#     foc.track <- track.list[[l]]
#     droppoints <- c(1, nrow(foc.track)) # drop endpoints
#     foc.track <- foc.track[-droppoints,] # remove first and last point
#     probabilitycorridor <- probabilitycorridorstack[[l]]
#     raster.percentiles <- ecdf_fun(values(probabilitycorridor),values(probabilitycorridor))
#     perc.ras <- probabilitycorridor
#     perc.ras[] <- 1-raster.percentiles
# #    plot(perc.ras)
#     vals.tracks <- raster::extract(perc.ras, foc.track) # landscape % at each point
#     track.percentiles <- ecdf_fun(vals.tracks,vals.tracks) # % of points in each landscape category
#     perc.df <- data.frame(prop.points= track.percentiles, prop.landscape=vals.tracks)
#     perc.df$trackid <- foc.track$trackid[1]
#     if(l==1){
#       perc.df.out <- perc.df
#     } else {
#       perc.df.out <- rbind(perc.df.out, perc.df)
#     }
#   }
#   return(perc.df.out)
# }    




calc_ROC <- function(probabilitycorridorstack, tracks){
  # depends on pROC
  library(pROC)
  track.list <- split(tracks, tracks$trackid)
  for(l in 1:length(track.list)){
    foc.track <- track.list[[l]]
    droppoints <- c(1, nrow(foc.track)) # drop endpoints
    foc.track <- foc.track[-droppoints,] # remove first and last point
    
    cases <- data.frame(vals=extract(probabilitycorridorstack[[l]], track.list[[l]]),
                        lab = 1, trackid=foc.track$trackid[1])
    controls <- data.frame(vals=values(probabilitycorridorstack[[l]]),
                           lab=0, trackid=foc.track$trackid[1])
    
    if(l==1){
      out.df <- rbind(cases,controls)
    } else {
      out.df <- rbind(out.df, rbind(cases, controls))
    }
  }

    roc.obj <-roc(cases = out.df$vals[out.df$lab==1], 
                   controls = out.df$vals[out.df$lab==0],
                   direction = "<")
    
    return(roc.obj)
}    

# test
# calc_ROC(probabilitycorridorstack=L_t.stack,
#         tracks=tracks)
# end test

PredictCorr <- function(points.list,
                        fittedbetas,
                        model.string,
                        covars,
                        directions) {

model.list.by.track <- lapply(X=points.list, 
         FUN=MakeModelList2,
         cropfactor = NULL,
         model.string=model.string,
         covars = covars)
  
  model.list.working <- model.list.by.track
  pars <- fittedbetas#*param.factor
  loglambda <- pars[1]
  for(t in 1:length(model.list.working)){
    for(i in 2:length(pars)){
      model.list.working[[t]]$focal.model.list[[i-1]][[4]]<-pars[i]
    }
  }
  trans.mat.by.track <- ModelListToTransMatListToCombinedTransMat2(model.list.by.track = model.list.working, directions=directions)  
  
  norm.corr.by.track <- lapply(X=trans.mat.by.track, FUN=CalcNormCorr2)
  
  L_t.list <- lapply(X=norm.corr.by.track, FUN=NormcorrListToLikelihoodSurf, loglambda=loglambda)
return(raster::stack(L_t.list))
}


