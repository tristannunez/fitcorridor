# Corridor Fitting Functions for Vignette Illustrating Statistical Corridor Modeling using Maximum Likelihood
#' To accompany "A statistical framework for predicting migration corridors using empirical data"
#' Nuñez, T., Hurley, M., Graves, T. Ortega, A. Sawyer, H., Fattebert, J., Merkle, J., and Kauffman, M.   
# 12 July 2022
# v1.0.0
# Tristan Nunez
# tristan.nunez@gmail.com

# # load needed libraries
library(gdistance)
library(parallel)
library(rgdal)
library(ggplot2)
library(ggspatial)
library(rgeos)

##################################################################
# CORE FUNCTIONS
##################################################################
#' Enumeration of functions:
#' FitSingleModel fits a single cost distance model (e.g., iso(HumanFootprint), which  
#' would be an isotropic model of the form conductance = exp(Beta*HumanFootprint) to a
#' set of movement tracks
#' 
#' 1. FitSingleModel depends on:
#'  2. MakeModelList
#'  3. FitModelTrackLevel.lap
#'  
#'  FitModelTrackLevel.lap depends on:
#'  4. TrackNLL.ObjFunc
#'  
#'  TrackNLL.ObjFunc depends on:
#'  5. CalcNormCorr
#'  6. NormCorrListtoProbList.lap
#'  7. ModelListtoTransMatlistToCombinedTransMat
#'  
#'  ModelListtoTransMatlistToCombinedTransMat depends on:
#'  8. ModelListtoTransMatlist

FitSingleModel <- function(tracks = tracks, 
                           covars = covars,
                           model.to.fit = NULL ,#"iso(SHRUB)+iso(HF)", # no spaces
                           neighbors = 8, # 8 or 16
                           covar.aggregate = 1 ,# 10, # factor by which to aggregate covariate rasters. run ?raster::aggregate for details
                           covar.standardize = NULL, #c(1,1), # same length as number of terms in model; 1 to standardize, 0 to not
                           param.factor = NULL, #c(1,1), # coefficient multiplier. 1 for no change. Generally 1 for iso() variables. For aniso() variables, try 10 or 100 for elevation (in meters), and 100 or 1000 for variables with Julian day units, such as snow-off date 
                           startval.string = NULL, #c(-15, 0, 0), # the first term is the start value for log(lambda); the others correspond to the model terms (in order)
                           dataset = "defaultherd", # name of track dataset 
                           out.folder.path = "./", #./fittedmods/"
                           n.inside.cores = 1, 
                           Lt.tiff = T 
                           ){
  model.string <- model.to.fit   # next 4 lines pull out info from model.string 
  model.split1 <- unlist(strsplit(as.character(model.string), split = "\\+"))
  ms2 <- sapply(strsplit(model.split1, "\\(") ,"[[", 2)
  variables <- sapply(strsplit(ms2, "\\)"),"[[", 1)
  nvars <- length(variables)
  npars <- nvars+1 # add 1 for lambda
  directions <- neighbors  # specify directions
  param_factor <-param.factor # extract scaling related parameters  
  covar_standardize <- covar.standardize
  # set up out files
  dataset <- as.character(dataset)
  if(!dir.exists(paste(out.folder.path, dataset, sep=""))){
    dir.create(paste(out.folder.path, dataset, sep=""))  
  }
  foldername <- model.to.fit
  folderpath <- paste(out.folder.path, dataset, "/", foldername, sep="") # might need to tweak
  dir.create(folderpath)
  params_nll_out_file_optim <- paste(folderpath, "/params_nll_out_optim.csv", sep="")
  covars <- aggregate(covars, fact=covar.aggregate)
  model.list.out <- MakeModelList(model.string = as.character(model.string),
                                  covars = covars, 
                                  covar_standardize = covar_standardize)
  focal.model.list <- model.list.out$focal.model.list
  track.list <- split(tracks, tracks$trackid)# split tracks, trackid as migration bout id or other separation...
  ntracks<- length(unique(tracks$trackid))
  track.startvals <- startval.string 
  cl.inside <- makeCluster(n.inside.cores) #
  Packages <- c("gdistance", "ggplot2", "ggspatial")
  clusterExport(cl.inside, "Packages",envir = environment())
  clusterEvalQ(cl.inside, lapply(Packages, library, character.only = TRUE)) # export needed library
  clusterExport(
    cl.inside,
    c("NormcorrListToProbList.lap",
      "ModelListToTransMatList",
      "ModelListToTransMatListToCombinedTransMat",
      "TrackNLL.ObjFunc",
      "CalcNormCorr",
      "LLProfileTrack",
      "PlotOneTrack",
      "NormcorrListToProbSurf.lap",
      "TrackSPDFtoTrackSPLFLines", 
      "GenerateL_t"
    ),
    envir = environment())  # this specifies that objects exported will come first within the function
  clusterExport(cl.inside, c("focal.model.list", "track.startvals", "variables", "param_factor", "directions"), envir=environment()) 
  dir.create(paste0(folderpath,"/likprofiles")) # create a folder for the likelihood profiles
  profile.optimize.out <- parallel::parLapplyLB(cl=cl.inside,
                                                X=track.list, 
                                                LLProfileOptimizeTrack, 
                                                startval.string=track.startvals,
                                                focal.model.list=focal.model.list,
                                                param_factor=param_factor,
                                                folderpath=folderpath,
                                                variables=variables, 
                                                model.split1 = model.split1,
                                                Lt.tiff = T)
  df.list.par <- profile.optimize.out
  short.out <- do.call(rbind,sapply(df.list.par, "[", 1))
  row.names(short.out) <- NULL
  long.out <- do.call(rbind,sapply(df.list.par, "[", 2))
  row.names(long.out)<-NULL
  parnames <- c("loglambda", variables)
  sds <- as.data.frame(t(sapply(short.out[,parnames], sd, na.rm = TRUE) ))
  means <- as.data.frame(t(sapply(short.out[,parnames], mean, na.rm = TRUE) ))
  tracknll <- sum(short.out$tracknll)
  allconverge <- sum(short.out$convergence)==0
  outvals <- c(tracknll, means$loglambda, sds$loglambda, allconverge)
  paramnames <-c("nll_h", "loglambda.mean", "loglambda.sd", "alltrackconverge")
  df.summary<-as.data.frame(t(as.data.frame(outvals, row.names = paramnames)))
  tracklevelall <- long.out
  tracklevelall$model <- model.string
  tracklevelall$dataset <- dataset
  tracklevelall.out.name <- paste(folderpath, "beta_ts.csv", sep="/")
  write.csv(tracklevelall, tracklevelall.out.name, row.names = F)
  reff.out <-rbind(means, sds)
  reff.out$descr <- c("beta_h", "stddev.of.beta_ts")
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
  sem <- sigma/(sqrt(ntracks))
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
            paste0(folderpath, "/", "beta_hs_w_95CIs.csv"),
            row.names = F)
  #### end CI calculation
  # now calculate AIC
  k <- npars # nvars plus one for the shape parameter 
  # 2k - 2 log(L) AIC formula
  # 2k + 2*nll AIC formula for nll, where nll is the negative log likelihood
  AIC <- 2*k + 2*tracknll
  #  AICc
  # AIC + (2k^2 + 2k / (n-k-1))
  AICc <- AIC + ((2*(k^2) + 2*k)/(ntracks-k-1))
  #### end AIC calculation
  # add AIC and AICc to df.summary
  df.summary$AIC <- AIC
  df.summary$AICc <- AICc
  df.summary$k <- k
  df.summary$ntracks <- ntracks
  df.summary$model <- model.string
  df.summary$dataset <- dataset
  df.summary$resolution <- res(covars)[1] # assumes square cells
  df.summary$neighbors <- neighbors
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
    Lt.stack <-stack(Lt.out)
    writeRaster(Lt.stack, paste0(folderpath, "/Lt_rasterstack.tif"), overwrite=TRUE)
    Lh.raster <-mean(Lt.stack)
    writeRaster(Lh.raster, paste0(folderpath, "/Lh_raster.tif"), overwrite=TRUE)
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
  }
  stopCluster(cl.inside)  
} # end wrapper function FitSingleModel

## 2. MakeModelList
MakeModelList <- function(model.string, covars, covar_standardize){  
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
    if(covar_standardize[n]==1){
      foc.covar <- scale(foc.covar, center=F)
    }
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

## 3.FitModelTrackLevel.lap
FitModelTrackLevel.lap <- function(X, startvals, model.list, variables, param_factor){ # where X is a list of tracks... # adding param_factor ... needed??? multiply startvals by param_factor???
  onetrack <- X
  trackfit <- optim(startvals, TrackNLL.ObjFunc, method="Nelder-Mead",
                    onetrack=onetrack, model.list=model.list, hessian=T,
                    param_factor=param_factor,
                    control = list(maxit=400)) # adding param_factor
  
  vc <-  tryCatch({solve(trackfit$hessian)}, error = function(e){vc <- NA})
  SE <-   tryCatch({SE <- sqrt(diag(vc))}, error = function(e){SE <- NA} )
  
  lower.raw <- trackfit$par-1.96*SE
  upper.raw <- trackfit$par+1.96*SE
  
  parnames <- c("loglambda", variables)
  df.pars <- t(as.data.frame(trackfit$par, row.names = parnames))
  df.pars <- as.data.frame(df.pars)
  rownames(df.pars)<-NULL
  
  df.out<-data.frame(trackid=unique(onetrack$trackid), tracknll=trackfit$value, convergence=trackfit$convergence)
  df.out <- cbind(df.out, df.pars)
  df.out
  
  # alternate arrangement:
  parnames <- c("loglambda", variables)
  df.e <- as.data.frame(trackfit$par, row.names = parnames)
  df.l <- as.data.frame(lower.raw, row.names = parnames)
  df.u <- as.data.frame(upper.raw, row.names = parnames)
  df.SE <- as.data.frame(SE, row.names = parnames)
  df.c <-cbind(df.l, df.e, df.u, df.SE)
  names(df.c)<- c("lCI", "est", "uCI", "SE")
  df.c$vars <- parnames
  df.c$trackid <- unique(onetrack$trackid)
  df.c$tracknll <- trackfit$value
  return(list(df.out, df.c))
}

## 4. TrackNLL.ObjFunc
# Objective function for calculating track likelihood for use in optim
TrackNLL.ObjFunc <- function(X, onetrack, model.list, param_factor){
  pars0 <- as.numeric(X)  # pars will be a vector of:c(lambda, coef1, coef2, ...) 
  print(pars0)
  pars <- pars0*param_factor
  lambda <- exp(pars[1])  
  for(i in 2:length(pars)){
    model.list[[i-1]][[4]]<-pars[i]
  }
  tm.corr <- ModelListToTransMatListToCombinedTransMat(model.list=model.list, directions=directions)  
  norm.corr <- CalcNormCorr(X=onetrack, trans.mat=tm.corr)
  tracknll <- NormcorrListToProbList.lap(X=norm.corr, rate=lambda)
  print(c(pars, tracknll))
  return(tracknll)
}

## 5.CalcNormCorr
# calculates a normalized corridor from a transition matrix, returns a list of tracks and the normalized corridor
CalcNormCorr <- function(X, trans.mat){ #X is list of tracks, trans.mat is transition matrix...
  startpoint <- X[1,] 
  endpoint <- X[length(X),]
  A.1 <- accCost(trans.mat, startpoint)
  B.1 <- accCost(trans.mat, endpoint)
  corr <- A.1 + B.1 # calculates least cost corridor
  nrml.corr <- corr -  min(values(corr), na.rm=T) 
  out.list <- list(X, nrml.corr) # returns a list of tracks with the normalized corridor...
  return(out.list)
}

## 6. NormcorrListToProbList.lap
NormcorrListToProbList.lap <- function(X, rate){ 
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

## 7. ModelListToTransMatListToCombinedTransMat
ModelListToTransMatListToCombinedTransMat <- function(model.list=model.list, directions=directions){
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
    tm.mult2 <-transition(ras, mean, directions = directions)
  }
  tm.corr <- geoCorrection(tm.mult2, type="c")
  tm.corr@transitionMatrix@x[tm.corr@transitionMatrix@x==0]<-.Machine$double.xmin
  return(tm.corr)
}

## 8. ModelListToTransMatList
ModelListToTransMatList <- function(X, directions=directions){
  covar <- X[[1]]
  form <- X[[2]]
  formtype <- X[[3]]
  if(formtype == "isotropic"){
    coef <- X[[4]]
    cond <- exp(eval(parse(text=form)))
    tr1 <- transition(cond, transitionFunction=mean, directions=directions)
    # geoCorrection for isotropic happens after transition matrices are combined
    # in first part of the TrackNLL.ObjFunc function
  } else if (formtype == "anisotropic"){
    coef <- X[[4]]
    transition.matrix <- transition(covar, transitionFunction=function(x){x[2] - x[1]}, directions=directions, symm=FALSE) # calculate a transition matrix
    delt <- geoCorrection(transition.matrix, type="c") 
    adj <- adjacent(covar, cells=1:ncell(covar), pairs=TRUE,  directions=directions)
    speed <- delt
    speed[adj] <-  exp(eval(parse(text=form))) # note that "form" is "(abs(delt[adj])*coef)"
    tr1 <- speed
  } 
}

GenerateL_t <- function(pars, onetrack, model.list, param_factor){
  pars <- pars*param_factor
  lambda<-pars[1] 
  for(i in 2:length(pars)){
    model.list[[i-1]][[4]]<-pars[i]
  }
  tm.corr <- ModelListToTransMatListToCombinedTransMat(model.list=model.list, directions=directions)
  norm.corr <- CalcNormCorr(X=onetrack, trans.mat=tm.corr)
  P_t <- NormcorrListToProbSurf.lap(X=norm.corr, rate=lambda)
  L_t <- P_t/(sum(P_t[]))
  return(L_t)
}

PlotOneTrack <- function(pars, onetrack, model.list, param_factor){
  pars <- pars*param_factor
  lambda<-pars[1] 
  for(i in 2:length(pars)){
    model.list[[i-1]][[4]]<-pars[i]
  }
  tm.corr <- ModelListToTransMatListToCombinedTransMat(model.list=model.list, directions=directions)
  norm.corr <- CalcNormCorr(X=onetrack, trans.mat=tm.corr)
  P_t <- NormcorrListToProbSurf.lap(X=norm.corr, rate=lambda)
  L_t <- P_t/(sum(P_t[]))
  onetrack <- onetrack[,"trackid"]
  names(onetrack) <- "id"
  l <- TrackSPDFtoTrackSPLFLines(onetrack) 
  Lcorr.ras <- ggspatial::df_spatial(L_t)
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
  L_t.plot.simp
}

LLProfileTrack <- function(onetrack=onetrack,
                           fittedpars=NULL,
                         model.list = focal.model.list, 
                         variables = variables,
                         formulas = model.split1,
                         plotme=0, 
                         param_factor=param_factor, 
                         params_nll_out_file=params_nll_out_file,
                         folderpath = folderpath,
                         searchgrid = NULL #exp(seq(-5,5,.25))
){
  npars <- length(fittedpars)
  variablenames <- c("loglambda", formulas)
  # generate estimate nll
  pars.list <- list(fittedpars)
  outprofile.fitted <- lapply(pars.list, TrackNLL.ObjFunc, 
                                #lambdalog = TRUE,       
                                model.list = model.list, 
                                onetrack=onetrack,
                                param_factor=param_factor)   
  fitted.nll <-unlist(outprofile.fitted)
  all.higher.ci <- c()
  all.lower.ci <- c()
  all.minpar <- c()
  all.minnll <- c()
  paramsdf.list.out <- list()
  grid <- searchgrid
  np <- 1
  for (np in 1:npars){ # loop through variables
    higher.ci <- NULL
    lower.ci <- NULL
    minpar <- NULL
    min.nll <- NULL
    foc.par <- fittedpars[np]
    #evaluate lower part of profile
    avec <- c(foc.par, foc.par-grid) # set up search vector
    if(foc.par > 0){ # if focal par is greater than 0, for lower wing set up search grid dense near 0
      avec <- c(avec, grid[grid < foc.par])
    }
    if(np == 1){
      avec[avec>0]<- -avec[avec>0]
    }
    avec <- sort(avec)  
    grid.params.df <- matrix(rep(fittedpars, each=length(avec)), nrow=length(avec)) # set up search dataframe
    grid.params.df[,np]<-avec
    pars.list <- split(grid.params.df, seq(nrow(grid.params.df)))
    
    outprofile <- lapply(pars.list, TrackNLL.ObjFunc, 
                                #lambdalog = TRUE,       
                                model.list = model.list, 
                                onetrack=onetrack,
                                param_factor=param_factor) 

    nlls.avec1<-as.numeric(unlist(outprofile)) # extract the nlls 
    print(nlls.avec1)
    
    grid.params.df <- cbind(grid.params.df,nlls.avec1)
    
    grid.params.df.all <- grid.params.df
    grid.params.df.lower <- grid.params.df.all

    low.val <- grid.params.df.lower[which.min(grid.params.df.lower[,npars+1]),np]
    prof.lower <-  grid.params.df.lower[grid.params.df.lower[,np]<low.val,npars+1]
    prof.avec <- grid.params.df.lower[grid.params.df.lower[,np]<low.val,np]
    lower.ci <- tryCatch({approx(prof.lower, prof.avec, xout=min(nlls.avec1)+qchisq(0.95,1)/2)}, # change second parameter of qchisq() for bivariate or multivariate CIs
                         error = function(e){data.frame(x=NA, y=NA)})
    print("lower")
    print(lower.ci)
    avec <- c(foc.par, foc.par+grid)   
    if(foc.par < 0){
      avec <- c(avec, -grid[(-grid)>foc.par]) # needs checking
    }
    # this forces loglambda (np==1) values to stay below 0 
    if(np == 1){
      avec[avec>0]<- -avec[avec>0]
    }
    avec <- sort(avec)    
    grid.params.df <- matrix(rep(fittedpars, each=length(avec)), nrow=length(avec))
    grid.params.df[,np]<-avec
    pars.list <- split(grid.params.df, seq(nrow(grid.params.df)))
    outprofile <- lapply(pars.list, TrackNLL.ObjFunc, 
                                model.list = model.list, 
                                onetrack=onetrack,
                                param_factor=param_factor) 
    nlls.avec1<-as.numeric(unlist(outprofile)) # extract the nlls 
    print(nlls.avec1)
    grid.params.df <- cbind(grid.params.df,nlls.avec1)
    grid.params.df.all <- grid.params.df
    grid.params.df.higher <- grid.params.df.all
    high.val <- grid.params.df.higher[which.min(grid.params.df.higher[,npars+1]),np]
    prof.higher <-  grid.params.df.higher[grid.params.df.higher[,np]>high.val,npars+1]
    prof.avec <- grid.params.df.higher[grid.params.df.higher[,np]>high.val,np]
    higher.ci <- tryCatch({approx(prof.higher, prof.avec, xout=min(nlls.avec1)+qchisq(0.95,1)/2)}, # change second parameter of qchisq() for bivariate or multivariate CIs
                          error = function(e){data.frame(x=NA, y=NA)})
    print("higher")
    print(higher.ci)
    grid.params.df.both <- rbind(grid.params.df.lower, grid.params.df.higher)
    grid.params.df.both <- data.frame(grid.params.df.both)
    min.nll <- min(grid.params.df.both$nlls.avec1, na.rm=T)[1]
    minpar.ind <- which.min(grid.params.df.both[,npars+1]) # troubleshooting
    minpar <- grid.params.df.both[minpar.ind,np]

    # plot all
    png(paste(folderpath,"/likprofiles/likprofile_",unique(onetrack$trackid),"_", variablenames[np],".png", sep=""),
        width=800, height=800, 
        units="px", pointsize =8,
        bg="transparent", res=100)
    par(mfrow=c(2,1))    
    plot(grid.params.df.both[,np], grid.params.df.both[,npars+1],
         xlab=paste("Beta_",variablenames[np], sep=""),
         ylab= "NLL") 
    min.nll.ind <- which.min(grid.params.df.both[,npars+1])
    points(grid.params.df.both[min.nll.ind,np], grid.params.df.both[min.nll.ind,npars+1], col="red", pch=16)
    points(fittedpars[np], fitted.nll, col="blue", pch="+", cex=2)
    abline(v=lower.ci$y, col="darkgreen")
    abline(v=higher.ci$y, col="darkblue")
    df <- data.frame(avec=grid.params.df.both[,np], 
                     nlls=grid.params.df.both[,npars+1])
    subsetdf.95 <-df[df$nlls < min(df$nlls, na.rm = T)+10,]
    
    plot(subsetdf.95$avec, subsetdf.95$nlls,ylab="NLL", xlab=paste("Beta_",variablenames[np], sep=""))  
    min.nll.ind <- which.min(subsetdf.95$nlls)
    points(subsetdf.95$avec[min.nll.ind], subsetdf.95$nlls[min.nll.ind], col="red", pch=16)
    points(fittedpars[np], fitted.nll, col="blue", pch="+", cex=2)
    abline(v=lower.ci$y, col="darkgreen")
    abline(v=higher.ci$y, col="darkblue")
    dev.off()
    
    l <- lower.ci$y
    h <- higher.ci$y
    
    # write to output
    all.lower.ci[np] <- l
    all.higher.ci[np] <- h
    all.minpar[np] <- minpar
    all.minnll[np] <- min.nll
  
    paramsdf.list.out[[np]] <- grid.params.df.both
    
    fittedpars[np] <- minpar
  } # end looping through parameters

  all.df <- do.call(rbind, paramsdf.list.out)
  globalmin.nll.ind <- which.min(all.df$nlls.avec1)
  global.minpar <- all.df[globalmin.nll.ind,1:npars]
  global.minpar <- as.numeric(global.minpar)
  global.minnll <- all.df$nlls.avec1[globalmin.nll.ind]
  
  outdf <- data.frame(variable = variablenames,
                      lower.ci = all.lower.ci, 
                      higher.ci = all.higher.ci, 
                      minpar = all.minpar, 
                      minnll = all.minnll,
                      global.minpar = global.minpar,
                      global.minnll = global.minnll)
  
  return(list(outdf, paramsdf.list.out))
  
}

LLProfileOptimizeTrack <- function (X,
                                startval.string,
                                focal.model.list,
                                param_factor,
                                cl.inside, 
                                folderpath, 
                                variables, 
                                model.split1, 
                                Lt.tiff = F) {
  onetrack <- X
  optim.start.par <- startval.string
  npars <- length(optim.start.par)
  
  repeatcounter <- 0
  repeat{  
    global.optim <- NULL
    
    repeatcounter <- repeatcounter + 1
    
    TrackNLL.ObjFunc(X=optim.start.par, 
                     onetrack=onetrack, 
                     model.list=focal.model.list,
                     param_factor=param_factor)
    
    global.optim <- optim(optim.start.par, TrackNLL.ObjFunc, 
                          method="Nelder-Mead",
                          model.list=focal.model.list, onetrack=onetrack, 
                          param_factor=param_factor,
                          hessian = T,
                          control = list(maxit=400))
    # calculate likelihood profile, see if can find a better nll than optim
    LL_list <- LLProfileTrack(fittedpars=global.optim$par,
                            model.list = focal.model.list, 
                            onetrack=onetrack, 
                            variables = variables,
                            formulas = model.split1,
                            param_factor=param_factor, 
                            folderpath = folderpath, 
                            searchgrid = exp(seq(-6,3,.5))) # third element of seq() adjusts grain, TODO make this a variable
    LLprofile.df <- do.call(rbind, LL_list[[2]])
    #write.csv(LLprofile.df, paste(folderpath, "/LLprofiledata_", unique(onetrack$trackid), "_", repeatcounter,".csv", sep=""))
    LL.CIs.df <-LL_list[[1]]
    LL.profile.nll <- min(LL.CIs.df$minnll, na.rm=T)
    
    globminind<-which.min(LLprofile.df$nlls.avec1)
    optim.start.par <- as.numeric(LLprofile.df[globminind,1:npars])
    if(LL.profile.nll > (global.optim$value - 0.1) | repeatcounter > 10){
      break
    }  
  }# end repeat loop  
  trackfit <- global.optim  
  vc <-  tryCatch({solve(trackfit$hessian)}, error = function(e){vc <- NA})
  SE <-   tryCatch({SE <- sqrt(diag(vc))}, error = function(e){SE <- NA} )
  
  lower.raw <- trackfit$par-1.96*SE
  upper.raw <- trackfit$par+1.96*SE
  
  parnames <- c("loglambda", variables)
  df.pars <- t(as.data.frame(trackfit$par, row.names = parnames))
  df.pars <- as.data.frame(df.pars)
  rownames(df.pars)<-NULL
  
  df.out<-data.frame(trackid=unique(onetrack$trackid), tracknll=trackfit$value, convergence=trackfit$convergence)
  df.out <- cbind(df.out, df.pars)
  df.out
  
  # alternate arrangement:
  parnames <- c("loglambda", variables)
  df.e <- as.data.frame(trackfit$par, row.names = parnames)
  df.l <- as.data.frame(lower.raw, row.names = parnames)
  df.u <- as.data.frame(upper.raw, row.names = parnames)
  df.SE <- as.data.frame(SE, row.names = parnames)
  df.c <-cbind(df.l, df.e, df.u, df.SE)
  names(df.c)<- c("lCI", "est", "uCI", "SE")
  df.c$vars <- parnames
  df.c$trackid <- unique(onetrack$trackid)
  df.c$tracknll <- trackfit$value
  trackplot <- PlotOneTrack(pars=as.numeric(df.pars), 
               onetrack=onetrack,
               model.list = focal.model.list,
               param_factor = param_factor)
  writeindivtracksdir=paste(folderpath, "indivtrackfits", sep="/")
   if(!dir.exists(writeindivtracksdir)){
      dir.create(writeindivtracksdir)}
    
    png(paste(writeindivtracksdir,"/trackplot_", unique(onetrack$trackid), ".png", sep=""),
        width=7, height=7, 
        units="in", pointsize =8,
        bg="transparent", res=300)
      print(trackplot)
    dev.off()
  
    if(Lt.tiff){ #Lt.tiff is a boolean
    L_t <- GenerateL_t(pars=as.numeric(df.pars), 
                onetrack=onetrack,
                model.list = focal.model.list,
                param_factor = param_factor)
    } else {
      L_t <- NULL
    }
    return(list(df.out, df.c, L_t))
} 

NormcorrListToProbSurf.lap <- function(X, rate){   
  norm.corr <- X[[2]]
  onetrack <- X[[1]]
  startpoint<- onetrack[1,]
  endpoint <- onetrack[nrow(onetrack),]
  corr.ys <- values(norm.corr)
  exp.prob <- dexp(corr.ys, rate=exp(rate))
  prob.corr <- norm.corr
  prob.corr[] <- exp.prob
  return(prob.corr)
}

TrackSPDFtoTrackSPLFLines <- function(tracks.spdf){
  x <- lapply(split(tracks.spdf, tracks.spdf$id), function(x) Lines(list(Line(coordinates(x))), x$id[1L]))
  lines <- SpatialLines(x)
  data <- data.frame(id = unique(tracks.spdf$id))
  rownames(data) <- data$id
  l <- SpatialLinesDataFrame(lines, data)
}
