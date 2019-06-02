#'@title PPM models fitted for GBIF occurrence data
#'
#'@authors: Payal Bal, Skipton Wooley
#'
#'@description Fit PPMs using GBIF data for IUCN listed species
#'GitHub: https://github.com/payalbal/gsdms/blob/master/R/6_ppms.R
#'
#'
#'

pacman::p_load(sp, raster, spatstat.data, nlme, rpart, spatstat, ppmlasso, 
       parallel, doParallel, doMC)

## Create 'output' folder
if(!dir.exists("./output")) {
  dir.create("./output")
}


## Load covariate data
covariates <- readRDS("./data/covariates.rds")
covariates_predict <- readRDS("./data/covariates_predict.rds")
raw_mask <- raster("./data/bio1.bil") 

## Load species data
add function...

gbif <- read.csv("./data/2019-05-14_gbif_iucnsp.csv", header = TRUE)


## Update mask based on NAs in covariate data
raw_mask[which(!is.na(raw_mask[]))] <- 1
source("./R/align.maskNA.R")
global_mask <- align.maskNA(covariates, raw_mask)
covariates <- mask(covariates, global_mask)
covariates_predict <- mask(covariates_predict, global_mask)


## Quadrature (background) points & associated covariate data
global_mask0 <- global_mask
global_mask0[which(is.na(global_mask0[]))] <- 0
rpts <- rasterToPoints(global_mask0, spatial=TRUE)
backxyz <- data.frame(rpts@data, X=coordinates(rpts)[,1], Y=coordinates(rpts)[,2])     
backxyz <- backxyz[,-1]
backxyz <- na.omit(cbind(backxyz, as.matrix(covariates)))
backxyz200k <- backxyz[sample(nrow(backxyz), 200000), ]
backxyz200k$Pres <- rep(0, dim(backxyz200k)[1])

  ## Checks
  summary(covariates)
  summary(backxyz200k)
  plot(global_mask0, legend = FALSE)
  plot(rasterFromXYZ(backxyz200k[,1:3]), col = "black", add = TRUE, legend=FALSE)
  

## Prediction points 
predxyz <- data.frame(rpts@data, X=coordinates(rpts)[,1], Y=coordinates(rpts)[,2])     
predxyz <- predxyz[,-1]
predxyz <- na.omit(cbind(predxyz, as.matrix(covariates_predict)))


## Define global ppm variables
  # the pred error you where getting was because there was different names between the fitted data and the predict data you had:
  ## "ppmxyz$X" not "X" so the model was just returning the fitted values from the ppmlasso model.
ppm_terms <- names(backxyz)[1:length(names(backxyz))-1]


## Specify covariates with interactions
interaction_terms <- c("X", "Y", "srtm")


## Fix n.fits = 20 and max.lambda = 100
lambdaseq <- round(sort(exp(seq(-10, log(100 + 1e-05), length.out = 20)), 
                        decreasing = TRUE),5)

## Initiate log files



## Estimate weights 
  ## calculate the area of the globe and then work out the weights based on the total area divided by number of points
ar <- raster::area(global_mask0)
ar <- mask(ar,global_mask)
totarea <- cellStats(ar,'sum')*1000 ## in meters^2
area_offset <- extract(ar, backxyz200k[,c('X','Y')], small = TRUE, fun = mean, na.rm = TRUE)*1000 ## in meters
bkgrd_wts <- c(totarea/area_offset)


  ## this was meant to calculate the species weights - but it's not quite right
  # spdat <- gbif
  # species_names <- levels(factor(spdat$species))
  # spwts <- list()
  # for(i in seq_along(species_names)){
  #       print(i)
  #       spxy <- spdat[spdat$species %in% species_names[i], c(4,3)]
  #       names(spxy) <- c("X", "Y")
  #       cellNo <- cellFromXY(ar,spxy)
  #       cellNoCounts <- table(cellNo)
  #       tmp_cell_area <- extract(ar, spxy, fun = mean, na.rm = TRUE)*1000
  #       tmp_dat <- data.frame(area=tmp_cell_area,cell_number=cellNo)
  #       tmp_wts <- ((tmp_dat$area*1000)/cellNoCounts[as.character(cellNo)])/totarea
  #       spwts[[i]] <- tmp_wts
  # }

    ## TESTING FOR CORRELATED VARIABLES
    ## There is no hard and fast rule about how many covariates to fit to data, and it will change depending on the data and the amount of information in each observation and how they vary with the covariates, how the covariates are correlated ect... But to start I'd only fit sqrt(n_observations) covariates. So if you have 20 occurrences that's 4-5 covariates (including the polynomials) so that's only two variables! You might have to identify the most important variable and start from there.
    
    ## let's do a PCA on the data to work out which are the most variable coefs.
    Xoriginal=t(as.matrix(backxyz200k))
    
    # Center the data so that the mean of each row is 0
    rowmns <- rowMeans(Xoriginal)
    X <-  Xoriginal - matrix(rep(rowmns, dim(Xoriginal)[2]), nrow=dim(Xoriginal)[1])
    
    # Calculate P
    A <- X %*% t(X)
    E <- eigen(A,TRUE)
    P <- t(E$vectors)
    
    dimnames(P) <- list(colnames(backxyz200k),paste0("PCA",1:ncol(P)))
    df <- as.data.frame(t(P[,1:5]))
    df$row.names<-rownames(df)
    library(reshape2)
    library(ggplot2)
    long.df<-reshape2::melt(df,id=c("row.names"))
    pca_plot <- ggplot2::ggplot(long.df,aes(x=row.names,y=variable,fill=value))+
      geom_raster()+
      scale_fill_viridis_c()+
      theme_minimal()
    
    ## not surprising that elevation, slope and aspect are all correlated (choose one).
    pca_plot


## Function: PPM model fitting        
fit_ppms_apply <- function(i, spdat, bkdat, bkwts, interaction_terms, ppm_terms, n.fits=50) {
   
  species_names <- levels(factor(spdat$species))
  cat('Fitting a ppm to', species_names[i],'\nThis is the', i,'^th model of',length(species_names),'\n')
  logfile <- paste0("./output/ppm_log_",gsub(" ","_",species_names[i]),"_",gsub("-", "", Sys.Date()), ".txt")
  writeLines(c(""), logfile)
  spxy <- spdat[spdat$species %in% species_names[i], c(4,3)]
  names(spxy) <- c("X", "Y")
  spxyz <- extract(covariates, spxy, buffer = 1000000, small = TRUE, 
                   fun = mean, na.rm = TRUE)
  
  spxyz <- cbind(spxy, spxyz)
  
  ## Remove remaining NA values & add check for < 20 records
  spxyz <- na.omit(spxyz)
  if (!(dim(spxyz)[1] < 20)) { 
    spxyz$Pres <- rep(1, dim(spxyz)[1])
    
    ## Specify weights for PPM data
    ppmxyz <- rbind(spxyz, bkdat)
    ppmxyz$wt <- NA
    ppmxyz[ppmxyz$Pres == 0,]$wt <- bkwts
    ppmxyz[ppmxyz$Pres == 1,]$wt <- 1e-6 #spwts[[i]]#
    
    ## Work out how many covariates to use in the model
    maxk <- 5 + length(ppm_terms[!(ppm_terms %in% interaction_terms)])*2
    nk <- ceiling(sqrt(nrow(spxy)))
    
    ## Specify ppm formula - independent with quadratic terms
    ## if low number of observations just use the covariates interaction terms for the spatial variables
    if(nk <= 5){   ppmform <- formula(paste0(" ~ poly(", paste0(interaction_terms, collapse = ", "),
                                             ", degree = 2, raw = FALSE)",collapse =""))
    } else  {
      extra_covar <- ceiling((nk - 5)/2)
      if(extra_covar > 8) extra_covar <- 8
      ppmform <- formula(paste0(paste0(" ~ poly(", paste0(interaction_terms, collapse = ", "),
                                       ", degree = 2, raw = FALSE)"), 
                                paste0(" + poly(", ppm_terms[!(ppm_terms %in% interaction_terms)][1:extra_covar],
                                       ", degree = 2, raw = FALSE)", 
                                       collapse = ""), collapse =""))
    }
    ## Fit ppm & save output
    cat(paste("\nFitting ppm model for",i , " @ ", Sys.time(), "\n"), 
        file = logfile, append = T)
    cat(paste("   # original records for",i , " = ", dim(spxy)[1], "\n"), 
        file = logfile, append = T)
    cat(paste("   # extracted records for",i , " = ", dim(spxyz)[1], "\n"), 
        file = logfile, append = T)
    ## data was being standardised twice!
    mod <- try(ppmlasso(formula = ppmform, data = ppmxyz, n.fits = n.fits, criterion = "bic",standardise = FALSE), silent=TRUE)
    gc()
    return(mod)
  } else {
    return(NULL)
  }
  
}

## Fit models and save output
spdat <- gbif
bkdat <- backxyz200k
bkwts <- bkgrd_wts
spp <- levels(factor(spdat$species))
mc.cores <- 1
seq_along(spp)
model_list <- parallel::mclapply(1:length(spp), fit_ppms_apply, spdat, #spwts,
                           bkdat, bkwts, interaction_terms, ppm_terms,
                           n.fits=100, mc.cores = mc.cores)

## Catch errors in model_list
errorfile <- paste0("./output/errorfile_", gsub("-", "", Sys.Date()), ".txt")
for(i in 1:length(model_list)){
  if(class(model_list[[i]]) == "try-error") {
    print(paste0("Model ",i, " for '", spp[i], "' has errors"))
    cat(paste(i, spp[i], "\n"),
        file = errorfile, append = T)
    model_list[[i]] <- NULL ## remove model from list
  }
}

saveRDS(model_list, file = paste0("./output/modlist_",  gsub("-", "", Sys.Date()), ".rds"))


## Function: PPM predictions
predict_ppms_apply <- function(i, models_list, newdata, bkdat, RCPs = c(26, 85)){
  
  cat('Predicting ', i,'^th model\n')
  
  if(class(models_list)== "try-error") {
    return(NULL)
  } else {
    
    predmu <- list()
    for (rcp in RCPs) {
      if (rcp == 26) {
        preddat <- newdata[which(names(newdata) %in% 
                                   names(newdata)[-grep('85', names(newdata))])]
        names(preddat) <- names(bkdat)[1:(length(names(bkdat))-2)]
        predmu$rcp26 <- predict.ppmlasso(models_list[[i]], newdata = preddat)
        # predmu <- 1-exp(-predmu) ## gives reative probabilities
        
      } else {
        preddat <- newdata[which(names(newdata) %in% 
                                   names(newdata)[-grep('26', names(newdata))])]
        names(preddat) <- names(bkdat)[1:(length(names(bkdat))-2)]
        predmu$rcp85 <- predict.ppmlasso(models_list[[i]], newdata = preddat)
        # predmu <- 1-exp(-predmu) ## gives reative probabilities
      }
      rm(preddat)
    }
    return(predmu)
  }
}

## Predict and save output
newdata <- predxyz
bkdat <- backxyz200k
RCPs <- c(26, 85)
mc.cores <- 1
prediction_list <- parallel::mclapply(seq_along(model_list), predict_ppms_apply,
                                      model_list, newdata, bkdat, RCPs, mc.cores = mc.cores)
saveRDS(prediction_list, file = paste0("./output/predlist_",  gsub("-", "", Sys.Date()), ".rds"))


## Prediction plots

for (i in 1:2){
  
  png(file = paste0("./output/plot_",tolower(gsub(" ","",spp[i])), ".png"),
      width=1500, height=600, res = 130, pointsize = 12)
  par(mfrow = c(1,2))
  par(mar=c(3,3,5,7))
  
  preds_sp <- rasterFromXYZ(cbind(predxyz[,1:2],prediction_list[[i]]))
  plot(preds_sp$rcp26, main = "rcp26", ext = extent(global_mask))
  points(spdat[spdat$species == spp[i], c(4,3)], pch=16, cex=.3)
  plot(preds_sp$rcp85, main = "rcp85", ext = extent(global_mask))
  points(spdat[spdat$species == spp[i], c(4,3)], pch=16, cex=.3)
  mtext(spp[i], side = 3, line = -2, outer = TRUE, cex = 1.5, font = 2)
  rm(preds_sp)
  dev.off()
}
