#'@title PPM models fitted for GBIF occurrence data
#'
#'@authors: Payal Bal, Skipton Wooley
#'
#'@description Fit PPMs using GBIF data for IUCN listed species
#'GitHub: https://github.com/payalbal/gsdms/blob/master/R/6_ppms.R



## SET UP WORKSPACE
pacman::p_load(sp, raster, spatstat.data, nlme, rpart, spatstat, ppmlasso, 
               parallel, doParallel, doMC)
# mc.cores <- 30 ## for boab
mc.cores <- 1

## Create 'output' folder
if(!dir.exists("./output")) {
  dir.create("./output")
}



## DATA
## Load covariate data
covariates <- readRDS("./output/covariates.rds")
covariates_predict <- readRDS("./output/covariates_predict.rds")
raw_mask <- raster("./data/bio1.bil") 

## Load species data
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

# ## Checks
# summary(covariates)
# summary(backxyz200k)
# plot(global_mask0, legend = FALSE)
# plot(rasterFromXYZ(backxyz200k[,1:3]), col = "black", add = TRUE, legend=FALSE)

## Prediction points 
predxyz <- data.frame(rpts@data, X=coordinates(rpts)[,1], Y=coordinates(rpts)[,2])     
predxyz <- predxyz[,-1]
predxyz <- na.omit(cbind(predxyz, as.matrix(covariates_predict)))



## MODEL PARAMETERS
## Define global ppm variables
ppm_terms <- names(backxyz200k)[1:(length(names(backxyz200k))-1)]

## Specify covariates with interactions
interaction_terms <- c("X", "Y")

## Fix n.fits = 20 and max.lambda = 100
lambdaseq <- round(sort(exp(seq(-10, log(100 + 1e-05), length.out = 20)), 
                        decreasing = TRUE),5)

## Estimate weights for background points
## calculate the area of the globe and then work out the weights based on the total area divided by number of points
ar <- raster::area(global_mask0)
ar <- mask(ar,global_mask)
totarea <- cellStats(ar,'sum')*1000 ## in meters^2
area_offset <- extract(ar, backxyz200k[,c('X','Y')], small = TRUE, fun = mean, na.rm = TRUE)*1000 ## in meters
bkgrd_wts <- c(totarea/area_offset)

## Estimate species weights - but it's not quite right
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



## SPECIFY MODEL
## Function: PPM model        
fit_ppms_apply <- function(i, spdat, bkdat, bkwts, interaction_terms, ppm_terms, species_names, n.fits=50, min.obs = 50) {
  
  ## Initialise log file
  cat('Fitting a ppm to', species_names[i],'\nThis is the', i,'^th model of',length(species_names),'\n')
  logfile <- paste0("./output/ppm_log_",gsub(" ","_",species_names[i]),"_",gsub("-", "", Sys.Date()), ".txt")
  writeLines(c(""), logfile)
  
  ## Definne species specific data & remove NA (if any)
  spxy <- spdat[spdat$species %in% species_names[i], c(4,3)]
  names(spxy) <- c("X", "Y")
  
  ## Extract covariates for presence points
  
  ## For landuse: Take the raster value with lowest distance to point AND non-NA value in the raster
  ## ---- takes very logn to run; fidn alternative ----
  r = covariates[[1]] ## landuse raster
  spxy_landuse <- apply(X = spxy, MARGIN = 1, FUN = function(X) r@data@values[which.min(replace(distanceFromPoints(r, X), is.na(r), NA))])
  ## spxy_landuse <- extract(r, spxy, method='simple', na.rm=TRUE, factor = TRUE) still gives NAs...
  # ## Check: Compare if non-NA values match with extract()
  # spxy_landuse_extract <- extract(covariates[[1]], spxy, method = 'simple') ## gives NAs
  # spxy_landuse_extract[which(!is.na(spxy_landuse_extract))] == spxy_landuse[which(!is.na(spxy_landuse_extract))]
  
  ## For other covariates: extract() using buffer
  spxyz <- extract(covariates[[-1]], spxy, buffer = 1000000, small = TRUE, 
                   fun = mean, na.rm = TRUE)
  
  spxyz <- cbind(spxy, spxy_landuse, spxyz)
  names(spxyz)[3] <- "landuse"
  spxyz <- na.omit(spxyz)
  
  ## Specify model based on number of observatios
  ## If number of observations < 20, do not fit a model.
  ## If number of observations >= 20, fit a model accordig to the followig rule:
  ##  1. If number of observations <= min.obs (deafult as 50), 
  ##    use only the spatial covariates i.e. lat and long. This gives 5 
  ##    terms: X, Y, X^2, Y^2, X:Y (linear, quadratic with interaction)
  ##  2. .... 50 - 90
  ##  2. If number of observatons > min.obs, 
  ##    use the one-in-ten rule (https://en.wikipedia.org/wiki/One_in_ten_rule)
  ##    where, an additional covariate is added for every 10 additonal obeservations.
  ##    Because we fit poly with degree 2, adding a covariate will add two 
  ##    terms x and x^2 to the model. Therefore, to operationalize this rule, 
  ##    we add 1 raw covariate for every 20 additional observations. All 10
  ##    covariates (in addition to lat long) are fit onlt at > 230 observations.
  ##    There are 12 covariates in all includign lat and log. 
  
  nk <- nrow(spxyz)
  if (!(nk < 20)) {
    
    ## Add presence column to species occurrence table
    spxyz$Pres <- rep(1, dim(spxyz)[1])
    
    ## Specify weights for PPM data
    ppmxyz <- rbind(spxyz, bkdat)
    ppmxyz$wt <- NA
    ppmxyz[ppmxyz$Pres == 0,]$wt <- bkwts
    ppmxyz[ppmxyz$Pres == 1,]$wt <- 1e-6 #spwts[[i]]
    
    ## Specify PPM formula based on number of observations
    if(nk <= min.obs){
      ppmform <- formula(paste0(" ~ poly(", paste0(interaction_terms, collapse = ", "),
                                ", degree = 2, raw = FALSE)",collapse =""))
    } else {
      if(nk > min.obs & nk <= min.obs + 40){
        ppmform <- formula(paste0(" ~ poly(", paste0(interaction_terms, collapse = ", "),
                                  ", degree = 2, raw = FALSE) + factor(landuse)",collapse ="")) 
      } else {
        extra_covar <- ceiling((nk - min.obs)/20)
        if(extra_covar > 10) extra_covar <- 10 ## Fit all 10 covariates
        ppmform <- formula(paste0(paste0(" ~ poly(", paste0(interaction_terms, collapse = ", "),
                                         ", degree = 2, raw = FALSE) + factor(landuse)"), 
                                  paste0(" + poly(", ppm_terms[!(ppm_terms %in% c(interaction_terms, "landuse"))][1:extra_covar],", degree = 2, raw = FALSE)", collapse = ""), collapse =""))
      }
    }
    
    ## Fit ppm & save output
    cat(paste("\nFitting ppm model for",i , " @ ", Sys.time(), "\n"), 
        file = logfile, append = T)
    cat(paste("   # original records for",i , " = ", dim(spxy)[1], "\n"), 
        file = logfile, append = T)
    cat(paste("   # extracted records for",i , " = ", dim(spxyz)[1], "\n"), 
        file = logfile, append = T)
    
    mod <- try(ppmlasso(formula = ppmform, data = ppmxyz, n.fits = n.fits, 
                        criterion = "bic", standardise = FALSE), silent=TRUE)
    # cat('Warnings for ', species_names[i],':\n')
    # warnings()
    
    # ## capture messages and errors to a file
    # sink(logfile, type = "message", append = TRUE, split = FALSE)
    # try(warnings())
    # ## reset message sink and close the file connection
    # sink(type="message")
    # close(logfile)
    
    gc()
    return(mod)
  } else {
    return(NULL)
  }
}

## Fit models
spdat <- gbif
bkdat <- backxyz200k
bkwts <- bkgrd_wts
spp <- levels(factor(spdat$species))
seq_along(spp)
ppm_models <- parallel::mclapply(1:length(spp), fit_ppms_apply, spdat, #spwts,
                                 bkdat, bkwts, interaction_terms, ppm_terms, 
                                 species_names = spp, n.fits=100, min.obs = 50, mc.cores = mc.cores)

## Catch errors in models & save outputs
error_list <- list()
model_list <- list()
n <- 1
m <- 1
errorfile <- paste0("./output/errorfile_1_", gsub("-", "", Sys.Date()), ".txt")

for(i in 1:length(ppm_models)){
  if(!class(ppm_models[[i]])[1] == "try-error") {
    model_list[[n]] <- ppm_models[[i]]
    n <- n+1
  }else{
    print(paste0("Model ",i, " for '", spp[i], "' has errors"))
    cat(paste(i, ",", spp[i], "\n"),
        file = errorfile, append = T)
    error_list[[m]] <- ppm_models[[i]]
    m <- m+1
  }
}

saveRDS(model_list, file = paste0("./output/modlist_",  gsub("-", "", Sys.Date()), ".rds"))
saveRDS(error_list, file = paste0("./output/errlist_",  gsub("-", "", Sys.Date()), ".rds"))



## PREDICTION
## Function: PPM predictions
predict_ppms_apply <- function(i, models_list, newdata, bkdat, RCPs = c(26, 85)){
  
  cat('Predicting ', i,'^th model\n')
  
  if(class(models_list)== "try-error") { ## redundant because error models are removed
    return(NULL)
  } else {
    
    predmu <- list()
    
    ## Current predictions
    preddat <- newdata[which(names(newdata) %in% 
                               names(newdata)[-grep('26|85', names(newdata))])]
    predmu$current <- predict.ppmlasso(models_list[[i]], newdata = preddat)
    
    ## Future predictions
    for (rcp in RCPs) {
      if (rcp == 26) {
        preddat <- newdata[which(names(newdata) %in% 
                                   names(newdata)[-grep('85|bio', names(newdata))])]
        names(preddat) <- names(bkdat)[1:(length(names(bkdat))-1)]
        predmu$rcp26 <- predict.ppmlasso(models_list[[i]], newdata = preddat)
        # predmu <- 1-exp(-predmu) ## gives reative probabilities
        
      } else {
        preddat <- newdata[which(names(newdata) %in% 
                                   names(newdata)[-grep('26|bio', names(newdata))])]
        names(preddat) <- names(bkdat)[1:(length(names(bkdat))-1)]
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
prediction_list <- parallel::mclapply(seq_along(model_list), predict_ppms_apply,
                                      model_list, newdata, bkdat, RCPs, mc.cores = mc.cores)
saveRDS(prediction_list, file = paste0("./output/predlist_",  gsub("-", "", Sys.Date()), ".rds"))



## LOCATE ERRORS AND RERUN ANALYSIS FOR SPECIES WITH ERROR
##  To be automated if error problem is not solved by species grouping
##  At the moment, it appears that error might be when species data is spatially restricted.
error_species <- read.table("./output/errorfile_1_20190725.txt", header = FALSE, sep = ",")
colnames(error_species) <- c("index", "species")
error_index <- error_species$index

names(model_list) <- tolower(gsub(" ","_", levels(factor(spdat$species))[-error_index]))
names(prediction_list) <- tolower(gsub(" ","_", levels(factor(spdat$species))[-error_index]))

spdat <- gbif
bkdat <- backxyz200k
bkwts <- bkgrd_wts
spp <- levels(factor(spdat$species))[error_index]
seq_along(spp)
mc.cores <- 1
error_models <- parallel::mclapply(1:length(spp), fit_ppms_apply, spdat, #spwts,
                                   bkdat, bkwts, interaction_terms, ppm_terms,
                                   species_names = spp, n.fits=100, min.obs = 50, mc.cores = mc.cores)
names(error_models) <- tolower(gsub(" ","_", levels(factor(spdat$species))[error_index]))

newdata <- predxyz
bkdat <- backxyz200k
RCPs <- c(26, 85)
error_pred <- parallel::mclapply(seq_along(error_models), predict_ppms_apply,
                                 error_models, newdata, bkdat, RCPs, mc.cores = mc.cores)
names(error_pred) <- tolower(gsub(" ","_", levels(factor(spdat$species))[error_index]))

errorfile <- paste0("./output/errorfile_2_", gsub("-", "", Sys.Date()), ".txt")
n <- length(model_list)+1
m <- length(errlist)+1
error_list <- list()
for (i in 1:length(error_models)){
  if(!class(error_models[[i]])[1] == "try-error") {
    model_list[[n]] <- error_models[[i]]
    n <- n+1
  }else{
    print(paste0("Model ",i, " for '", spp[i], "' has errors"))
    cat(paste(i, ",", spp[i], "\n"),
        file = errorfile, append = T)
    error_list[[m]] <- error_models[[i]]
    m <- m+1
  }
}
names(model_list)[((length(model_list)-length(error_models))+1):length(model_list)] <- names(error_models)

n <- length(prediction_list) + 1
for (i in 1:length(error_pred)) {
  prediction_list[n] <- error_pred[i]  
  n <- n + 1
}
names(prediction_list)[((length(prediction_list)-length(error_pred))+1):length(prediction_list)] <- names(error_pred)






## -----------------------------------------------------------------------------
## EXTRAS
## -----------------------------------------------------------------------------

## FUNCTION: Catch errors in model outputs
##  (1) List of outputs for models without errors, 
##  (2) list of outputs for models with errors, 
##  (3) text file with list of models with errors (species indices and species names) 
catch_errors <- function(i, ppm_models, species_names, errorfile) {
  cat('Checking model for ', species_names[i],'for errors\n')
  for(i in 1:length(ppm_models)){
    if(!class(ppm_models[[i]])[1] == "try-error") {
      model_list[[n]] <- ppm_models[[i]]
      n <- n+1
    }else{
      print(paste0("Model ",i, " for '", spp[i], "' has errors"))
      cat(paste(i, ",", spp[i], "\n"),
          file = errorfile, append = T)
      error_list[[m]] <- ppm_models[[i]]
      m <- m+1
    }
  }
}
catch_errors(seq_along(error_models), ppm_models = model_list, species_names = spp, errorfile = errorfile)


## TESTING FOR CORRELATED VARIABLES
## There is no hard and fast rule about how many covariates to fit to data, and it will change depending on the data and the amount of information in each observation and how they vary with the covariates, how the covariates are correlated ect... But to start I'd only fit sqrt(n_observations) covariates. So if you have 20 occurrences that's 4-5 covariates (including the polynomials) so that's only two variables! You might have to identify the most important variable and start from there.

## see slide 14 & 26 in: http://www.bo.astro.it/~school/school09/Presentations/Bertinoro09_Jasper_Wall_3.pdf

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


##ALTERNATIVE - Specify ppm formula based on number of observations for a species
nk <- ceiling(sqrt(nrow(spxy)))
## Specify ppm formula - independent with quadratic terms
## if low number of observations, (i.e. =< 25), just use the covariates interaction terms for the spatial variables
if(nk <= 5){ ## because linear and quad for X ad Y gives 5 terms
  ppmform <- formula(paste0(" ~ poly(", paste0(interaction_terms, collapse = ", "),
                            ", degree = 2, raw = FALSE)",collapse =""))
  ## including land-use change for species with low number of observations
  # paste0(" ~ poly(", paste0(interaction_terms, collapse = ", "),
  #        ", degree = 2, raw = FALSE) + poly(landuse, degree = 2, raw = FALSE)",collapse ="")
  
} else  {
  ## if number of observations, (i.e. > 25), fit independent with linear and quadratic terms startign with landuse
  extra_covar <- ceiling((nk - 5)/2) ## fit additional covariates based on 
  if(extra_covar > 10) extra_covar <- 10 ## if nrow(spxy) > 630 and nk > 25, only then all parameters are used!! ,...fix
  ppmform <- formula(paste0(paste0(" ~ poly(", paste0(interaction_terms, collapse = ", "),
                                   ", degree = 2, raw = FALSE)"), 
                            paste0(" + poly(", ppm_terms[!(ppm_terms %in% interaction_terms)][1:extra_covar],
                                   ", degree = 2, raw = FALSE)", 
                                   collapse = ""), collapse =""))
}

## Explore relationship between a method to specify cut-off and 3obs
n.obs <- seq(20,1000, 10)
cutoff <- ceiling(sqrt(n.obs))
add.params <- ceiling((sqrt(n.obs) - 5)/2)
max.params <- length(names(covariates))
obs.cutoff <- n.obs[which(add.params > max.params)[1]]
add.params[which(add.params > max.params)] <- max.params
plot(n.obs,add.params, type ="l")
abline(v=obs.cutoff)
text(obs.cutoff - 20, 4, paste0("max obs for params cut off: ", obs.cutoff), srt=90)
