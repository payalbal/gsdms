## MBH design - Scott D.Foster

library(sp)
library(raster)
# install.packages( "MBHdesign")
library(MBHdesign)

data_raw <- "/Volumes/discovery_data/data/raw"
data_gsdms <- "/Volumes/discovery_data/data/gsdms"


## Load mask
global_mask <- raster(paste0(data_raw, "/climate/wc10/bio1.bil"))
global_mask[which(!is.na(global_mask[]))] <- 1
global_mask[which(is.na(global_mask[]))] <- NA


## Define sampling grid
## Get latitude and longitude values for raster
rpts <- rasterToPoints(global_mask, spatial=TRUE)
X <- data.frame(rpts@data, long=coordinates(rpts)[,1],
                lat=coordinates(rpts)[,2])     
colnames(X)[1] <- "mask"
## Get covariatres values
cov_values <- as.matrix(covariates)
  ## Checks
  summary(cov_values)
  summary(covariates)

X <- cbind(X, complete.cases(cov_values))
  ## Check - that 0 mask values align with NA values in covariates: yes
X <- X[,-1]
rm(cov_values, rpts)


## Define inclusion probabilities as 0 for NA and 1 for non-NA values in global mask
global_mask[which(is.na(global_mask[]))] <- 0
InclProb <- as.matrix(global_mask)
## Checks: 
sum(is.na(values(global_mask))) # to check for NA values.
print(paste0("# ones in global mask = ", sum(values(global_mask) == 1)))
print(paste0("# zeros in global mask = ", sum(values(global_mask) == 0)))
print(paste0("check that total # cells in global mask is equal to the sum of ones and zeros : ",
             length(global_mask@data@values) == sum(values(global_mask) == 1)
             + sum(values(global_mask) == 0)))
  
## Qausi random sample ... Gives error (see details below). FIX LATER
  # library(rgdal)
  # ## Number of samples
  # n <- 200000 # Renner et al 2015 
  # ## Number of points to sample from (= ncells in rasters)
  # N <- dim(X)[1]
  # 
  # samp <- quasiSamp(n=n, dimension = dim(X)[2], potential.sites=X, inclusion.probs=InclProb)
  # head(samp, row.names=FALSE)
  # 
  # ## visualise
  # as.raster(samp...)
  # plot(inclprobraster)
  # points(samp[,1:2], pch=21, bg=grey(0.75), cex=1.5)
  # 
  # 
  ## Errors with MBH:
  # When running the quasiSamp with potential.sites just as a lat long matrix, I'm told:
  # Error in quasiSamp(n = n, dimension = dim(X)[2], potential.sites = X,  : 
  #   Failed to find a design. It is possible that the inclusion probabilities are 
  #   very low and uneven OR that the sampling area is very irregular (e.g. long and skinny) 
  #   OR something else. Please try again (less likely to work) OR make inclusion probabilities 
  #   more even (more likely but possibly undesireable) OR increase the number of sites considered 
  #   (likely but computationally expensive).
  #
  # If I give the potential matrix as a complete matrix with all covariates in 
  # addition to lat-long, this is the error i get: 
  # Error in class::knn1(potential.sites, samp[, 1:dimension], 1:nrow(potential.sites)) : 
  #   no missing values are allowed


## Random sample - FOR NOW
Xsamp <- X[sample(nrow(X), 100000), ]
## OR
# raster::sampleRandom(global_mask, 100000)
summary(X)
summary(Xsamp)
# prop of NAs similar in X and Xsamp ~ 0.3


# saveRDS(Xsamp, file = "./output/designmatrix.rds")





