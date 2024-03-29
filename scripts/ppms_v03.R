## Example global PPMs

## Set up work environment ####
rm(list = ls())
gc()
# system("ps")
# system("pkill -f R")

x <- c('data.table', 'sp', 'raster', 
       'sense', 'tools', 'bitops', 'RCurl', 
       'rgdal', 'gdalUtils', 'usethis', 'rgeos',
       'spatstat.data', 'nlme', 'rpart', 'spatstat', 'ppmlasso', 
       'parallel', 'doParallel', 'doMC')
lapply(x, require, character.only = TRUE)
rm(x)

## Folder paths
data_dir <- "/home/payalb/gsdms_r_vol/tempdata/research-cifs/uom_data/gsdms_data"
input_dir <- file.path(data_dir, "layers_10k")
output_dir <- file.path(data_dir, "outputs_10k")
if(!dir.exists(output_dir)) {
  dir.create(output_dir)
}

## Functions
'%!in%' <- function(x,y)!('%in%'(x,y))



## Data ####
## >> Biodiversity data
gbif_wgs <- as.data.frame(fread(file.path(data_dir, "gbif/2019-05-14_gbif_iucnsp.csv")))
gbif_wgs[, .N, species]

## Project data points in equal area projection
wgs_crs  <-  "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
eqarea_crs <- "+proj=eqearth +datum=WGS84 +ellps=WGS84 +wktext" ## check conversion; looks finne in QGIS
spdf <- SpatialPointsDataFrame(coords = gbif_wgs[c("decimallongitude", "decimallatitude")],
                               data =gbif_wgs, proj4string = CRS(wgs_crs))
spdf <- spTransform(spdf, CRSobj = CRS(eqarea_crs))

gbif <- as.data.table(spdf)
gbif <- gbif[,c(6,7,2,1,5)]
names(gbif)[1:2] <- c("X", "Y")
fwrite(gbif, file = file.path(output_dir, "2019-05-14_gbif_iucnsp_EA.csv"))


## >> Model data: Quadrature (background) points ####
backxyz <- fread(file.path(input_dir, "covariates_model.csv"))
backxyz200k <- backxyz[sample(nrow(backxyz), 200000), ] ## ~1/5th
backxyz200k$Pres <- rep(0, dim(backxyz200k)[1])

## Checks
summary(backxyz200k)
global.mask <- raster(file.path(input_dir, "globalmask_10k_ee_minNA.tif"))

plot(global.mask, legend = FALSE)
plot(rasterFromXYZ(backxyz[,1:3]), col = "grey", add = TRUE, legend=FALSE)
plot(rasterFromXYZ(backxyz200k[,1:3]), col = "red", add = TRUE, legend=FALSE)
points(gbif[,1:2], col = "yellow")
  ## not sure why it doesn't plot properly in R, in QGIS is fine...


## >> Prediction data ####
## RCP 45 and 85
predxyz <- file.path(input_dir, "covariates_predict_rcp45.csv")
# predxyz <- file.path(input_dir, "covariats_predict_rcp85.csv")


## >> Define PPM parameters ####
## Define global ppm variables
ppm_terms <- names(backxyz200k) [1:(length(names(backxyz200k))-1)]

## Specify covariates with interactions
interaction_terms <- c("X", "Y")

## Fix n.fits = 20 and max.lambda = 100
lambdaseq <- round(sort(exp(seq(-10, log(100 + 1e-05), length.out = 20)), 
                        decreasing = TRUE),5)

## >> Estimate weights for background points ####
## Current code using raster::area() only works in WGS projection
## calculate the area of the globe and then work out the weights based on the total area divided by number of points
mask_file <- file.path(output_dir, "globalmask_10k_ee_minNA.tif")
global.mask0 <- global.mask <- raster(mask_file)
global.mask0[which(is.na(global.mask0[]))] <- 0
ar <- raster::area(global.mask)
ar <- mask(ar,global.mask)
totarea <- cellStats(ar,'sum')*10000 ## in meters^2
area_offset <- extract(ar, backxyz200k[,c('X','Y')], small = TRUE, fun = mean, na.rm = TRUE)*10000 ## in meters
bkgrd_wts <- c(totarea/area_offset)

## Estimate species weights - not quite right yet - SW
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


## >> Define model function ####
fit_ppms_apply <- function(i, spdat, bkdat, bkwts, interaction_terms, ppm_terms, species_names, n.fits=10, min.obs = 50) {
  
  ## Initialise log file
  cat('Fitting a ppm to', species_names[i],'\nThis is the', i,'^th model of',length(species_names),'\n')
  logfile <- paste0("./ppm_log_",gsub(" ","_",species_names[i]),"_",gsub("-", "", Sys.Date()), ".txt")
  writeLines(c(""), logfile)
  
  ## Definne species specific data & remove NA (if any)
  spxy <- spdat[spdat$species %in% species_names[i], c(1,2)]
  names(spxy) <- c("X", "Y")
  
  # ## Check if points fall outside the mask
  plot(global.mask)
  points(spxy)
  
  ## Move species occurrence points falling off the mask to nearest 'land' cells
  ## TO BE USED WITH CAUTION!!
  ## find points which fall in NA areas on the raster
  vals <- extract(global.mask, spxy)
  outside_mask <- is.na(vals)
  if(sum(outside_mask) > 0){
    outside_pts <- spxy[outside_mask, ]
    ## find the nearest land within 5 decimal degrees of these
    land <- nearestLand(outside_pts, global.mask, 1000000)
    ## replace points falling in NA with new points on nearest land
    spxy[outside_mask, ] <- land
    ## count how many were moved
    sum(!is.na(land[, 1]))
  }
  # ## PROBLEM: We lose data again i.e. number of unique locations is reduced. This can be problematic for ppms...
  # nrow(unique(outside_pts))
  # nrow(unique(land))
  
  
  ## Extract covariates for presence points
  
  ## WITH LAND USE AS A COVARIATE (IN COLUMN 1)
  ## For landuse: Take the non-NA value at shortest distance from point
  ## -- takes very logn to run; find alternative
  r = cov.mod[[1]] ## landuse raster
  spxy_landuse <- apply(X = spxy, MARGIN = 1, FUN = function(X) r@data@values[which.min(replace(distanceFromPoints(r, X), is.na(r), NA))])
  ## spxy_landuse <- extract(r, spxy, method='simple', na.rm=TRUE, factor = TRUE) still gives NAs...
  # ## Check: Compare if non-NA values match with extract()
  # spxy_landuse_extract <- extract(cov.mod[[1]], spxy, method = 'simple') ## gives NAs
  # spxy_landuse_extract[which(!is.na(spxy_landuse_extract))] == spxy_landuse[which(!is.na(spxy_landuse_extract))]

  ## For other covariates: extract() using buffer
  spxyz <- extract(cov.mod[[-1]], spxy, buffer = 1000000, small = TRUE,
                   fun = mean, na.rm = TRUE)
  spxyz <- cbind(spxy, spxy_landuse, spxyz) ## if land-use is used
  names(spxyz)[3] <- "landuse"
  
  
  # ## WITHOUT LAND USE
  # spxyz <- extract(cov.mod, spxy, buffer = 1000000, small = TRUE,
  #                  fun = mean, na.rm = TRUE)
  # spxyz <- cbind(spxy, spxyz)
  # spxyz <- na.omit(spxyz)
  
  ## Specify model based on number of observatios
  ## If number of observations < 20, do not fit a model.
  ## If number of observations >= 20, fit a model accordig to the followig rule:
  ##  1. If number of observations <= min.obs (default as 50), 
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
  min.obs <- 50
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

## >> Fit models ####
spdat <- as.data.frame(gbif)
bkdat <- backxyz200k
bkwts <- bkgrd_wts
spp <- levels(factor(spdat$species))[1:2]
mc.cores <- 1
seq_along(spp)
now <- Sys.time()
ppm_models <- parallel::mclapply(1:length(spp), fit_ppms_apply, spdat, #spwts,
                                 bkdat, bkwts, interaction_terms, ppm_terms, 
                                 species_names = spp, n.fits=100, min.obs = 50, mc.cores = mc.cores)
ppm_mod_time <- Sys.time()-now


## Catch errors in models & save outputs ####
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



## Prediction function ####
predict_ppms_apply <- function(i, models_list, newdata, bkdat, RCPs = c(26, 85)){
  
  cat('Predicting ', i,'^th model\n')
  
  if(class(models_list)== "try-error") { ## redundant because error models are removed
    return(NULL)
  } else {
    
    predmu <- list()
    
    ## Current predictions
    preddat <- newdata[which(names(newdata) %in% 
                               names(newdata)[-grep('45|85', names(newdata))])]
    predmu$current <- predict.ppmlasso(models_list[[i]], newdata = preddat)
    
    ## Future predictions
    for (rcp in RCPs) {
      if (rcp == 45) {
        preddat <- newdata[which(names(newdata) %in% 
                                   names(newdata)[-grep('85|bio', names(newdata))])]
        names(preddat) <- names(bkdat)[1:(length(names(bkdat))-1)]
        predmu$rcp26 <- predict.ppmlasso(models_list[[i]], newdata = preddat)
        # predmu <- 1-exp(-predmu) ## gives reative probabilities
        
      } else {
        preddat <- newdata[which(names(newdata) %in% 
                                   names(newdata)[-grep('45|bio', names(newdata))])]
        names(preddat) <- names(bkdat)[1:(length(names(bkdat))-1)]
        predmu$rcp85 <- predict.ppmlasso(models_list[[i]], newdata = preddat)
        # predmu <- 1-exp(-predmu) ## gives reative probabilities
      }
      rm(preddat)
    }
    return(predmu)
  }
}

## Predict and save output ####
newdata <- predxyz
bkdat <- backxyz200k
RCPs <- c(45, 85)
now <- Sys.time()
prediction_list <- parallel::mclapply(seq_along(model_list), predict_ppms_apply,
                                      model_list, newdata, bkdat, RCPs, mc.cores = mc.cores)
pred_time <- Sys.time()-now
saveRDS(prediction_list, file = paste0("./output/predlist_",  gsub("-", "", Sys.Date()), ".rds"))



## Locate errors and rerun analysis for species with errors ###
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


## Catch remianing errors ####
n <- 1
m <- 1
catch_errors(seq_along(error_models), ppm_models = model_list, species_names = spp, errorfile = errorfile)





## ------ EXTRAS ----------------

## Testing for corerelations in covariates ####
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
max.params <- length(names(cov.mod))
obs.cutoff <- n.obs[which(add.params > max.params)[1]]
add.params[which(add.params > max.params)] <- max.params
plot(n.obs,add.params, type ="l")
abline(v=obs.cutoff)
text(obs.cutoff - 20, 4, paste0("max obs for params cut off: ", obs.cutoff), srt=90)


## 
## >> Quadrature (background) points ####
cov.mod <- stack(raster(infile))
names(cov.mod) <- c("bio1", "srtm")

global.mask <- raster(mask_file)
global.mask[which(is.na(global.mask[]))] <- 0
rpts <- rasterToPoints(global.mask, spatial=TRUE)
backxyz <- data.frame(rpts@data, X=coordinates(rpts)[,1], Y=coordinates(rpts)[,2])     
## if NAs in mask
backxyz <- backxyz[,-1]
backxyz <- na.omit(cbind(backxyz, as.matrix(cov.mod)))
## if 0s in mask
# backxyz <- backxyz[backxyz$globalmask_ee_minNA == 1,]
# backxyz <- backxyz[,-1]
backxyz200k <- backxyz[sample(nrow(backxyz), 200000), ]
backxyz200k$Pres <- rep(0, dim(backxyz200k)[1])

