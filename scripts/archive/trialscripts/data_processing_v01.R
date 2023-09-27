## Set working environment ####
setwd("...")
# devtools::install_github('skiptoniam/sense')
# devtools::install_github('smwindecker/gdaltools')
## For gdal to work on Mac in R, see: https://github.com/smwindecker/gdaltools
x <- c("sp", "raster", "rgdal", "sense", "tools", "bitops", "RCurl")
lapply(x, require, character.only = TRUE)
data_raw <- "/Volumes/discovery_data/gsdms_data" # ./data/raw/..." # server
data_processed <- "/Volumes/discovery_data/gsdms_data/processed/"


## Load covariate data ####
## WorldClim data
## Source: https://www.worldclim.org/version1
bio_current <- list.files(paste0(data_raw, "/bio_30s"), pattern = "bio_current*", full.names = TRUE)
bio_rcp26 <- list.files(file.path(data_raw, "bio_30s"), pattern = "*bc26*", full.names = TRUE)
bio_rcp85 <- list.files(file.path(data_raw, "bio_30s"), pattern = "*bc85*", full.names = TRUE)
bioclim <- c(bio_current, bio_rcp26, bio_rcp85)

## SRTM
## Source: https://topotools.cr.usgs.gov/gmted_viewer/gmted2010_global_grids.php
srtm <- file.path(data_raw, "srtm/mn30_grd/srtm.adf")

## Soil
## Source: https://daac.ornl.gov/SOILS/guides/igbp-surfaces.html
soil <- list.files(file.path(data_raw,"orders"), pattern = "*.dat", full.names = T, recursive = T)

## Landuse
## Source: https://code.earthengine.google.com/?scriptPath=users%2Fpayalbal%2Fglobal_layers%3Afraclu_mcloor_sklu
## *** Here, copy Mar's stiched landuse layer (/tempdata/workdir/troubleshooting/outputs/lu_world.tif) and provide file path..
landuse <- ... 


## Create Mask from WorldClim layer ####
## *** Here reproject using the Equal Area projectionro
global_mask <- raster(bio_current[1])
global_mask[which(!is.na(global_mask[]))] <- 1
crs(global_mask) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0" 
  #ref: https://worldclim.org/data/v1.4/formats.html
writeRaster(global_mask, filename = file.path(data_processed))

gdalwarp() using reproj_ras() # see https://github.com/kapitzas/gdaltools/blob/master/R/reproj_ras.R


## Data processing ####
## sense::gdal_crop
file_in <- c(bioclim, srtm, soil, landuse)
file_out <- paste0(data_processed, paste0(tools::file_path_sans_ext(basename(file_in)), "_treated.", tools::file_ext(file_in)))
e <- c(-180,180,-60,90)
reso <- res(global_mask)
mapply(gdal_crop, inpath = file_in, outpath = file_out, MoreArgs = list(extent=e,res=reso)) 

## sense::gdal_mask
mapply(gdal_mask, inpath = file_out, outpath = file_out, MoreArgs = list(mask = file.path(data_raw, "globalmask.tif")))

## Find min non-NA set values across mask and covariates and sync NAs
for(j in 1:length(file_out)){
  print(paste0("processing file = ", j, " of ", length(file_out) ,":", file_out[j]))
  gdal_mask(inpath = file.path(data_raw, "globalmask.tif"), outpath = file.path(data_raw, "globalmask.tif"), MoreArgs = list(mask = file_out[j]))
}

for(j in 1:length(file_out)){
  print(paste0("processing file = ", j, " of ", length(file_out) ,":", file_out[j]))
  gdal_mask(inpath = file_out[j], outpath = file_out[j], MoreArgs = list(mask = file.path(data_raw, "globalmask.tif")))
}

  ## Check
  summary(file.path(data_processed, "globalmask.tif"))
  summary(stack(file_out))


## Check for correlations, subset and save as .rds ####
## Stack covariates
file_out <- list.files(data_processed, full.names = TRUE)
covariates <- stack(file_out[-grep("srtm", file_out)])
elevation <- raster(file_out[grep("srtm", file_out)])
aspect <- terrain(elevation, opt = "aspect")
slope <- terrain(elevation, opt = 'slope')
roughness <- terrain(elevation, opt = "roughness")
covariates <- stack(covariates, elevation, aspect, slope, roughness)

## Idenitfy correlated variables
cov_df <- getValues(covariates)
cov_df <- na.omit(cov_df)
preds <- colnames(correlations(cov_df))
saveRDS(preds, file = file.path(data_processed, "preds_aus.rds"))

## Remove highly correlated covariates (> 0.7)
'%!in%' <- function(x,y)!('%in%'(x,y))
covariates <- covariates[[which(names(covariates) %!in% preds)]]  
saveRDS(covariates, file = file.path(data_processed, "covariates.rds"))

## Covariates for model
cov_mod <- names(covariates)[-grep('bc26|bc85', names(covariates))]
cov_pred_rcp26 <- names(covariates)[-grep('bio_current|bc85', names(covariates))]
cov_pred_rcp85 <- names(covariates)[-grep('bio_current|bc26', names(covariates))]


## Create data frame for model building and prediction (2 scennarios rcp 26 and rcp 85) ####
covs <- covs_mod ## replace by the .rds you are interested in, i.e. for model building versus predictions
global_mask0 <- global_mask
global_mask0[which(is.na(global_mask0[]))] <- 0
rpts <- rasterToPoints(global_mask0, spatial=TRUE)
model_data <- data.frame(rpts@data, X=coordinates(rpts)[,1], Y=coordinates(rpts)[,2])     
model_data <- model_data[,-1]
model_data <- na.omit(cbind(model_data, as.matrix(covs)))


data <- as.data.frame(covs)
lu_names <- names(lu)
lu_cols <- grep(paste0(lu_names, collapse = "|"), colnames(data))
k <- length(lu_names)
data <- na.omit(data)
inds <- as.numeric(rownames(data))





