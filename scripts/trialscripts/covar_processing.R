## Download covariate data for gsdms


## Set working environment ####
# devtools::install_github('skiptoniam/sense')
x <- c("sp", "raster", "rgdal", "sense", "tools", "bitops", "RCurl")
## For gdal to work on Mac in R, see: https://github.com/smwindecker/gdaltools

data_covs <- "/Volumes/discovery_data/gsdms_data" # ./data/..." # server


## Load covariate data ####

## WorldClim data
## Source: https://www.worldclim.org/version1
bio_current <- list.files(paste0(data_covs, "/bio_30s"), pattern = "bio_current*", full.names = TRUE)
bio_rcp26 <- list.files(file.path(data_covs, "bio_30s"), pattern = "*bc26*", full.names = TRUE)
bio_rcp85 <- list.files(file.path(data_covs, "bio_30s"), pattern = "*bc85*", full.names = TRUE)
bioclim <- c(bio_current, bio_rcp26, bio_rcp85)


## SRTM
## Source: https://topotools.cr.usgs.gov/gmted_viewer/gmted2010_global_grids.php
srtm <- file.path(data_covs, "srtm/mn30_grd/srtm.adf")


## Soil
## Source: https://daac.ornl.gov/SOILS/guides/igbp-surfaces.html
soil <- list.files(file.path(data_covs,"orders"), pattern = "*.dat", full.names = T, recursive = T)


## Landuse
## Source: https://code.earthengine.google.com/?scriptPath=users%2Fpayalbal%2Fglobal_layers%3Afraclu_mcloor_sklu
## LU classes: 
landuse <- list.files(file.path(data_covs, "copernicus/global_frac", ""))


 
## Create Mask from WorldClim layer ####
global_mask <- raster(bio_current[1])
global_mask[which(!is.na(global_mask[]))] <- 1
writeRaster(global_mask, filename = paste0(data_covs, "/globalmask.tif"))



## Data processing ####
## sense::gdal_crop
file_in <- c(bioclim, srtm, soil, landuse)
file_out <- paste0(data_covs, "/processed/", paste0(tools::file_path_sans_ext(basename(file_in)), "_treated.", tools::file_ext(file_in)))

e <- c(-180,180,-60,90)
reso <- res(global_mask)
mapply(gdal_crop, inpath = file_in, outpath = file_out, MoreArgs = list(extent=e,res=reso)) 


## SRTM variables
file_out <- list.files(data_gsdms_in, full.names = TRUE)
elevation <- raster(file_out[grep("srtm", file_out)])
aspect <- terrain(elevation, opt = "aspect")
slope <- terrain(elevation, opt = 'slope')
roughness <- terrain(elevation, opt = "roughness")


    # ## Landuse reclassification - lulcc analysis
    #   ## Original land use classes in GLCC data: https://www.usgs.gov/media/images/global-land-cover-characteristics-data-base-version-20
    #   ##    # 0	Water
    #   ##    # 1	Evergreen Needle leaf Forest
    #   ##    # 2	Evergreen Broadleaf Forest
    #   ##    # 3	Deciduous Needle leaf Forest
    #   ##    # 4	Deciduous Broadleaf Forest
    #   ##    # 5	Mixed Forests
    #   ##    # 6	Closed Shrublands
    #   ##    # 7	Open Shrublands
    #   ##    # 8	Woody Savannas
    #   ##    # 9	Savannas
    #   ##    # 10	Grasslands
    #   ##    # 11	Permanent Wetland
    #   ##    # 12	Croplands
    #   ##    # 13	Urban and Built-Up
    #   ##    # 14	Cropland/Natural Vegetation Mosaic
    #   ##    # 15	Snow and Ice
    #   ##    # 16	Barren or Sparsely Vegetated
    # landuse <- raster(file_out[grep("landuse", file_out)])
    # landuse[landuse[]%in%c(12,14)] <- 100
    # landuse[landuse[]%in%c(6,8,7,9,10)] <- 200
    # landuse[landuse[]%in%c(1,2,3,4,5)] <- 300
    # landuse[landuse[]%in%c(13)] <- 400
    # landuse[landuse[]%in%c(0, 11,15,16,17)] <- 500
    # landuse <- landuse/100 - 1 ## these steps ensure that substituted alues do not overlap wth existing values
    # # check
    # unique(values(landuse))


## Stack, mask & save as .rds ####
covariates_all <- stack(setdiff(file_out, file_out[grep("srtm|landuse", file_out)]), 
                        elevation, slope, roughness, aspect, landuse) 
## load and stack covariates files, excluding files for srtm and landuse
##  which were updated and therefore these rasters are loaded from memory.



rm(list=setdiff(ls(), c("data_gsdms", "covariates_all", "global_mask")))
crs(covariates_all) <- crs(global_mask)
covariates_all <- mask(covariates_all, global_mask)


## Update mask based on NAs in covariate data - to get min set NAs
source("./R/align.maskNA.R")

raw_mask <- global_mask
raw_mask[which(!is.na(raw_mask[]))] <- 1
global_mask <- align.maskNA(covariates_all, raw_mask)
covariates_all <- mask(covariates_all, global_mask)
names(covariates_all) <- sub("_treated","", names(covariates_all))

## Check to see same number of NAS in covariates as in global mask
summary(global_mask)
summary(covariates_all)

# saveRDS(covariates_all, file = paste0(data_gsdms, "/covariates_all.rds"))

## -----------------------------------------------------------------------------


## SUBSET COVARIATES
## Subset variables 
covariates_all <- readRDS(paste0(data_gsdms, "/covariates_all.rds"))
cov_keep <- c("bio1", "bio4","bio12", "bio15", "bulkdens","pawc","soilcarb","totaln",
              "srtm","slope","roughness","aspect","landuse")
covariates <- covariates_all[[which(names(covariates_all) %in% cov_keep)]]

## Test correlation in covariate
corr1 <- layerStats(covariates, stat = 'pearson', na.rm = TRUE)
cov_values <- getValues(covariates)
corr2 <- cor(cov_values, use = 'complete.obs', method = 'pearson') 
rm(cov_values)

## Visulaisation
# library(corrplot)
# library(mnormt); library(psych)
# library(reshape); library(GGally)
corrplot::corrplot(corr1$`pearson correlation coefficient`, type = "upper")
corrplot::corrplot(corr1$`pearson correlation coefficient`, type = "upper", method = "number")
corrplot::corrplot(corr2, type = "upper", method = "number")
psych::pairs.panels(corr1$`pearson correlation coefficient`, scale = TRUE)
# GGally::ggpairs(as.data.frame(corr1$`pearson correlation coefficient`)) # don't like this. Slow and affiliated to ggplot.


## Remove highly correlated covariates (> 0.8)
'%!in%' <- function(x,y)!('%in%'(x,y))
covariates <- covariates[[which(names(covariates) %!in% c("bio4", "totaln", "roughness"))]]  
## Bring landuse up in the stack
covariates <- subset(covariates, c(length(names(covariates)), 1:length(names(covariates))-1))
saveRDS(covariates, file = "./output/covariates.rds")



## FUTURE CLIMATE LAYERS
rm(cov_keep)
cov_keep <- names(covariates_all)[grep('26|85', names(covariates_all))]
cov_keep <- cov_keep[grep('701$|704$|7012$|7015$', cov_keep)]
cov_keep <- c(names(covariates), cov_keep)
covariates_predict <- covariates_all[[which(names(covariates_all) %in% cov_keep)]]
## Remove same variables as for covariates because prediction matrix and design matrix shoudl be the same for ppms
covariates_predict <- covariates_predict[[which(names(covariates_predict) %!in% c("bc26bi704" ,"bc85bi704","totaln", "roughness"))]]
## Bring landuse up in the stack
covariates_predict <- subset(covariates_predict, c(length(names(covariates_predict)), 1:(length(names(covariates_predict))-1)))
saveRDS(covariates_predict, file = "./output/covariates_predict.rds")






# -----------------------------------------------------------------------------

# ## OTHER 'SENSE' FUNCTIONS
# ## gdal_mask - NOT REQUIRED. BUT WORKS. 
# global_mask[which(is.na(global_mask[]))] <- 0   #gdal_mask requires mask with 0s for mask and 1s for keep 
# writeRaster(global_mask, paste0(data_processed, "/global_mask01.tif"))
# file_in <- list.files(data_gsdms_in, full.names = TRUE)
# file_out <- sub("_treated", "_v2", file_in)
# mapply(gdal_mask, inpath = file_in, outpath = file_out, MoreArgs = list(mask = paste0(data_processed, "/global_mask01.tif")))
# 
# ## gdal_resample - NOT REQUIRED, resampling is done in crop (but gdal_crop cannot specify menthod of resampling)
# gdal_resample(inpath = file_in[i], outpath = file_out[i], resolution = reso,
#               method = 'near', bigtif = FALSE, return_raster = FALSE)
# 
# # Downloading - DOESN'T WORK: large files get corrupted; is a generic issue
# download.file("http://biogeo.ucdavis.edu/data/worldclim/v2.0/tif/base/wc2.0_30s_bio.zip", destfile = "./wc2.0_30s_bio.zip")
# if (!dir.exists("wc2.0_30s_bio")) {
#   dir.create("wc2.0_30s_bio")
# }
# unzip(".wc2.0_30s_bio.zip", exdir = "./wc2.0_30s_bio")

## -----------------------------------------------------------------------------



## BIAS LAYER
## Distance to  roads, built up/population density, coastlines, pas, water
## gdal_distance()

## Roads
## source: ...
## .gdb converetd to .shp in ArcGIS
## Distance-to-feature raster in QGIS...
## Filter by ("EXS" = '0'  OR  '1')  AND ("FCLASS" = 0 OR 1 OR 2 OR 3 OR 4 OR 5)

## Protected areas
## source: ...
# pa <- shapefile(paste0(data_path, "raw/", "boundaries/WDPA_Nov2018-shapefile/WDPA_Nov2018-shapefile-polygons.shp"))

## Population density (2005)
## source: http://sedac.ciesin.columbia.edu/data/set/gpw-v4-population-density-rev10/data-download
# pop <- raster(paste0(data_path, "raw/", "population/NASA_SEDAC/gpw-v4-population-density-rev10_2005_30_sec_tif/gpw_v4_population_density_rev10_2005_30_sec.tif"))
# pop <- crop(pop, global_mask)
# names(pop) <- "pop"
#   ## NOT running extent(pop) <- extent(global_mask) because diff resolution gives pixels with unequal length and width.
# pop <- projectRaster(pop, global_mask)
# pop <- mask(pop, global_mask)
# pop <- round(pop, 2)
# writeRaster(pop, paste0(data_path, "processed/", "population.tif"), format = "GTiff"))
# rm(pop)
# gc()
# pop <- raster(paste0(data_path, "processed/","population.tif"))



## xxxx .. not working
## Alternatively, .. not working
## Get the ftp site for the data from: https://edcftp.cr.usgs.gov/project/glcc/globdoc2_0.html
## Go >> https://edcftp.cr.usgs.gov/project/glcc/globdoc2_0.html#anony
system("wget -r -np -nH --reject 'index.html*' -e robots=off https://edcftp.cr.usgs.gov/project/glcc/globe") #corrupted...
system("curl https://edcftp.cr.usgs.gov/project/glcc/globe/gigbp1_2.img.gz -o data/gigbp1_2.img.gz") ##corrupted...
R.utils::gunzip(file.path(data_covs, "gigbp1_2.img.gz"), destname = file.path(data_covs, "landcover.img"), remove = FALSE)
landuse <- file.path(data_covs, "gigbp1_2.img")
landuse_out <- file.path(".data/gigbp1_2.img")
gdal_crop(inpath = landuse, outpath = landuse_out, extent=e, res=reso)
