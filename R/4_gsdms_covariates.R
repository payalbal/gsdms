
## PREPARE SPATIAL COVARIATE DATA FOR GSDMs
##
## Notes:
## 1. For gdal to work on Mac in R, see: https://github.com/smwindecker/gdaltools
## 2. Create Bias layer using gdal_distance
## 3. Add direct download code for landuse, soil, srtm, (bias layers for later)



## SET WORKING ENVIRONMENT
library(sp)
library(raster)
library(rgdal)
library(sense) # devtools::install_github('skiptoniam/sense'): 
library(tools)
# library(rgeos)
# library(sf)
# library(dplyr)
# library(progress)

## Data paths
data_raw <- "/Volumes/discovery_data/data/raw"
data_gsdms_in <- "/Volumes/discovery_data/data/gsdms_data/covariates"
data_gsdms <- "/Volumes/discovery_data/data/gsdms_data"


## DOWNLOAD DATA
## WorldClim: http://www.worldclim.org/
  ## download by tile functionality is only required fro 30 arc sec resolution, which will need to be followed by tile merging 
  ## See SK's WorldClimTiles package
  ## WorldClim tilenames: tile_names <- c(paste0(0, 0:11), paste0(1, 0:11), paste0(2, 0:11), paste0(3, 0:11), paste0(4, 0:11))

  # name <- c("worldclim", "CMIP5")
  # var <- "Bio"
  # rcps <- c(26,85)
  # model <- "BC"
  # year = 70
  # path <- paste0(data_raw, "climate/")
  # 
  # if(is.null(path)) path <- getwd()
  # 
  # rs <- raster(nrows = 5, ncols = 12, xmn = -180, xmx = 180,
  #              ymn = -60, ymx = 90)
  # rs[] <- 1:length(rs)
  # points <- xyFromCell(rs, rs[])
  # tileout <- list()
  # 
  # for (i in 1:nrow(points)){
  #   if(tolower(name) == "worldclim"){
  #     tileout[[i]] <- getData(name = name,  var = var, res = 10, path = path, 
  #                             lon = points[i,1], lat = points[i,2])
  #     
  #   }else if(name == "CMIP5"){
  #     for(j in 1:length(rcps)){
  #       tileout[[i]] <- getData(name = name,  var = var, res = 10, rcp = rcps[j], 
  #                               year = year, model= model, path = path, 
  #                               lon = points[i,1], lat = points[i,2])
  #     }
  #   }
  # }
  # rm(list=setdiff(ls(), c("data_raw", "data_gsdms_in")))

bio_current <- list.files(paste0(data_raw, "/climate/wc10"), pattern = '.bil$', full.names = TRUE)
bio_future <- list.files(paste0(data_raw, "/climate/cmip5/10m"),full.names = TRUE)
bio_rcp26 <- bio_future[grepl(paste0("(?=.*26)(?=.*tif)"), bio_future, perl = TRUE)]
bio_rcp85 <- bio_future[grepl(paste0("(?=.*85)(?=.*tif)"), bio_future, perl = TRUE)]
bioclim <- c(bio_current, bio_rcp26, bio_rcp85)


## SRTM: http://srtm.csi.cgiar.org/srtmdata/ 
  ## SRTM tile names: tile_names <- as.vector(unlist(as.data.frame(outer(paste0(1:72, "_"), 1:24, FUN = "paste0"))))
  
  # name <- "srtm"
  # path <- paste0(data_raw, "srtm/")
  # if(is.null(path)) path <- getwd()
  # rs <- raster(nrows = 24, ncols = 72, xmn = -180, xmx = 180,
  #              ymn = -60, ymx = 60)
  # rs[] <- 1:length(rs)
  # points <- xyFromCell(rs, rs[])
  # tileout <- list()
  # for (i in 1:nrow(points)){
  #   # tileout[[i]] <- getData(name = toupper(name), path = path, lon = points[i,1], lat = points[i,2])
  #   tileout[[i]] <- SRTM(lon = points[i,1], lat = points[i,2])
  # }
srtm <- paste0(data_raw, "/srtm/NASA/sdat_10008_1_20181118_172630534_01d.tif") # source: https://webmap.ornl.gov/ogc/dataset.jsp?dg_id=10008_1


## Soil: https://daac.ornl.gov/SOILS/guides/igbp-surfaces.html
## Subset to "bulkdens","pawc","soilcarb","totaln":
  ## bulkdens - bulk density (weight of soil in a given volume)
  ## pawc - profile available water capacity (amount of water that can be stored in a soil profile)
  ## soilcarb - soil carbon density 
  ## totaln - total nitrogen density

## Add direct download
soil <- list.files(paste0(data_raw,"/soil/NASA/IGBPDIS_SURFPRODS_569/data"), pattern = "*.dat", full.names = T)[c(1,3,4,6)]


## Landuse: https://earthexplorer.usgs.gov/
## Add direct download
landuse <- paste0(data_raw,"/landcover/USGS_GLCC/glccgbe20_tif_geoproj/gbigbpgeo20.tif")



## DATA PROCESSING
## Create Mask from WorldClim layer
global_mask <- raster(paste0(data_raw, "/climate/wc10/bio1.bil")) # 10min resolution ~ 345km2
global_mask[which(!is.na(global_mask[]))] <- 1
# writeRaster(global_mask, filename = paste0(data_gsdms, "/globalmask10m.tif"))


## sense::gdal_crop
file_in <- c(bioclim, srtm, soil, landuse)
file_out <- paste0(data_gsdms_in, "/", paste0(tools::file_path_sans_ext(basename(file_in)), "_treated.", tools::file_ext(file_in)))

e <- c(-180,180,-60,90)
reso <- res(global_mask)
mapply(gdal_crop, inpath = file_in, outpath = file_out, MoreArgs = list(extent=e,res=reso)) 


## Rename files
srtmfile <- file_out[grep(file_path_sans_ext(basename(file_in[grep("srtm", file_in)])), file_out)]
file.rename(srtmfile, sub(file_path_sans_ext(basename(srtmfile)),"srtm_treated", srtmfile))
landusefile <- file_out[grep(file_path_sans_ext(basename(file_in[grep("landcover", file_in)])), file_out)]
file.rename(landusefile, sub(file_path_sans_ext(basename(landusefile)),"landuse_treated", landusefile))


## SRTM variables
file_out <- list.files(data_gsdms_in, full.names = TRUE)
elevation <- raster(file_out[grep("srtm", file_out)])
aspect <- terrain(elevation, opt = "aspect")
slope <- terrain(elevation, opt = 'slope')
roughness <- terrain(elevation, opt = "roughness")


## Landuse reclassification (as per SK's code)
## Original land use classes in GLCC data: https://www.usgs.gov/media/images/global-land-cover-characteristics-data-base-version-20
##    # 0	Water
##    # 1	Evergreen Needle leaf Forest
##    # 2	Evergreen Broadleaf Forest
##    # 3	Deciduous Needle leaf Forest
##    # 4	Deciduous Broadleaf Forest
##    # 5	Mixed Forests
##    # 6	Closed Shrublands
##    # 7	Open Shrublands
##    # 8	Woody Savannas
##    # 9	Savannas
##    # 10	Grasslands
##    # 11	Permanent Wetland
##    # 12	Croplands
##    # 13	Urban and Built-Up
##    # 14	Cropland/Natural Vegetation Mosaic
##    # 15	Snow and Ice
##    # 16	Barren or Sparsely Vegetated
landuse <- raster(file_out[grep("landuse", file_out)])
landuse[landuse[]%in%c(12,14)] <- 100
landuse[landuse[]%in%c(6,8,7,9,10)] <- 200
landuse[landuse[]%in%c(1,2,3,4,5)] <- 300
landuse[landuse[]%in%c(13)] <- 400
landuse[landuse[]%in%c(0, 11,15,16,17)] <- 500
landuse <- landuse/100 - 1 ## these steps ensure that substituted alues do not overlap wth existing values


## COMPARE COVARIATES FOR NAS AND REMOVE NAS
covariates_all <- stack(setdiff(file_out, file_out[grep("srtm|landuse", file_out)]), 
                        elevation, slope, roughness, aspect, landuse) 
  ## load and stack covariates files, excluding files for srtm and landuse
  ##  which were updated and therefore these rasters are loaded from memory.


## SAVE COVARIATES AS RDS
rm(list=setdiff(ls(), c("data_gsdms_in", "covariates_all", "global_mask")))
crs(covariates_all) <- crs(global_mask)
covariates_all <- mask(covariates_all, global_mask)
names(covariates_all) <- sub("_treated","", names(covariates_all))

  ## Check to see same number of NAS in covariates as in global mask
summary(global_mask)
summary(covariates_all)

# saveRDS(covariates_all, file = paste0(data_gsdms, "/covariates_all.rds"))

## -----------------------------------------------------------------------------


## SUBSET COVARIATES
## Subset  bvariables 
covariates_all <- readRDS(paste0(data_gsdms, "/covariates_all.rds"))
cov_keep <- names(covariates_all)[grep('bio1$|bio4$|bio12$|bio15$', names(covariates_all))]
cov_keep <- c(cov_keep, c("bulkdens","pawc","soilcarb","totaln",
                          "srtm","slope","roughness","aspect","landuse"))
covariates <- covariates_all[[which(names(covariates_all) %in% cov_keep)]]

## Test correlation in covariate
corr1 <- layerStats(covariates, stat = 'pearson', na.rm = TRUE)
cov_values <- getValues(covariates)
corr2 <- cor(cov_values, use = 'complete.obs', method = 'pearson') 
rm(cov_values)

## Visulaisation
library(corrplot)
library(mnormt); library(psych)
library(reshape); library(GGally)
corrplot::corrplot(corr1$`pearson correlation coefficient`, type = "upper")
corrplot::corrplot(corr1$`pearson correlation coefficient`, type = "upper", method = "number")
corrplot::corrplot(corr2, type = "upper", method = "number")
psych::pairs.panels(corr1$`pearson correlation coefficient`, scale = TRUE)
GGally::ggpairs(as.data.frame(corr1$`pearson correlation coefficient`)) # don't like this. Slow and affiliated to ggplot.


## Remove highly correlated covariates (> 0.8)
'%!in%' <- function(x,y)!('%in%'(x,y))
covariates <- covariates[[which(names(covariates) %!in% c("bio4", "totaln", "roughness"))]]  
# saveRDS(covariates, file = "./output/covariates.rds")



## FUTURE CLIMATE LAYERS
rm(cov_keep)
cov_keep <- names(covariates_all)[grep('26|85', names(covariates_all))]
cov_keep <- cov_keep[grep('701$|704$|7012$|7015$', cov_keep)]
cov_keep <- c(cov_keep, c("bulkdens","pawc","soilcarb","totaln",
                          "srtm","slope","roughness","aspect","landuse"))
covariates_predict <- covariates_all[[which(names(covariates_all) %in% cov_keep)]]
## Remove same variables as for covariates because prediction matrix and design matrix shoudl be the same for ppms
covariates_predict <- covariates_predict[[which(names(covariates_predict) %!in% c("bc26bi704" ,"bc85bi704","totaln", "roughness"))]]  
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
## gdal_distance...
  
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

