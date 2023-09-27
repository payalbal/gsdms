
## GLOBAL - STICKY TAPE ANALYSIS @ 1km2 resolution
## -------------------------------------------------

## Prepare spatial predictor layers for use in landuse models and SDMs
## Processing done on boab and resulting rasters saved.

## Note for gdal to work on Mac in R, see: https://github.com/smwindecker/gdaltools

library(sp)
library(raster)
library(rgdal)
library(rgeos)
library(sf)
library(dplyr)
library(progress)
# devtools::install_github('skiptoniam/sense')
library(sense)

# library(pacman)
# library(data.table,sp,raster,rgdal,rgeos,sf,dplyr,progress, sense)

## Data paths
data_raw <- "/Volumes/payal_umelb/data/raw/"
data_processed <- "/Volumes/payal_umelb/data/processed/"

## DOWNLOAD DATA

## WorldClim: http://www.worldclim.org/
  ## download by tile functionality is only required fro 30 arc sec resolution, which will need to be followed by tile merging 
  ## See SK's WorldClimTiles package
  ## WorldClim tilenames: tile_names <- c(paste0(0, 0:11), paste0(1, 0:11), paste0(2, 0:11), paste0(3, 0:11), paste0(4, 0:11))

name <- c("worldclim", "CMIP5")
var <- "Bio"
rcps <- c(26,85)
model <- "BC"
year = 70
path <- paste0(data_raw, "climate/")

if(is.null(path)) path <- getwd()

rs <- raster(nrows = 5, ncols = 12, xmn = -180, xmx = 180,
             ymn = -60, ymx = 90)
rs[] <- 1:length(rs)
points <- xyFromCell(rs, rs[])
tileout <- list()

for (i in 1:nrow(points)){
  if(tolower(name) == "worldclim"){
    tileout[[i]] <- getData(name = name,  var = var, res = 10, path = path, 
                            lon = points[i,1], lat = points[i,2])
    
  }else if(name == "CMIP5"){
    for(j in 1:length(rcps)){
      tileout[[i]] <- getData(name = name,  var = var, res = 10, rcp = rcps[j], 
                              year = year, model= model, path = path, 
                              lon = points[i,1], lat = points[i,2])
    }
  }
}
rm(rs,points,tileout)

## MASK
## source: DEM.tif from IUCN basedata: https://www.iucnredlist.org/resources/spatialtoolsanddata
  # global_mask <- raster(paste0(data_path, "raw/", "biodiversity/IUCN/basedata/DEM/DEM.tif"))
  # global_mask[which(!is.na(global_mask[]))] <- 0
  # e <- c(-180,180,-60,90)
  # global_mask <- crop(global_mask, e)
  # writeRaster(global_mask, paste0(data_path, "processed/", "global_mask"), format = "GTiff", options="COMPRESS=LZW", overwrite = TRUE)
global_mask <- raster(paste0(data_path, "processed/", "global_mask.tif"))

## mask with 1 values
  # global_mask <- raster("DEM.tif")
  # global_mask[which(!is.na(global_mask[]))] <- 1
  # e <- c(-180,180,-60,90)
  # global_mask <- crop(global_mask, e) # to remove Antarctica
  # writeRaster(global_mask, "global_mask1", format = "GTiff", options="COMPRESS=LZW", overwrite = TRUE)
global_mask1 <- raster(paste0(data_path, "processed/", "global_mask1.tif"))


## BIODIVERSITY DATA
## source: GBIF.org (17 September 2018) GBIF Occurrence Download https://doi.org/10.15468/dl.7kootk

# dat <- fread(paste0(data_path, "raw/", "biodiversity/GBIF/tytoalba/0012585-180824113759888.csv"), quote = "")
# backbone <-fread(paste0(data_path, "raw/", "biodiversity/GBIF/backbone/Taxon.tsv"))
# source("R/filter_gbif_data.R")
# dat <- filter_gbif_data(dat, backbone, subset.gbifnubtaxonomy.byclass = "Aves",
#                         output_folder = paste0(data_path, "processed/"), domain.mask = global_mask)
# rm(backbone)
# gc()

dat <- fread(paste0(data_path, "processed/", "filtereted_data_2018-11-19.csv"))
dat <- dat[,.(species, decimallatitude,decimallongitude)]


## DEM VARIABLES
## source: NASA https://webmap.ornl.gov/ogc/dataset.jsp?dg_id=10008_1
## to do: alternate source - IIASA: http://webarchive.iiasa.ac.at/Research/LUC/External-World-soil-database/HTML/global-terrain-slope.html?sb=6
  # srtm <- raster(paste0(data_path, "raw/", "srtm/NASA/sdat_10008_1_20181118_172630534_01d.tif"))  
  # srtm <- crop(srtm, global_mask)
  # names(srtm) <- "srtm"
  # extent(srtm) <- extent(global_mask)
  # elevation <- mask(srtm, global_mask)
  # slope <- terrain(srtm, opt = "slope")
  # roughness <- terrain(srtm, opt = "roughness")
  # aspect <- terrain(srtm, opt = "aspect")
  # terrain_vars <- stack(elevation, slope, roughness, aspect)
  # writeRaster(terrain_vars, paste0(data_path, "processed/", "terrain.tif"), format = "GTiff")
  # rm(terrain_vars, srtm, elevation, slope, roughness, aspect)
  # gc()
  terrain <- raster(paste0(data_path, "processed/", "terrain.tif"))


## SOIL
## source: NASA https://daac.ornl.gov/SOILS/guides/igbp-surfaces.html
  ## see notes re extent() and aggreate() here
  # soil <- stack(list.files(paste0(data_path, "raw/","soil/NASA/IGBPDIS_SURFPRODS_569/data"), pattern = "*.dat", full.names = T)) ## $ specifies end of string
  # soil <- soil[[c(1,3,4,6)]] ## subset data 
  # crs(soil) <- crs(global_mask)
  # soil <- crop(soil, global_mask) ## NOT extent (soil) <- extent(global_mask) because this gives pixels with unequal length and width.
  # soil <- projectRaster(soil, global_mask) # NOT aggregate() because cannot scale up to 0.1 degree with an interger
  # soil <- mask(soil, global_mask)
  # names(soil) <- c("bulk", "pawc", "carb",  "nitro")
  # writeRaster(soil, paste0(data_path, "processed/", "soil.tif"), format = "GTiff")
  # rm(soil)
  # gc()
  soil <- raster(paste0(data_path, "processed/", "soil.tif"))


## TOPOGRAPHIC DIVERSITY... [download from Google Earth Engine]
## source: https://developers.google.com/earth-engine/datasets/catalog/JAXA_ALOS_AW3D30_V1_1
## layer incomplete. Don't use.  
  
## LANDUSE
## source: https://earthexplorer.usgs.gov/
## Recalssified land use classes as per SK's code
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
  # landuse <- raster(paste0(data_path, "raw/","landcover/USGS_GLCC/glccgbe20_tif_geoproj/gbigbpgeo20.tif"))
  # landuse <- crop(landuse, global_mask)
  # landuse <- projectRaster(landuse, global_mask, method = "ngb")
  # landuse <- mask(landuse, global_mask)
  # names(landuse) <- "landuse"
  # landuse[landuse[]%in%c(12,14)] <- 100
  # landuse[landuse[]%in%c(6,8,7,9,10)] <- 200
  # landuse[landuse[]%in%c(1,2,3,4,5)] <- 300
  # landuse[landuse[]%in%c(13)] <- 400
  # landuse[landuse[]%in%c(0, 11,15,16,17)] <- 500
  # landuse <- landuse/100 - 1 ## these steps ensure that substituted alues do not overlap wth existing values
  # table(landuse[])
  # writeRaster(landuse, paste0(data_path, "processed/", "landuse.tif"), format = "GTiff")
  # rm(landuse)
  # gc()
  landuse <- raster(paste0(data_path, "processed/", "landuse.tif"))

  
## CLIMATE - on BOAB
## source: http://worldclim.org/
## code author: Skip - https://github.com/skiptoniam/sense
  
  setwd("./discovery_globalSDMs")
  data_path <- ""
  
  library(sp)
  library(raster)
  library(rgdal)
  devtools::install_github('skiptoniam/sense')
  library(sense)
  
  e <- c(-180,180,-60,90)
  reso <- res(global_mask)
  
  ## Future 2.6, model BCC-CSM1-1, Year 2070
  file_in <- list.files(paste0(data_path, "raw/", "climate/bioclim/bc26bi70/"), pattern = '.tif',full.names = TRUE)
  file_out <-  sub(".tif","crop.tif",file_in)
  mapply(gdal_crop, inpath = file_in, outpath = file_out, MoreArgs = list(extent=e,res=reso))  

    ## Future 8.5, model BCC-CSM1-1, Year 2070
  file_in <- list.files(paste0(data_path, "raw/", "climate/bioclim/bc85bi70/"), pattern = '.tif',full.names = TRUE)
  file_out <-  sub(".tif","crop.tif",file_in)
  mapply(gdal_crop, inpath = file_in, outpath = file_out, MoreArgs = list(extent=e,res=reso))  
  
  ## Current
  download.file("http://biogeo.ucdavis.edu/data/worldclim/v2.0/tif/base/wc2.0_30s_bio.zip", destfile = "./wc2.0_30s_bio.zip")
  if (!dir.exists("wc2.0_30s_bio")) {
    dir.create("wc2.0_30s_bio")
  }
  unzip(".wc2.0_30s_bio.zip", exdir = "./wc2.0_30s_bio")
  
  file_in <- list.files(paste0(data_path, "./wc2.0_30s_bio/"), pattern = '.tif',full.names = TRUE)
  file_out <-  sub(".tif","crop.tif",file_in)
  mapply(gdal_crop, inpath = file_in, outpath = file_out, MoreArgs = list(extent=e,res=reso))  


## BIAS LAYER
## Distance to  roads, built up/population density, coastlines, pas, water
## options - Roozbeh's function using sf package (too long); Casey SQL function in postGIS
  
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

