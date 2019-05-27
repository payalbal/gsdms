
## ****** To update see discovery_paper_1 and SK's birds_ccimpacts scripts on GitHib ****** ##


## BARN OWL ANALYSIS @ 100km2 resolution
## ---------------------------------------
## @description: Data prepapration for global SDMs for tyto alba at 100 km2 resolution
## base year for GTAP trajectories, population and landuse = ??


library(data.table)
library(sp)
library(raster)
library(rgdal)
library(rgeos)
library(sf)
library(dplyr)
library(progress)


## Specify data path
data_path <- "/Volumes/payal_umelb/data/"


## MASK
## source: DEM.tif from IUCN basedata: https://www.iucnredlist.org/resources/spatialtoolsanddata
## to do: add ice mask: see Hoskins et al 

# global_mask <- raster(paste0(data_path, "raw/", "biodiversity/IUCN/basedata/DEM/DEM.tif"))
# global_mask[which(!is.na(global_mask[]))] <- 0
# writeRaster(global_mask, paste0(data_path, "processed/", "global_mask"), format = "GTiff", options="COMPRESS=LZW", overwrite = TRUE)
# global_mask <- raster(paste0(data_path, "processed/", "global_mask.tif"))
# global_mask100 <- aggregate(global_mask, fact=12, fun=mean, na.rm=TRUE, filename=paste0(data_path, "processed/", "global_mask100.tif"))
#  ## note: factor 12 in aggregate() to increase resolution from 0.008... to 0.1 degrees i.e. approx 10 km s.t. area of a pixel is 100 km2
# writeRaster(global_mask100, paste0(data_path, "processed/", "global_mask100"), format = "GTiff", options="COMPRESS=LZW", overwrite = TRUE)

global_mask <- raster(paste0(data_path, "processed/", "global_mask100.tif"))

## BIODIVERSITY DATA
## source: GBIF.org (17 September 2018) GBIF Occurrence Download https://doi.org/10.15468/dl.7kootk

# dat <- fread(paste0(data_path, "raw/", "biodiversity/GBIF/tytoalba/0012585-180824113759888.csv"), quote = "")
# backbone <-fread(paste0(data_path, "raw/", "biodiversity/GBIF/backbone/Taxon.tsv"))
# source("R/filter_gbif_raw.R")
# dat <- filter_gbif_raw(dat, backbone, subset.gbifnubtaxonomy.byclass = "Aves",
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
## to do - alternate source - IIASA: http://webarchive.iiasa.ac.at/Research/LUC/External-World-soil-database/HTML/index.html?sb=1
##          local downloaded copy at /payal_unimelb/data/raw/soil/IIASA/HWSD_RASTER
## notes: This dataset needs to be linked to the accompanying .mdb database in QGIS - look up r packages to do this (Casey for help)
## Or install mdbtools via terminal using 'brew install mdbtools'...not sure what this does 
## see https://brewinstall.org/Install-mdbtools-on-Mac-with-Brew/

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



## POPULATION DENSITY (2005)
## source: http://sedac.ciesin.columbia.edu/data/set/gpw-v4-population-density-rev10/data-download
## to do: alternate source - use urban settlement map?

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

pop <- raster(paste0(data_path, "processed/","population.tif"))



## LANDUSE
## source: https://earthexplorer.usgs.gov/
## issue: this is land cover and not land use. See Chris Ware's email dated 20 Nov 2018 
##  + Marco notes going from land cover to land use is problematic
## Use Hoskins (5 classes) or Hurtt's LUH2 (12 classes) instead to align work with CSIRO's 
## CSIRO link: https://data.csiro.au/dap/landingpage?pid=csiro:15276&v=3&d=true 
## description of primary vs secondary forests: http://luh.umd.edu/faq.shtml 
## as per Hoskins, can add a 'permanent ice' mask and do away with the barren class 
##  see Olson in Hoskins paper (local copy of Olson data on payal_unimelb) for ice mask
## [alternate datasets to explore]: http://webarchive.iiasa.ac.at/Research/LUC/External-World-soil-database/HTML/LandUseSharesDesc.html?sb=8

# landuse <- raster(paste0(data_path, "raw/","landcover/USGS_GLCC/glccgbe20_tif_geoproj/gbigbpgeo20.tif"))
# landuse <- crop(landuse, global_mask)
# landuse <- projectRaster(landuse, global_mask, method = "ngb")
# landuse <- mask(landuse, global_mask)
# names(landuse) <- "landuse"

## Reclassify land use classes as per SK's code
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

## Landuse plot
# # devtools::install_github("ropenscilabs/ochRe")
# library(ochRe)
# class.palette <- colorRampPalette(c(ochre_palettes[["lorikeet"]][3], ochre_palettes[["parliament"]][c(2,3)], "red", ochre_palettes[["parliament"]][6]))
# par(xpd=FALSE)
# plot(landuse, col=class.palette(5), legend = NULL, axes=FALSE, box=FALSE)
# par(xpd=TRUE)
# class.palette <- colorRampPalette(c(ochre_palettes[["lorikeet"]][3], ochre_palettes[["parliament"]][c(2,3)], "red4", ochre_palettes[["parliament"]][6]))
# legend(-190,-20, legend = c("crop", "grass/shrub", "forest", "urban", "bare"), fill = class.palette(5), bty = "n", cex = 0.7)
# 
# plot(NULL)
# legend("center", legend = c("crop", "grass/shrub", "forest", "urban", "bare"), fill = class.palette(5), bty = "n", cex = 3)
# 
# # mapview::mapview(landuse)
# tyto <- SpatialPointsDataFrame(dat[,c(3,2)], data = dat, proj4string = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
# tyto <- spTransform(tyto, proj4string(global_mask))
# points(tyto , pch=".", cex=0.1, col="maroon")



## DISTANCE TO FEATURE LAYERS
## options - Roozbeh's function using sf package (too long); Casey SQL function in postGIS

## built up areas
## coastlines


## Roads
## FOR NOW
## Processing in ArcGIS and QGIS
## 1. Load .gdb in ArcGIS and filter by ("EXS" = '0'  OR  '1')  AND ("FCLASS" = 0 OR 1 OR 2 OR 3 OR 4 OR 5)
##    Save subset as a vector layer
## 2. In QGIS create a distance-to-feature raster: Toolbox > GDAL > Raster analysis > Proximity

## FOR LATER 
## to read .gdb in r see: https://gis.stackexchange.com/questions/184013/read-a-table-from-an-esri-file-geodatabase-gdb-using-r
## layer name for gROADS_v1.gdb is Global_Roads (obtained from reading the file in QGIS)
## subset the road attributes (see gROADSv1-ReadMe.txt) 
## 1. FClass = Functional class (1=Highway, 2=Primary, 3=Secondary, 4=Tertiary, 5=Local/ Urban, 6=Trail, 7=Private, 0=Unspecified)
## 2. EXS = Existence Category (1=Definite, 2=Doubtful, 0=Unspecified)

roads <- sf::st_read(dsn = paste0(data_path, "raw/", "roads/groads-v1-global-gdb/gROADS_v1.gdb"), layer = "Global_Roads")
sub_road <- roads %>% 
  filter(FCLASS == (0:5) & (EXS == c(0,1)))
plot(global_mask)
plot(st_geometry(sub_road))
source("R/rasterDistance.R")
...
distroads <- rasterDistance(sub_road, global_mask)
# writeRaster(distroads, paste0(data_path, "processed/", "distroads.tif"), format = "GTiff")
# rm(distroads, roads, sub_road)
# gc()



## Protected areas
pa <- shapefile(paste0(data_path, "raw/", "boundaries/WDPA_Nov2018-shapefile/WDPA_Nov2018-shapefile-polygons.shp"))

... 


## Water bodies (including coastline, lakes, rivers)
## source: https://www.ngdc.noaa.gov/mgg/shorelines/
ras1 <- shapefile(paste0(data_path, "raw/", "protectedareas/GSHHG/gshhg-shp-2.3.7/GSHHS_shp/f/GSHHS_f_L5.shp"))
require(rgdal)
require(sf)
shape <- read_sf(dsn = paste0(data_path, "raw/", "protectedareas/GSHHG/gshhg-shp-2.3.7/GSHHS_shp/f/GSHHS_f_L5.shp"), layer = "SHAPEFILE")








## ------------------------------------------------------------------------------------
##
## NOTES & ISSUES
## install.packages("data.table",type="source") # to fix issue of data.table not updating in the enviroment. See https://stackoverflow.com/questions/26704837/r-does-not-update-data-table
## reading in .asc or .tif files (e.g. in loading srtm data) with raster makes them similar sized raster objectes even though the size of .asc >> .tif files.
## Projection used and the implications of equal-area assumption for global analysis. possible solution - land_area spatial layer may be used (??) available at https://daac.ornl.gov/ISLSCP_II/guides/global_population_xdeg.html OR http://sedac.ciesin.columbia.edu/data/set/gpw-v4-land-water-area-rev10


## NOTES FOR 1km2 ANALYSIS
## Increasing speed of raster processing
## link: https://www.gis-blog.com/increasing-the-speed-of-raster-processing-with-r-part-13/
## link: https://stackoverflow.com/questions/38368939/rasteroptions-difference-between-chunksize-and-maxmemory
## ALSO: think of converting float to interger or redusing the number of decimals.
##
## Spatial data processing notes for function
## standard practice: 
## download > import raw > check projection from source
## 1. first specify projection of input raster using crs(raster)
##    if projection is not WGS84 use projectRaster to reproject to desired projection
## 2: crop > extent > mask
## 3. to increase resolution of conitnuous data, use aggregate() or resample() or projectRaster()
##      aggregate() is faster BUT can only scale up by interger values
##      resample(method="bilenear") is same as aggregate()
##      see: https://gis.stackexchange.com/questions/255150/using-resample-vs-aggregate-extend-in-r-to-have-rasters-of-matching-resolutio
##      alternately, projectRaster can be used if mask with desired resolution is provided; 
##        uses resample(method="bilenear") by default
##        output size might be bigger...
##    to increase resolution of categorical data, use resample(method="ngb")
##
##
## calculate distance ro roads
##  https://gis.stackexchange.com/questions/233443/finding-distance-between-raster-pixels-and-line-features-in-r
##  
# roads <- crop(temp, global_mask)
# dd = gDistance(as(road_sub, "Spatial"), as(global_mask,"SpatialPoints"), byid=TRUE)
# dist.roads <- st_distance(roads, as(global_mask,"SpatialPoints"), by_element = TRUE)


