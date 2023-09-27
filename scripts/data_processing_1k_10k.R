## Envirobnmental data processing for global analyses @1k @10k
## NOTE: Data sources recorded in data_downloads.R
## Note: -ot in gdal funcrions specifies Byte vs Float32 depending on input data



## ------------------------------------------------------------------------- ##
## Set working environment ####
## ------------------------------------------------------------------------- ##
rm(list = ls())
gc()
# system("ps")
# system("pkill -f R")


## Load libraries
# install.packages("remotes")
# remotes::install_github("skiptoniam/geostuff")
x <- c('data.table', 'sp', 'raster', 'terra',
       'rgdal', 'geostuff', 'tools', 'bitops', 
       'ncdf4', 'processNC',
       'RCurl', 'gdalUtils', 'usethis',
       'parallel', 'doMC',
       'reticulate')
lapply(x, require, character.only = TRUE)
rm(x)

mc.cores = future::availableCores()-2


## >> Specify global variables ####
proj.res.km <- 10
proj_res <- proj.res.km*1000
# proj_res <- c(1000,1000) ## 1km res

## >> Projections ####
## >> For mapping
# map_crs <- "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"

## >> For analyses: Equal Earth
## Ref for EE: https://proj.org/operations/projections/eqearth.html#id2
## jgarber note the '+proj=eqearth' is a new format and wont work wth gdalUtils
## instead we use '+proj=eqearth +ellips=WGS84 +wktext' which give equal earth and works here
## see: https://en.wikipedia.org/wiki/Equal_Earth_projection
equalearth_crs = '+proj=eqearth +ellips=WGS84 +wktext'

## >> For analyses: WGS
wgs_crs = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0" 


## >> Folder paths ####
skip_dir <- "/home/payalb/gsdms_r_vol/tempdata/mediaflux/proj-2200_nature_futures21-1128.4.411/Alex/aus-ppms_alex/outputs/"

gsdms_dir <- "/tempdata/research-cifs/6300-payalb/uom_data/gsdms_data"
output_dir <- file.path(gsdms_dir, "outputs", sprintf("layers_%sk", proj.res.km))

if(!dir.exists(output_dir)) {
  dir.create(output_dir)
}



## >> Functions ####
# source("/home/payalb/gsdms_r_vol/tempdata/workdir/gsdms/scripts/0_functions.R")
# source("/home/payalb/gsdms_r_vol/tempdata/workdir/landuse_projects/land-use/gdal_calc.R") # by jgarber
# ## To pull directly from jgarber's repo:
# library(RCurl)
# url <- "https://gitlab.unimelb.edu.au/garberj/gdalutilsaddons/-/blob/master/gdal_calc.R"
# url.exists(url)
# source(url)
# devtools::source_url(url)
# url <- "https://gitlab.unimelb.edu.au/garberj/gdalutilsaddons.git"
# devtools::install_git(url = url)

## Change NoData values to -9999
## ref: https://gis.stackexchange.com/questions/298230/change-no-data-value-geotif-file-with-qgis-or-gdal
change_values = "/home/payalb/gsdms_r_vol/tempdata/workdir/gsdms/scripts/change_values.sh"



## ------------------------------------------------------------------------- ##
## Create mask ####
## ------------------------------------------------------------------------- ##
## >> crs for Worldclim layers: https://worldclim.org/data/v1.4/formats.html

## >> Processing in repo: https://github.com/daarias331/dp_repo/tree/main




## ------------------------------------------------------------------------- ##
## Land use layers ####
## ------------------------------------------------------------------------- ## 
## >> Hurtt data: https://luh.umd.edu/
## >> TIme steps available: 2015 - 2100 (86 time steps)


## >> Processing in repo: https://github.com/daarias331/dp_repo/tree/main
## >> Layers processed for ime steps : 2015 2020 2025 2030 2035 2040 2045 2050 2055 2060 2065 2070



## ------------------------------------------------------------------------- ##
## Prepare future Worldclim climate layers ####
## ------------------------------------------------------------------------- ##
## >> Worlclim data: https://www.worldclim.org/data/cmip6/cmip6_clim30s.html
## >> Time steps available: 2021-2040, 2041-2060, 2061-2080, 2081-2100

## >> Processing in repo: https://github.com/daarias331/dp_repo/tree/main




## ------------------------------------------------------------------------- ##
## Calculate distance layers ####
## ------------------------------------------------------------------------- ##

## ------------------------------------------------------------------------------------------ ##
## NOTE: Layers provided to Alejandro for pre-processing
# Builtup (90.5 MB)
# “/research-cifs/6300-payalb/uom_data/gsdms_data/builtup/GHS_BUILT_LDS2014_GLOBE_R2018A_54009_1K_V2_0.tif”
# 
# Roads (7.57 GB)
# “/research-cifs/6300-payalb/uom_data/gsdms_data/groads/groads.tif”
# 
# Pop density (369.7 MB)
# “/research-cifs/6300-payalb/uom_data/gsdms_data/population/SEDAC/gpw-v4-population-density-rev11_2020_30_sec_tif/gpw_v4_population_density_rev11_2020_30_sec.tif”
# 
# Wetlands (933.3 MB)
# “/research-cifs/6300-payalb/uom_data/gsdms_data/lakesrivers/GLWD/wetlands.tif”
# 
# Lakes & rivers (933.3 MB)
# “/research-cifs/6300-payalb/uom_data/gsdms_data/lakesrivers/GLWD/lakesrivers.tif”
# 
# Protected areas (9.9 GB)
# “/research-cifs/6300-payalb/uom_data/gsdms_data/protectedareas/wdpa.tif”
# 
# 
# Idea is to tun your pre-processing pipeline and provide me with tifs for Australia as per the Australia Albers mask (“/research-cifs/6300-payalb/uom_data/aus-ppms_data/outputs/mask/ausmask_nesp_250mAlbersEA.tif”) and I will thereafter run the gdal_proximity.py function on them.
## ------------------------------------------------------------------------------------------ ##


dst_folder <- file.path(gsdms_dir, "outputs", "distance_layers_wgs")
if(!dir.exists(dst_folder)) {
  dir.create(dst_folder)
}

# 1k: /Volumes/proj-2200_nature_futures21-1128.4.411/Alex/gsdms_alex/outputs/processed/1km/dist_layers
# 10k: /Volumes/proj-2200_nature_futures21-1128.4.411/Alex/gsdms_alex/outputs/processed/10km/dist_layers



### Builtup ####
#### Explore raw data ####
## >> Values are expressed as decimals (Float) from 0 to 100 (see p14 in  GHSL_Data_Package_2019_light.pdf)
infile <- file.path(gsdms_dir, "builtup", "GHS_BUILT_LDS2014_GLOBE_R2018A_54009_1K_V2_0.tif")
gdalinfo(infile)
gdalsrsinfo(infile)


#### Pre-processing ####
## >> Processing in repo: https://github.com/daarias331/dp_repo/tree/main


#### Calculate distance layer ####
## >> Distance is calculated from all non 0 values
infile <- ... # specify pre-processing output file
outfile <- file.path(dst_folder, "dst_builtup.tif")

job::job({
  system(sprintf("gdal_proximity.py %s %s -of GTiff -distunits GEO -nodata %d -co compress=LZW -co BIGTIFF=YES", infile, outfile, -9999))
})

gdalinfo(outfile)
# system(sprintf("gdalinfo %s", outfile))


#### Mask ? ####



### Roads ####
#### Explore raw data ####
## >> https://gis.stackexchange.com/questions/151613/reading-feature-class-in-file-geodatabase-using-r
roadsgdb <- file.path(gsdms_dir, "groads/groads-v1-global-gdb/gROADS_v1.gdb")

## List all feature classes in a file geodatabase 
subset(ogrDrivers(), grepl("GDB", name))
ogrListLayers(roadsgdb)


#### Subset the road attributes ####
## >> see gROADSv1_documentation.pdf
## >>   1. FClass = Functional class (1=Highway, 2=Primary, 3=Secondary, 4=Tertiary, 5=Local/ Urban, 6=Trail, 7=Private, 0=Unspecified)
## >>   2. EXS = Existence Category (1=Definite, 2=Doubtful, 0=Unspecified) ## do not use this caregory, north canada and UK excluded

roads <- readOGR(dsn=roadsgdb,layer="Global_Roads")
summary(roads)
unique(roads$FCLASS)
unique(roads$EXS)

# roadsub <- roads[roads$FCLASS == (0:3),] ## groads03 files
roadsub <- roads[roads$FCLASS == (0:2),]
unique(roadsub$FCLASS)
unique(roadsub$EXS)

## Drop columns
roadsub <- roadsub[, 4]


#### Write shapefile ####
writeOGR(roadsub, dsn = file.path(gsdms_dir, "groads"), layer = "groads", driver="ESRI Shapefile")


#### Rasterise @ .0025 ~ 250m ####
system(paste0("gdal_rasterize -at -burn 1 -ot Byte -tr .0025 .0025 ",
              file.path(gsdms_dir, "groads", "groads.shp "),
              file.path(gsdms_dir, "groads", "groads.tif")))
gdalinfo(file.path(gsdms_dir, "groads", "groads.tif"))
rm(roads, roadsub, roadsgdb)


#### Pre-processing ####
## >> Processing in repo: https://github.com/daarias331/dp_repo/tree/main


#### Calculate distance layer ####
## >> Distance is calculated from all values of 1
infile <- ... # specify pre-processing output file
gdalinfo(infile)
outfile <- file.path(dst_folder, "dst_roads.tif")

job::job({
  system(sprintf("gdal_proximity.py %s %s -values %d -of GTiff -distunits GEO -nodata %d -co compress=LZW -co BIGTIFF=YES", infile, outfile, 1, -9999))
})

gdalinfo(outfile)


#### Mask ? ####


### Population density ####
## >> Note: Population density rasters were created by dividing the population count raster for a given target year by the land area raster (Float): see https://sedac.ciesin.columbia.edu/data/set/gpw-v4-population-density-rev11

#### Explore raw data ####
infile <- file.path(gsdms_dir, "population/SEDAC/gpw-v4-population-density-rev11_2020_30_sec_tif", "gpw_v4_population_density_rev11_2020_30_sec.tif")
gdalinfo(infile)


#### Pre-processing ####
## >> Processing in repo: https://github.com/daarias331/dp_repo/tree/main


#### Calculate distance layer ####
## >> Distance is calculated from all values for pop density of 60
infile <- ... # specify pre-processing output file
outfile <- file.path(dst_folder, "dst_pop.tif")

job::job({
  system(sprintf("gdal_proximity.py %s %s -values %d -of GTiff -distunits GEO -nodata %s -co compress=LZW -co BIGTIFF=YES", infile, outfile, 60, -9999))
})

gdalinfo(outfile)


#### Mask ? ####



### Water bodies ####
## >> Source: https://www.worldwildlife.org/publications/global-lakes-and-wetlands-database-lakes-and-wetlands-grid-level-3
## >> Note: Projection of original dataset is lat/long (see p7 of GLWD_Data_Documentation.pdf)
infile <- file.path(gsdms_dir, "lakesrivers/GLWD/glwd_3")
gdalinfo(infile)


#### Lakes & rivers ####
##### Extract rivers and lakes from data ####
outfile <- file.path(gsdms_dir, "lakesrivers/GLWD", "lakesrivers.tif")

job::job({
  system(paste0("gdal_calc.py -A ", infile,
                " --calc='((A==1) + (A==2) + (A==3) + (A==4) + (A==5))*1 + ( + (A==6) + (A==7) + (A==8) + (A==9) + (A==10) + (A==11) + (A==12))*0' --NoDataValue=-9999",
                " --outfile=", outfile))
})

gdalUtils::gdalinfo(outfile)


##### Pre-processing ####
## >> Processing in repo: https://github.com/daarias331/dp_repo/tree/main


##### Calculate distance layer ####
## >> Distance is calculated from all values of 1
infile <- ... # specify pre-processing output file
outfile <- file.path(dst_folder, "dst_lakesrivers.tif")

job::job({
  system(sprintf("gdal_proximity.py %s %s -values %d -of GTiff -distunits GEO -nodata %d -co compress=LZW -co BIGTIFF=YES", infile, outfile, 1, -9999))
})

gdalinfo(outfile)


##### Mask ? ####



#### Wetlands ####
##### Extract wetlands from data ####
infile <- file.path(gsdms_dir, "lakesrivers/GLWD/glwd_3")
outfile <- file.path(gsdms_dir, "lakesrivers/GLWD", "wetlands.tif")

job::job({
  system(paste0("gdal_calc.py -A ", infile,
                " --calc='((A==1) + (A==2) + (A==3) + (A==4) + (A==5))*0 + ( + (A==6) + (A==7) + (A==8) + (A==9) + (A==10) + (A==11) + (A==12))*1' --NoDataValue=-9999",
                " --outfile=", outfile))
})

gdalUtils::gdalinfo(outfile)

##### Pre-processing ####
## >> Processing in repo: https://github.com/daarias331/dp_repo/tree/main


##### Calculate distance layer ####
## >> Distance is calculated from all values of 1
infile <- ... # specify pre-processing output file
outfile <- file.path(dst_folder, "dst_wetlands.tif")

job::job({
  system(sprintf("gdal_proximity.py %s %s -values %d -of GTiff -distunits GEO -nodata %d -co compress=LZW -co BIGTIFF=YES", infile, outfile, 1, -9999))
})

gdalinfo(outfile)


##### Mask ? ####



### Protected areas ####
## >> Input data format: Byte (0/1 data)
## >> Output data format: Floast32
## >> Distance is calculated from all non 0 values (only 1 and 0 values in this layer)

#### Rasterise ####
## >> Tif file created in data_downloads.R

#### Pre-processing ####
## >> Processing in repo: https://github.com/daarias331/dp_repo/tree/main


#### Calculate distance layer ####
## >> Distance is calculated from all values of 1
infile <- ... # specify pre-processing output file
outfile <- file.path(dst_folder, "dst_wdpa.tif")

job::job({
  system(sprintf("gdal_proximity.py %s %s -values %d -of GTiff -distunits GEO -nodata %s -co compress=LZW -co BIGTIFF=YES", infile, outfile, 1, -9999))
})

gdalinfo(outfile)


#### Mask ? ####


### Check all dst layers
list.files(file.path(gsdms_dir, "distance_layers_wgs"),
           full.names = TRUE)
lapply(list.files(file.path(gsdms_dir, "distance_layers_wgs"),
                  full.names = TRUE), gdalsrsinfo)




## ------------------------------------------------------------------------- ##
## Prepare soil layers ####
## ------------------------------------------------------------------------- ##

## >> Processing in repo: https://github.com/daarias331/dp_repo/tree/main




## ------------------------------------------------------------------------- ##
## Prepare SRTM ####
## ------------------------------------------------------------------------- ##

## >> Processing in repo: https://github.com/daarias331/dp_repo/tree/main




## ------------------------------------------------------------------------- ##
## Prepare Hansen tree cover ####
## ------------------------------------------------------------------------- ##

## >> Processing in repo: https://github.com/daarias331/dp_repo/tree/main



## ------------------------------------------------------------------------- ##
## Check layers ####
## ------------------------------------------------------------------------- ##
# .. use code from aus-ppms