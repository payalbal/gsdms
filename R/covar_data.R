## Download and prepare covariate data for gsdms


## Load libraries
## For gdal to work on Mac in R, see: https://github.com/smwindecker/gdaltools
# devtools::install_github('skiptoniam/sense')
pacman::p_load(sp, raster, rgdal, sense, tools, bitops, RCurl)


## Set environment
setwd()
# data_covs <- "./data" # boab
data_covs <- "/Volumes/discovery_data/gsdms_data" # local


## WorldClim data - download tiles and stitch to make global layer
## Source: https://www.worldclim.org/version1
## author: Chris Ware

## set path to download each zip archive to (will be deleted once used)
## must have .zip ext
zipdst = file.path(data_covs, 'bio_30s.zip')

## set dir to extract raster files to
rasterdst = file.path(data_covs, 'bio_30s/')
if(!dir.exists(rasterdst)) {
  dir.create(rasterdst)
} # end if !dir.exists

## set path to write the global raster for each bioclim variable to
bioclim_dst = paste0(rasterdst, 'bio_current_') # or whatever

## Looking at the getData function, it finds which tile a given lat lon  falls witin. 
## The tiles are then indexed by row/col which maps to a url. 
## Can cut to the chase and just create urls with all combos of row/cols (i.e. tile ids)
urls_bioclim = NULL
for (r in 1:5){
  for (c in 1:12){
    
    rc = paste0(r-1, c-1)
    fn = paste0('bio_', rc, '.zip')
    thisurl = paste0('https://biogeo.ucdavis.edu/data/climate/worldclim/1_4/tiles/cur/', fn)
    urls_bioclim = c(urls_bioclim, thisurl)
    
  }
}

## check all the urls exist (which hopefully means all the tiles can be downloaded)
require(bitops)
require(RCurl)
test = lapply(urls_bioclim, url.exists)
all(unlist(test))

## Then, it's possible to loop over the urls, downloading, and unzipping them. 
download_from_url <- function (urls, zipdst, rasterdst) {
  
  for (url in urls){
    response = tryCatch({
      download.file(url, destfile = zipdst)
    }, 
    error = function(e){e}) # shouldn't get called
    
    print(response) # should be 0
    unzip(zipdst, exdir = rasterdst)
    file.remove(zipdst)
  }
} 

download_from_url(urls = urls_bioclim, zipdst = zipdst, rasterdst = rasterdst)

## rasterdst should now be full of bioclim 1-19 tiles (60 tiles for each var)
## loop over each bioclim variable, collect all tiles, and mosaic. 
for (i in 1:19){
  f = list.files(rasterdst, pattern = paste0('bio', i, '_'), full.names = TRUE)
  f = f[grep('.bil$', f)]
  
  f = lapply(f, raster)
  
  ## specify a function for any overlapping cells (of which there won't be any, 
  ## but the function requires a function...)
  f$fun = mean
  mos = do.call(mosaic, f)
  
  this_dst = paste0(bioclim_dst, i, '.tif') # or whatever
  writeRaster(mos, this_dst)
}


## Future Bioclim for Year: 2070, GCM: BCC-CSM1-1
url_bc85bi70 = "http://biogeo.ucdavis.edu/data/climate/cmip5/30s/bc85bi70.zip"
url_bc26bi70 = "http://biogeo.ucdavis.edu/data/climate/cmip5/30s/bc26bi70.zip"

download_from_url(urls = c(url_bc26bi70, url_bc85bi70), zipdst, rasterdst)


bio_current <- list.files(paste0(data_covs, "/bio_30s"), pattern = "bio_current*", full.names = TRUE)
bio_rcp26 <- list.files(file.path(data_covs, "bio_30s"), pattern = "*bc26*", full.names = TRUE)
bio_rcp85 <- list.files(file.path(data_covs, "bio_30s"), pattern = "*bc85*", full.names = TRUE)
bioclim <- c(bio_current, bio_rcp26, bio_rcp85)


## SRTM
## Source: https://topotools.cr.usgs.gov/gmted_viewer/gmted2010_global_grids.php
system("wget https://edcintl.cr.usgs.gov/downloads/sciweb1/shared/topo/downloads/GMTED/Grid_ZipFiles/mn30_grd.zip -O data/srtm.zip")
unzip(file.path(data_covs, "srtm.zip"), list = TRUE)
unzip(file.path(data_covs, "srtm.zip"), exdir = file.path(data_covs, "srtm"))
file.remove(file.path(data_covs, "srtm.zip"))
  ## all files within the mn30_grd folder seem to be the same when in as raster
srtm <- file.path(data_covs, "srtm/mn30_grd/w001000.adf")

    ## Alternate SRTM source: http://srtm.csi.cgiar.org/srtmdata/ (extent= -180,180,-60,60)
    # tile_names <- as.vector(unlist(as.data.frame(outer(paste0(1:72, "_"), 1:24, FUN = "paste0"))))
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
    
    ## Alterate SRTM source: https://webmap.ornl.gov/ogc/dataset.jsp?dg_id=10008_1
    ## download by lower resolution and resample or download by smaller extents and stitch tiles together


## Soil
## Source: https://daac.ornl.gov/SOILS/guides/igbp-surfaces.html
## Click 'Get Data' on above link
## Select data, e.g. bulkdens,pawc,soilcarb,totaln and 'Add checked items'
  ## bulkdens - bulk density (weight of soil in a given volume)
  ## pawc - profile available water capacity (amount of water that can be stored in a soil profile)
  ## soilcarb - soil carbon density 
  ## totaln - total nitrogen density
## Click 'order checked items'
## Get the link to download data from email, which will look like
## this:  https://daac.ornl.gov/orders/<some_number>/
## Run the system(wget...) command with the retrieved link
system("wget -r -np -nH --reject 'index.html*' -e robots=off https://daac.ornl.gov/orders/71a04c738fde75ee64f57118ed466530/")
soil <- list.files(file.path(data_covs,"orders"), pattern = "*.dat", full.names = T, recursive = T)


## Landuse: https://earthexplorer.usgs.gov/
## Get the ftp site for the data: https://edcftp.cr.usgs.gov/project/glcc/globdoc2_0.html
## Go >> https://edcftp.cr.usgs.gov/project/glcc/globdoc2_0.html#anony
## Readme: https://www.usgs.gov/media/images/global-land-cover-characteristics-data-base-version-20
system("wget -r -np -nH --reject 'index.html*' -e robots=off https://edcftp.cr.usgs.gov/project/glcc/globe") #corrupted...


system("curl https://edcftp.cr.usgs.gov/project/glcc/globe/gigbp1_2.img.gz -o data/gigbp1_2.img.gz") ##corrupted...
R.utils::gunzip(file.path(data_covs, "gigbp1_2.img.gz"), destname = file.path(data_covs, "landcover.img"), remove = FALSE)
landuse <- file.path(data_covs, "gigbp1_2.img")
landuse_out <- file.path(".data/gigbp1_2.img")
gdal_crop(inpath = landuse, outpath = landuse_out, extent=e, res=reso)


unzip(file.path(data_covs, "gigbp1_2.img.gz"), exdir = file.path(data_covs, "landcover"))
landuse <- list.files(file.path(data_covs,"landcover"), pattern = c("*igbp", "*.tif"), full.names=T)
landuse <- landuse[3] 

  ## Hurtt dataset...

r <- brick("/Volumes/discovery_data/gsdms_data/gigbp1_2.img")

## DATA PROCESSING
## Create Mask from WorldClim layer
global_mask <- raster(bio_current[1])
global_mask[which(!is.na(global_mask[]))] <- 1
writeRaster(global_mask, filename = paste0(data_gsdms, "/globalmask10m.tif"))


## sense::gdal_crop
file_in <- c(bioclim, srtm, soil, landuse)
file_out <- paste0(data_covs, "/processed/", paste0(tools::file_path_sans_ext(basename(file_in)), "_treated.", tools::file_ext(file_in)))
file_out[58] <- sub(file_path_sans_ext(basename(file_out[58])), "srtm_treated", file_out[58]) # rename srtm
file_out[63] <- sub(file_path_sans_ext(basename(file_out[63])), "landuse_treated", file_out[63]) #renname landuse

e <- c(-180,180,-60,90)
reso <- res(global_mask)
mapply(gdal_crop, inpath = file_in, outpath = file_out, MoreArgs = list(extent=e,res=reso)) 


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
# check
unique(values(landuse))

## STACK, MASK & SAVE COVARIATES AS RDS
covariates_all <- stack(setdiff(file_out, file_out[grep("srtm|landuse", file_out)]), 
                        elevation, slope, roughness, aspect, landuse) 
## load and stack covariates files, excluding files for srtm and landuse
##  which were updated and therefore these rasters are loaded from memory.



rm(list=setdiff(ls(), c("data_gsdms", "covariates_all", "global_mask")))
crs(covariates_all) <- crs(global_mask)
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

