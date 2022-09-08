## Covariate data processing for global analyses
## NOTE: Data sources recorded in data_downloads.R
## NOTE: -ot in gdal functions specifies Byte vs Float32 depending on input data

## NOTES FOR ALEJANDRO:
## Search $$ for specific notes where something needs to be fixed or noted
## I have included outputs of sprtinf inn the script as i use this fuction to create my function call i.e. text string (ref for sprintf: https://www.rdocumentation.org/packages/base/versions/3.6.2/topics/sprintf)
## I have commented out the 'Calculate distance layers' section. Skip this. I have left it in, in case we need to recreate the distance layers
## Run script till L556
## This script lists the steps needed for processing of global layers; use it to understand the steps but you needn't stick to the exact text strings and the solutions found therein. There are mutoiple ways to skin a cat they say, so feel free to use the appropriate approach as long as we get tot he desired result :)
## Task 3 will be to then go onto the landuse layer processing, i.e. L566 onwards. We'll discuss this when we come to it. 


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

## >> step one_create mask from WorldClim layer ####
## crs for Worldclim layers: https://worldclim.org/data/v1.4/formats.html

  ## $$$ Note for Alejandro: terra is an R function form the terra package, you can convert this to a python function. Basically I am taking one bioclim layer and creating a new layer based on that, such that the new layer has value = 1 wherever bioclim layer has a data value and the new layer has value = NA (as -9999 to be consistent with out previous scripts? or NA if you feel this is better) wherever bioclim layer has a no_data value. Makes sense?

bio_current <- list.files(paste0(gsdms_dir, "/bio_30s"), pattern = "bio_current*", full.names = TRUE)
global_mask <- raster(bio_current[1])
global_mask[which(!is.na(global_mask[]))] <- 1
crs(global_mask) <- wgs_crs

terra::writeRaster(global_mask, 
                   filename = file.path(output_dir, "mask", "globalmask_wgs_30s.tif"),
                   NAflag = -9999)

unique(getValues(global_mask)) ## NA 1
gdalinfo(file.path(output_dir, "mask", "globalmask_wgs_30s.tif")) ## NoData Value=-9999
rm(bio_current, global_mask)


## >> step two_clip mask to desired extent ####
infile <- file.path(output_dir, "mask", "globalmask_wgs_30s.tif")
outfile <- file.path(output_dir, "mask", "globalmask_wgs_30s_clip.tif")
new_extent <- "-180 -60 180 90"

system(sprintf("gdalwarp -overwrite -ot Byte -te %s %s %s", 
               new_extent, infile, outfile))
  ## Note for Alejandro: output of sprintf below
  ## "gdalwarp -overwrite -ot Byte -te -180 -60 180 90 /tempdata/research-cifs/6300-payalb/uom_data/gsdms_data/outputs/layers_10k/mask/globalmask_wgs_30s.tif /tempdata/research-cifs/6300-payalb/uom_data/gsdms_data/outputs/layers_10k/mask/globalmask_wgs_30s_clip.tif"

raster(infile)
raster(outfile)

gdalinfo(outfile) ## NoData Value=-9999
unique(values(raster(outfile))) ## 0 1


## >> step three_reproject mask to Equal Earth and desired resolution ####
infile <- outfile
outfile <- gsub("wgs_30s_clip", paste0("ee_", proj.res.km, "k"), infile)
new_res <- proj_res
new_crs = equalearth_crs

system(sprintf("gdalwarp -overwrite -ot Byte -tr %d %d -t_srs '%s' %s %s",
               new_res, new_res, new_crs, infile, outfile))
  ## Note for Alejandro: output of sprintf below
  ## "gdalwarp -overwrite -ot Byte -tr 10000 10000 -t_srs '+proj=eqearth +ellips=WGS84 +wktext' /tempdata/research-cifs/6300-payalb/uom_data/gsdms_data/outputs/layers_10k/mask/globalmask_wgs_30s_clip.tif /tempdata/research-cifs/6300-payalb/uom_data/gsdms_data/outputs/layers_10k/mask/globalmask_ee_10k.tif"

raster(infile)
raster(outfile)

gdalinfo(outfile) ## NoData Value=-9999
unique(values(raster(outfile))) ## 0 1


## >> step four_mask to contain only 1 as min/max values ####
## $$ TO FIX: NoData value in mask (globalmask_ee_10k_nodata) is 0. Should this be changed to -9999 (Alejandro to decide what is best here). Changing nodata value to -9999 using gdal_cals or gdal_translate or change_values makes mask values = 0 1

infile <- outfile
outfile <- gsub("ee_10k", "ee_10k_nodata", infile)
system(sprintf("gdalwarp -dstnodata 0 %s %s", infile, outfile))
  ## Note for Alejandro: output of sprintf below
  ## "gdalwarp -dstnodata 0 /tempdata/research-cifs/6300-payalb/uom_data/gsdms_data/outputs/layers_10k/mask/globalmask_ee_10k.tif /tempdata/research-cifs/6300-payalb/uom_data/gsdms_data/outputs/layers_10k/mask/globalmask_ee_10k_nodata.tif"

gdalinfo(outfile) ## NoData Value=0
unique(values(raster(outfile))) ## NA  1

gdalinfo(file.path(output_dir, "mask", 
                   sprintf("globalmask_ee_%s_nodata.tif", proj.res.km)))
unique(values(raster(
  file.path(output_dir, "mask", 
            sprintf("globalmask_ee_%s_nodata.tif", proj.res.km)))))


## Note: Alternate mask file created using change_values function: globalmask_ee_10k_nodata9999.tif with NoData Value=-9999 and unique values of 0 1




# ## ------------------------------------------------------------------------- ##
# ## Calculate distance layers ####
# ## ------------------------------------------------------------------------- ##
#  ## Note for Alejandro: These steps are not needed UNLESS we fidn an issue with the
#  ##  files in "/tempdata/research-cifs/6300-payalb/uom_data/gsdms_data/outputs/distance_layers_wgs"
#
# dst_folder <- file.path(gsdms_dir, "outputs", "distance_layers_wgs")
# if(!dir.exists(dst_folder)) {
#   dir.create(dst_folder)
# }
# 
# 
# ## >> Builtup ####
# ## Values are expressed as decimals (Float) from 0 to 100 (see p14 in  GHSL_Data_Package_2019_light.pdf)
# infile <- file.path(gsdms_dir, "builtup", "GHS_BUILT_LDS2014_GLOBE_R2018A_54009_1K_V2_0.tif")
# gdalinfo(infile)
# gdalsrsinfo(infile)
# 
# ## >> >> Convert Projection
# new_crs <- wgs_crs
# old_crs <- "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"
# outfile <- gsub(".tif$", "_wgs.tif", infile)
# system(sprintf("gdalwarp -overwrite -ot Float32 -r bilinear -s_srs '%s' -t_srs '%s' %s %s",
#                old_crs, new_crs, infile, outfile))
# gdalinfo(outfile)
# raster(outfile)
# 
# ## >> >> Change NoData values to -9999
# system(paste0("bash ", change_values, " ", outfile, " -9999"))
# file.remove(outfile)
# file.rename(gsub(".tif", "_edt.tif", outfile),
#             outfile)
# file.remove(gsub(".tif", "_temp.tif", outfile))
# gdalinfo(outfile)
# 
# ## >> >> Calculate distance layer
# ## Note: Distance is calculated from all non 0 values
# infile <- outfile
# outfile <- file.path(dst_folder, "dst_builtup.tif")
# system(sprintf("gdal_proximity.py %s %s -of GTiff -distunits GEO -nodata %d -co compress=LZW -co BIGTIFF=YES", infile, outfile, -9999))
# gdalinfo(outfile)
# # system(sprintf("gdalinfo %s", outfile))
# 
# 
# ##  >> Roads ####
# ## >> >> The input file geodatabase
# ## https://gis.stackexchange.com/questions/151613/reading-feature-class-in-file-geodatabase-using-r
# roadsgdb <- file.path(gsdms_dir, "groads/groads-v1-global-gdb/gROADS_v1.gdb")
# 
# ## >> >> List all feature classes in a file geodatabase
# subset(ogrDrivers(), grepl("GDB", name))
# ogrListLayers(roadsgdb)
# 
# ## >> >> Read the feature class
# roads <- readOGR(dsn=roadsgdb,layer="Global_Roads")
# summary(roads)
# unique(roads$FCLASS)
# unique(roads$EXS)
# roadstemp <- roads[, 4]
# 
# ## >> >> Subset the road attributes (see gROADSv1_documentation.pdf) 
# ## 1. FClass = Functional class (1=Highway, 2=Primary, 3=Secondary, 4=Tertiary, 5=Local/ Urban, 6=Trail, 7=Private, 0=Unspecified)
# ## 2. EXS = Existence Category (1=Definite, 2=Doubtful, 0=Unspecified) ## do not use this caregory, nnorth canada and UK excluded
# # roadsub <- roads[roads$FCLASS == (0:3),] ## groads03 files
# roadsub <- roads[roads$FCLASS == (0:2),]
# unique(roadsub$FCLASS)
# unique(roadsub$EXS)
# 
# ## >> Drop columns
# roadsub <- roadsub[, 4]
# 
# ## >> >> Write shapefile
# writeOGR(roadsub, dsn = file.path(gsdms_dir, "groads"), layer = "groads", driver="ESRI Shapefile")
# 
# ## >> >> Rasterise @ .0025 ~ 250m
# system(paste0("gdal_rasterize -at -burn 1 -ot Byte -tr .0025 .0025 ",
#               file.path(gsdms_dir, "groads", "groads.shp "),
#               file.path(gsdms_dir, "groads", "groads.tif")))
# gdalinfo(file.path(gsdms_dir, "groads", "groads.tif"))
# rm(roads, roadsub, roadsgdb)
# 
# ## >> >> Calculate distance raster
# infile <- file.path(gsdms_dir, "groads", "groads.tif")
# gdalinfo(infile)
# 
# outfile <- file.path(dst_folder, "dst_roads.tif")
# 
# job::job({system(sprintf("gdal_proximity.py %s %s -values %d -of GTiff -distunits GEO -nodata %d -co compress=LZW -co BIGTIFF=YES", infile, outfile, 1, -9999))
#   gdalinfo(outfile)
# })
# 
# ## >> Population density ####
# ## The population density rasters were created by dividing the population count raster for a given target year by the land area raster (Float): see https://sedac.ciesin.columbia.edu/data/set/gpw-v4-population-density-rev11
# infile <- file.path(gsdms_dir, "population/SEDAC/gpw-v4-population-density-rev11_2020_30_sec_tif", "gpw_v4_population_density_rev11_2020_30_sec.tif")
# gdalinfo(infile)
# 
# ## >> >> Change NoData values to -9999
# outfile <- file.path(dst_folder, basename(infile))
# system(sprintf("gdalwarp -dstnodata -9999 %s %s", outfile, gsub(".tif", "2.tif", outfile)))
# 
# ## >> >> Calculate distance layer
# ## Distance is calculated from all values for pop density of 60
# infile <- outfile
# outfile <- file.path(dst_folder, "dst_pop.tif")
# system(sprintf("gdal_proximity.py %s %s -values %d -of GTiff -distunits GEO -nodata %s -co compress=LZW -co BIGTIFF=YES", infile, outfile, 60, -9999))
# gdalinfo(outfile)
# file.remove(infile)
# 
# 
# ## >> Water bodies ####
# ## Source: https://www.worldwildlife.org/publications/global-lakes-and-wetlands-database-lakes-and-wetlands-grid-level-3
# infile <- file.path(gsdms_dir, "lakesrivers/GLWD/glwd_3")
# gdalinfo(infile)
# 
# 
# ## >> >> Lakes & rivers ####
# ## >> >> Extract rivers and lakes from data
# outfile <- file.path(gsdms_dir, "lakesrivers/GLWD", "lakesrivers.tif")
# system(paste0("gdal_calc.py -A ", infile,
#               " --calc='((A==1) + (A==2) + (A==3) + (A==4) + (A==5))*1 + ( + (A==6) + (A==7) + (A==8) + (A==9) + (A==10) + (A==11) + (A==12))*0' --NoDataValue=-9999",
#               " --outfile=", outfile))
# gdalUtils::gdalinfo(outfile)
# 
# ## >> >> Caluclate distance raster
# infile <- outfile
# outfile <- file.path(dst_folder, "dst_lakesrivers.tif")
# system(sprintf("gdal_proximity.py %s %s -values %d -of GTiff -distunits GEO -nodata %d -co compress=LZW -co BIGTIFF=YES", infile, outfile, 1, -9999))
# gdalinfo(outfile)
# 
# ## >> >> Specify projection
# ## Projection of original dataset is lat/long (see p7 of GLWD_Data_Documentation.pdf)
# infile <- outfile
# outfile <- gsub(".tif$", "_wgs.tif", infile)
# system(sprintf("gdalwarp -overwrite -ot Float32 -t_srs '%s' %s %s",
#                wgs_crs, infile, outfile))
# 
# file.remove(infile)
# file.rename(outfile, infile)
# 
# ## >> >> Wetlands ####
# ## >> >> Extract wetlands from data
# infile <- file.path(gsdms_dir, "lakesrivers/GLWD/glwd_3")
# outfile <- file.path(gsdms_dir, "lakesrivers/GLWD", "wetlands.tif")
# system(paste0("gdal_calc.py -A ", infile,
#               " --calc='((A==1) + (A==2) + (A==3) + (A==4) + (A==5))*0 + ( + (A==6) + (A==7) + (A==8) + (A==9) + (A==10) + (A==11) + (A==12))*1' --NoDataValue=-9999",
#               " --outfile=", outfile))
# gdalUtils::gdalinfo(outfile)
# 
# ## >> >> Calculate distance raster
# infile <- outfile
# outfile <- file.path(dst_folder, "dst_wetlands.tif")
# system(sprintf("gdal_proximity.py %s %s -values %d -of GTiff -distunits GEO -nodata %d -co compress=LZW -co BIGTIFF=YES", infile, outfile, 1, -9999))
# gdalinfo(outfile)
# 
# ## >> >> Specify projection
# infile <- outfile
# outfile <- gsub(".tif$", "_wgs.tif", infile)
# system(sprintf("gdalwarp -overwrite -ot Float32 -t_srs '%s' %s %s",
#                wgs_crs, infile, outfile))
# 
# file.remove(infile)
# file.rename(outfile, infile)
# 
# 
# ## >> Protected areas ####
# ## Input data format: Byte (0/1 data)
# ## Output data format: Floast32
# ## Distance is calculated from all non 0 values (only 1 and 0 values in this layer)
# infile <- file.path(gsdms_dir, "protectedareas/wdpa.tif")
# gdalinfo(infile)
# outfile <- file.path(dst_folder, "dst_wdpa.tif")
# 
# system(sprintf("gdal_proximity.py %s %s -values %d -of GTiff -distunits GEO -nodata %s -co compress=LZW -co BIGTIFF=YES", infile, outfile, 1, -9999))
# gdalinfo(outfile)
# 
# 
# ## >> Check all dst layers
# list.files(file.path(gsdms_dir, "distance_layers_wgs"),
#            full.names = TRUE)
# lapply(list.files(file.path(gsdms_dir, "distance_layers_wgs"),
#                   full.names = TRUE), gdalsrsinfo)




## ------------------------------------------------------------------------- ##
## Processing environmental layers ####
## ------------------------------------------------------------------------- ##
## $$ Note for alejandro: I have updated the steps as per biofuture_aus_processing_Alejandro.R
##
## ADD STEP: Reduce file size by compression? + reducing precision to 4 decimal places 
##  [-co DECIMAL_PRECISION=4 in gdalwarp/gdal_translate does nnot work for GTiffs; 
##  alternate hacky way: multiply values by 100 and specify "Int16" 
##  as output type. Divide values by 100 when loading rasters to 
##  get the decimals again]


## >> specify file names for data ####
## >> >> Elevation layer
srtm <- file.path(gsdms_dir, "srtm/mn30_grd/srtm.adf")

## >> >> Soil layers
soil <- list.files(file.path(gsdms_dir,"soil"), pattern = "_wgs.tif", 
                   full.names = TRUE, recursive = TRUE)

## >> >> Bias (distance) layers
dst_layers <- list.files(file.path(gsdms_dir, "distance_layers_wgs"),
                         full.names = TRUE)

## >> >> Current climate (worldclim) layers
bio_current <- list.files(file.path(gsdms_dir, "worldclim", "current", "wc2.1_30s_bio"), full.names = TRUE)

## >> >> Future climate (worldclim) layers ####
ssps <- c("ssp1", "ssp3", "ssp5")
bio_future <- list.files(file.path(skip_dir, "exp_alex/out"), 
                         pattern = paste0(ssps, collapse = "|"), 
                         recursive = TRUE,
                         full.names = TRUE)
bio_future <- grep(".tif$", bio_future, value = TRUE) ## remove file types that are not .tif
bio_future_ssp1 <- grep("ssp1", bio_future, value = TRUE) ## files for ssp1 only
bio_future_ssp3 <- grep("ssp3", bio_future, value = TRUE) ## files for ssp3 only 
bio_future_ssp5 <- grep("ssp5", bio_future, value = TRUE) ## files for ssp5 only 

## >> >> Tree cover (Hansen) layer 
## $$ Note for Alejandro: This is a LARGE file, you might have to run this independently.
tree <- list.files(file.path(gsdms_dir, "treecover_hansen"), 
                   pattern = ".tif$", full.names = TRUE)

## >> >> Ignore files
  # chelsa_current <- list.files(file.path(gsdms_dir, "CHELSA/envicloud/chelsa/chelsa_V2/GLOBAL/climatologies/1981-2010/bio"), full.names = TRUE, recursive = TRUE)
  # kgvars <- grep("CHELSA_kg", chelsa_current, value = TRUE)
  # chelsa_current <- grep("CHELSA_bio", chelsa_current, value = TRUE)


## >> specify input file names #### 
## $$ Note to Alejandro: Specify a subset of files here, e.g. c(soil) or c(srtm,soil), etc.
cov_files <- c(srtm)
all(lapply(cov_files, file.exists))

if(length(cov_files) > future::availableCores()-2) {
  mc.cores = future::availableCores()-2
}else{
  mc.cores = length(cov_files)
}


## >> Step one_Check layer CRS and specify input projection IF this is missing ####
## NOTE: Only run step IF CRS is missing, skip step otherwise
infiles <- c(cov_files)
lapply(infiles, gdalsrsinfo)
new_crs = "EPSG:4326"

outfiles <- file.path(output_subdir, paste0(tools::file_path_sans_ext(basename(infiles)), "_wgs.tif"))

parallel::mclapply(seq_along(infiles),
                   function(x) system(sprintf("gdal_translate -ot Float32 -a_srs '%s' %s %s",
                                              new_crs, infiles[x], outfiles[x])),
                   mc.cores = mc.cores, mc.preschedule = TRUE)

length(list.files(output_subdir, pattern = "_wgs")) ## should be equal to number of input files
lapply(outfiles, gdalsrsinfo)


## >> Step two_Reduce layer extent, as specified in WGS 84 ####
new_extent <- "-180 -60 180 90"
# ## IF step one is run, uncomment this:
# infiles <- outfiles
# outfiles = gsub("_wgs", "_clip", infiles)

## IF step one is not run, uncomment this:
infiles <- cov_files
outfiles <- file.path(output_subdir, paste0(tools::file_path_sans_ext(basename(infiles)), "_clip.tif"))

parallel::mclapply(seq_along(infiles),
                   function(x) system(paste0("gdalwarp -overwrite -ot Float32",
                                             " -te ", new_extent, " ",
                                             infiles[x], " ", outfiles[x])),
                   mc.cores = mc.cores, mc.preschedule = TRUE)

length(list.files(output_subdir, pattern = "_clip.tif")) # should be equal to number of input files
gdalinfo(outfiles[1])


## >> Step three_Reproject layer to Equal earth ####
## resamplig method: bilinear  (https://support.esri.com/en/technical-article/000005606)
new_res <- proj_res
new_crs = equalearth_crs

infiles <- outfiles
outfiles <- gsub("_clip", "_reproj", infiles)

parallel::mclapply(seq_along(infiles),
                   function(x) system(sprintf("gdalwarp -overwrite -ot Float32 -r bilinear -tr %s -s_srs '%s' -t_srs '%s' %s %s",
                                              paste(new_res, collapse = " "),
                                              wgs_crs, new_crs,
                                              infiles[x], outfiles[x])),
                   mc.cores = mc.cores, mc.preschedule = TRUE)

length(list.files(output_subdir, pattern = "_reproj")) # should be equal to number of input files
lapply(outfiles, gdalsrsinfo)
gdalinfo(outfiles[1])
raster(outfiles[1])


## >> Step four_Clip layer s.t. #cells in layer equal #cells in mask ####
new_extent <- extent(raster(mask_file))

infiles <- sort(outfiles)
outfiles <- sub("_reproj", "_clip2", outfiles)

parallel::mclapply(seq_along(infiles),
                   function(x) system(sprintf("gdalwarp -overwrite -ot Float32 -te %s %s %s %s %s %s",
                                              new_extent@xmin,
                                              new_extent@ymin,
                                              new_extent@xmax,
                                              new_extent@ymax,
                                              infiles[x], outfiles[x])),
                   mc.cores = mc.cores, mc.preschedule = TRUE)


length(list.files(output_subdir, pattern = "_clip2")) # should be equal to number of input files
raster(infiles[1])
raster(mask_file)
raster(outfiles[1])
gdalinfo(outfiles[1])


## >> Step five_Change nodata values to -9999 ####
## To remove nodata value of -3.4e+38
## This is not needed for bio_future layers (check)
infiles <- outfiles
outfiles <- sub("_clip2", "_nodataval", outfiles)

parallel::mclapply(seq_along(infiles),
                   function(x) system(sprintf("gdalwarp -dstnodata -9999 %s %s", infiles[x], outfiles[x])),
                   mc.cores = mc.cores, mc.preschedule = TRUE)


length(list.files(output_subdir, pattern = "_nodataval")) # should be equal to number of input files
gdalinfo(infiles[1])
gdalinfo(outfiles[1])


## >> Step six_Mask ####
## https://gitlab.unimelb.edu.au/garberj/gdalutilsaddons/-/blob/master/gdal_calc.R
infiles <- outfiles
outfiles <- sub("_nodataval", "", outfiles)
parallel::mclapply(seq_along(infiles),
                   function(x) system(paste0("gdal_calc.py -A ",
                                             infiles[x], " -B ", mask_file,
                                             " --calc='((B==1)*A)+(-9999*(B!=1))' --NoDataValue=-9999",
                                             " --outfile=", outfiles[x])),
                   mc.cores = mc.cores, mc.preschedule = TRUE)

## check number of files manually
gdalUtils::gdalinfo(mask_file)
gdalUtils::gdalinfo(outfiles[1])


## >> Remove excess files ####
list.files(output_subdir, pattern = "_clip.tif")
file.remove(list.files(output_subdir, pattern = "_clip.tif", full.names = TRUE))
list.files(output_subdir, pattern = "_reproj.tif")
file.remove(list.files(output_subdir, pattern = "_reproj.tif", full.names = TRUE))
list.files(output_subdir, pattern = "_clip2.tif")
file.remove(list.files(output_subdir, pattern = "_clip2.tif", full.names = TRUE))
list.files(output_subdir, pattern = "_nodataval.tif")
file.remove(list.files(output_subdir, pattern = "_nodataval.tif", full.names = TRUE))




## ------------------------------------------------------------------------- ##
## Create seperate tifs for SRTM variables  ####
## ------------------------------------------------------------------------- ##
infile <- list.files(output_dir, pattern = "srtm_eqar.tif$", full.names = TRUE)
gdalinfo(infile)

outfile <- gsub("srtm", "slope", infile)
system(paste0("gdaldem slope ", infile, " ", outfile))

outfile <- gsub("srtm", "aspect", infile)
system(paste0("gdaldem aspect ", infile, " ", outfile))

outfile <- gsub("srtm", "terrain", infile)
system(paste0("gdaldem TRI ", infile, " ", outfile))

outfile <- gsub("srtm", "roughness", infile)
system(paste0("gdaldem roughness ", infile, " ", outfile))

file.rename(infile, gsub("srtm", "elevation", infile))








## ------------------------------------------------------------------------- ##
## $$ TO FIX - Prepare land use layers ####
## ------------------------------------------------------------------------- ## 
## Hurtt data: https://luh.umd.edu/

## $$ Note for Alejandro: IGNORE for now, we will tackle this as TASK 3 :)

## >> Future: 2015 - 2100 (86 time steps) ####
hurtt_out <- file.path(gsdms_dir, "outputs", "landuse_hurtt_future")
if(!dir.exists(hurtt_out)) {
  dir.create(hurtt_out)
}
# file.remove(list.files(hurtt_out, full.names = TRUE))

## Explore Hurtt input files
infiles = list.files(file.path(gsdms_dir, "landuse/Hurtt"), 
                     pattern = ".nc$", full.names = TRUE)

gdalUtils::gdalinfo(infiles[1])
lu <- ncdf4::nc_open(infiles[1], write=FALSE, readunlim=FALSE, 
                     verbose=FALSE,auto_GMT=TRUE, 
                     suppress_dimvals=FALSE)
names(lu)
names(lu$var)

temp <- ncvar_get(lu, varid = "primf", 
                  start = NA, count = NA, 
                  verbose=FALSE,
                  signedbyte=TRUE, 
                  collapse_degen=TRUE, 
                  raw_datavals=TRUE )
str(temp)

## Specify land use categories in data
luvars <- c("primf", "primn", "secdf", "secdn", "urban", "pastr", "range", "crop" )

## Specify time steps and index
t.all <- seq(2015, 2100, 1)
t.steps <- seq(2015, 2070, 5)
t.idx <- which(t.all %in% t.steps)
t.idx.py <- t.idx - 1 ## as per python indices


## >> >> step one_Extract layers by ssp, land use and year ####
system("bash /tempdata/workdir/gsdms/scripts/extract_nc.sh")

## Check files
length(list.files(hurtt_out))
x <- list.files(hurtt_out, full.names = TRUE)
raster(x[1])
plot(raster(x[1]))


## >> >> step two_Rename files by year of time slice ####
for (i in 1:length(t.idx.py)){
  x <- list.files(hurtt_out, full.names = TRUE, pattern = paste0("t", t.idx.py[i], ".nc$"))
  y <- gsub(paste0("t", t.idx.py[i], ".nc$"), paste0("t", t.steps[i], ".nc"), x)
  file.rename(x, y)
}


## >> >> step three_Aggregate crop layers ####
system("bash /tempdata/workdir/gsdms/scripts/max_nc.sh")


## >> >> step four_Reproject to Equal Earth ####
infile <- list.files(hurtt_out, 
                     pattern = paste0(luvars, collapse = "|"), 
                     full.names = TRUE)
length(infile) == length(luvars)*12*3 ## number of variables *  number of time steps * number of scenarios
gdalUtils::gdalinfo(infile[1])
raster(infile[1])

outfile <- file.path(output_dir, paste0("hurtt_", basename(tools::file_path_sans_ext(infile)), ".tif"))
new_res <- proj_res
new_crs = equalearth_crs

parallel::mclapply(seq_along(infile),
                   function(x) system(paste0("gdalwarp -overwrite -ot Float32 -r bilinear -tr ",
                                             paste(new_res, collapse = " "),
                                             " -s_srs '", wgs_crs, "'",
                                             " -t_srs '", new_crs, "' ",
                                             infile[x], " ", outfile[x])),
                   mc.cores = mc.cores, mc.preschedule = TRUE)
## Note warning: Warning 1: No UNIDATA NC_GLOBAL:Conventions attribute

raster(infile[1])
raster(outfile[1])
raster(file.path(output_dir, sprintf("globalmask_%sk_ee.tif", proj.res.km)))
## Check extent with maskfile


## >> >> step five_Mask layers ####
## https://gitlab.unimelb.edu.au/garberj/gdalutilsaddons/-/blob/master/gdal_calc.R
infile <- outfile
outfile <- sub(".tif", "_ee.tif", outfile)
mask_file <- file.path(output_dir, sprintf("globalmask_%sk_ee.tif", proj.res.km))

gdalUtils::gdalinfo(infile[1])
gdalUtils::gdalinfo(mask_file)
parallel::mclapply(seq_along(infile),
                   function(x) system(paste0("gdal_calc.py -A ", infile[x], " -B ", mask_file,
                                             " --calc='((B==1)*A)+(-9999*(B!=1))' --NoDataValue=-9999",
                                             " --outfile=", outfile[x])),
                   mc.cores = mc.cores, mc.preschedule = TRUE)
gdalUtils::gdalinfo(outfile[1])
file.remove(infile)



## >> Current: 2015 ####
hurtt_out <- file.path(gsdms_dir, "outputs", "landuse_hurtt_2015")
if(!dir.exists(hurtt_out)) {
  dir.create(hurtt_out)
}

# file.remove(list.files(hurtt_out, full.names = TRUE))
infile = list.files(file.path(gsdms_dir, "landuse/hurtt/LUH2/LUH2_v2h"), 
                    pattern = ".nc$", full.names = TRUE)
lu <- ncdf4::nc_open(infile, write=FALSE, readunlim=FALSE, 
                     verbose=FALSE,auto_GMT=TRUE, 
                     suppress_dimvals=FALSE)
temp <- ncvar_get(lu, varid = "primf", 
                  start = NA, count = NA, 
                  verbose=FALSE,
                  signedbyte=TRUE, 
                  collapse_degen=TRUE, 
                  raw_datavals=TRUE )
str(temp)

## Specify time steps and index
t.all <- seq(850, 2015, 1)
t.steps <- 2015
t.idx <- which(t.all %in% "2015")
t.idx.py <- t.idx - 1 ## as per python indices


## >> >> step one_Extract layers year (2015 only) ####
system("bash /tempdata/workdir/gsdms/scripts/extract_nc_hurtt_2015.sh")


## >> >> step two_Rename files by year of time slice ####
x <- list.files(hurtt_out, full.names = TRUE, pattern = paste0("t", t.idx.py, ".nc$"))
y <- gsub(paste0("t", t.idx.py, ".nc$"), paste0("t", t.steps, ".nc"), x)
file.rename(x, y)


## >> >> step three_Aggregate crop layers ####
system("bash /tempdata/workdir/gsdms/scripts/max_nc_hurtt_2015.sh")


## >> >> step four_Reproject to Equal Earth ####
infile <- list.files(hurtt_out, 
                     pattern = paste0(luvars, collapse = "|"), 
                     full.names = TRUE)
length(infile) == length(luvars) ## number of variables

outfile <- file.path(output_dir, paste0("hurtt_", basename(tools::file_path_sans_ext(infile)), ".tif"))
new_res <- proj_res

parallel::mclapply(seq_along(infile),
                   function(x) system(paste0("gdalwarp -overwrite -ot Float32 -r bilinear -tr ",
                                             paste(new_res, collapse = " "),
                                             " -s_srs '", wgs_crs, "'",
                                             " -t_srs '", new_crs, "' ",
                                             infile[x], " ", outfile[x])),
                   mc.cores = mc.cores, mc.preschedule = TRUE)
file.exists(outfile)


## >> >> step five_Mask layers ####
## https://gitlab.unimelb.edu.au/garberj/gdalutilsaddons/-/blob/master/gdal_calc.R
infile <- outfile
outfile <- sub(".tif", "_ee.tif", outfile)
mask_file <- file.path(output_dir, sprintf("globalmask_%sk_ee.tif", proj.res.km))

parallel::mclapply(seq_along(infile),
                   function(x) system(paste0("gdal_calc.py -A ", infile[x], " -B ", mask_file,
                                             " --calc='((B==1)*A)+(-9999*(B!=1))' --NoDataValue=-9999",
                                             " --outfile=", outfile[x])),
                   mc.cores = mc.cores, mc.preschedule = TRUE)
file.exists(outfile)
file.remove(infile)
