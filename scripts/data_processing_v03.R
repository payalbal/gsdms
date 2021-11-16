## Covariate data processing for global analyses
## NOTE: Data sources recorded in data_downloads.R


## --------------------- ##
## Set working environment ####
## --------------------- ##
rm(list = ls())
gc()
# system("ps")
# system("pkill -f R")

## Ensure ncdf4{} works
## In terminal ssh into ubuntu@global-diversity/or in R terminal
## type nc-config --libs
Sys.getenv()
Sys.getenv("LD_LIBRARY_PATH")
Sys.setenv(LD_LIBRARY_PATH = paste(Sys.getenv("LD_LIBRARY_PATH"), "/usr/lib/x86_64-linux-gnu", "/usr/lib/x86_64-linux-gnu/hdf5/serial", sep = ":"))

## Install processNC{}
## https://github.com/RS-eco/processNC
if(!"remotes" %in% installed.packages()[,"Package"]) install.packages("remote")
if(!"processNC" %in% installed.packages()[,"Package"]) remotes::install_github("RS-eco/processNC", build_vignettes=T)

## cdo 
## https://code.mpimet.mpg.de/projects/cdo/


## Load libraries
# devtools::install_github('skiptoniam/sense')
x <- c('data.table', 'sp', 'raster', 
       'rgdal', 'sense', 'tools', 'bitops', 
       'ncdf4', 'processNC',
       'RCurl', 'gdalUtils', 'usethis',
       'parallel', 'doMC',
       'reticulate')
lapply(x, require, character.only = TRUE)
rm(x)

## Folder paths
data_dir <- "/home/payalb/gsdms_r_vol/tempdata/research-cifs/uom_data/gsdms_data"
output_dir <- file.path(data_dir, "processed_data_10k")
if(!dir.exists(output_dir)) {
  dir.create(output_dir)
}

# ## To copy files from land-use
# lufiles <- list.files("/home/payalb/gsdms_r_vol/tempdata/workdir/landuse_projects/land-use/processed_layers/", 
#                      full.names = TRUE, all.files = TRUE)
# basename(tools::file_path_sans_ext(lufiles))
# x <- grep("landuse|lu", basename(tools::file_path_sans_ext(lufiles)), value = FALSE)
# lufiles <- lufiles[-x]
# 
# file.copy(lufiles, output_dir,
#           overwrite = FALSE, recursive = TRUE,
#           copy.mode = TRUE, copy.date = TRUE)

## Functions
source("/home/payalb/gsdms_r_vol/tempdata/workdir/gsdms/scripts/0_functions.R")
source("/home/payalb/gsdms_r_vol/tempdata/workdir/landuse_projects/land-use/gdal_calc.R") # by jgarber
# ## To pull directly from jgarber's repo:
# library(RCurl)
# url <- "https://gitlab.unimelb.edu.au/garberj/gdalutilsaddons/-/blob/master/gdal_calc.R"
# url.exists(url)
# source(url)
# devtools::source_url(url)
# url <- "https://gitlab.unimelb.edu.au/garberj/gdalutilsaddons.git"
# devtools::install_git(url = url)



## ------------------------- ##
## Create mask ####
## ------------------------- ##

## step one_create mask from WorldClim layer 
## crs for Worldclim layers: https://worldclim.org/data/v1.4/formats.html
bio_current <- list.files(paste0(data_dir, "/bio_30s"), pattern = "bio_current*", full.names = TRUE)
global_mask <- raster(bio_current[1])
global_mask[which(!is.na(global_mask[]))] <- 1
wgs_crs = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0" 
crs(global_mask) <- wgs_crs
writeRaster(global_mask, filename = file.path(output_dir, "globalmask_wgs_1degree.tif"))

## step two_clip mask to desired extent
## Note: Clipping has to be done before the next step for proj and reso
infile <- file.path(output_dir, "globalmask_wgs_1degree.tif")
outfile <- file.path(output_dir, "globalmask_clip.tif")
new_extent <- "-180 -60 180 90"
system(paste0("gdalwarp -overwrite -ot Byte",
              " -te ", new_extent, " ",
              infile, " ", outfile))
raster(infile)
raster(outfile)

## step three_reproject mask to Equal Earth and desired resolution
## Ref for EE: https://proj.org/operations/projections/eqearth.html#id2
## jgarber note the '+proj=eqearth' is a new format and wont work wth gdalUtils
## instead we use '+proj=eqearth +ellips=WGS84 +wktext' which give equal earth and works here
## see: https://en.wikipedia.org/wiki/Equal_Earth_projection
infile <- outfile
outfile <- gsub("_clip", "_ee", infile)
new_res <- c(10000,10000) ## 10km res
# new_res <- c(1000,1000) ## 1km res
new_crs = '+proj=eqearth +ellips=WGS84 +wktext'
system(paste0("gdalwarp -overwrite -ot Byte -tr ",
              paste(new_res, collapse = " "),
              " -t_srs '", new_crs, "' ",
              infile, " ", outfile))
# gdalwarp(outfile, outfile2, s_srs = wgs_crs, t_srs = new_crs, tr = reso_ee) #, output_Raster =TRUE, verbose=TRUE) #these extra arguments were used to see what was going on and the output
raster(infile)
raster(outfile)
unique(values(raster(outfile)))
plot(raster(outfile))

file.remove(infile)

## step four_change NoData values to -9999
## ref: https://gis.stackexchange.com/questions/298230/change-no-data-value-geotif-file-with-qgis-or-gdal
change_values = "/home/payalb/gsdms_r_vol/tempdata/workdir/gsdms/scripts/change_values.sh"
system(paste0("bash ", change_values, " ", outfile, " -9999"))

gdalUtils::gdalinfo(outfile)
gdalUtils::gdalinfo(file.path(output_dir, "globalmask_ee_edt.tif"))
unique(values(raster(file.path(output_dir, "globalmask_ee_edt.tif"))))
plot(raster(file.path(output_dir, "globalmask_ee_edt.tif")))

file.rename(file.path(output_dir, "globalmask_ee_edt.tif"), 
            file.path(output_dir, "globalmask_10k_ee.tif"))

file.remove(file.path(output_dir, "globalmask_ee_temp.tif"))
file.remove(outfile)




## ----------------------------------------------------------------------------- ##
## Processing CATEGORICAL covariate layers (land use) - Equal Area projection ####
## ----------------------------------------------------------------------------- ## 

## Land use 
## >> I. Source: ESA (fractional land use) ####
## http://maps.elie.ucl.ac.be/CCI/viewer/download.php 
landuse <- list.files("/home/payalb/gsdms_r_vol/tempdata/research-cifs/uom_data/gsdms_data/landuse/CCI",
                      pattern = "landuse.tif$", full.names = TRUE)
outfile <- gsub("ESA_landuse", "ESA_landuse_reclass", landuse)

## >> >> step one_reclassify ####
## Recalssify according to flutes classes - See landuse_classifications.xlsx
gdalUtils::gdalinfo(landuse)
system(paste0("gdal_calc.py -A ", landuse,
              " --calc='((A==10) + (A==11) + (A==12) + (A==20) + (A==30))*1 + ((A==40))*2 + ((A==50) + (A==60) + (A==61) + (A==62) + (A==70) + (A==71) + (A==72) + (A==80) + (A==81) + (A==82) + (A==90) + (A==100) + (A==160) + (A==170))*3 + ((A==110) + (A==130))*4 + ((A==120) + (A==121) + (A==122))*5 + ((A==180))*6 + ((A==190))*7 + ((A==140) + (A==150) + (A==151) + (A==152) + (A==153) + (A==200) + (A==201) + (A==202) + (A==220))*8 + ((A==210))*9' --NoDataValue=0",
              " --outfile=", outfile))
gdalUtils::gdalinfo(outfile)


## >> >> step two_change NoData values to -9999 ####
change_values = "/home/payalb/gsdms_r_vol/tempdata/workdir/gsdms/scripts/change_values.sh"
system(paste0("bash ", change_values, " ", outfile, " -9999"))
gdalUtils::gdalinfo(gsub(".tif", "_edt.tif", outfile))
file.copy(gsub(".tif", "_edt.tif", outfile), file.path(output_dir, "ESA_landuse_reclass.tif"))
file.remove(c(gsub(".tif", "_edt.tif", outfile), gsub(".tif", "_temp.tif", outfile)))


## >> >> step three_create separate tif for each landuse class ####
## NOTES: get unique values: https://gis.stackexchange.com/questions/33388/python-gdal-get-unique-values-in-discrete-valued-raster
## NOTES: run python code: https://rstudio.github.io/reticulate/
infile <- file.path(output_dir, "ESA_landuse_reclass.tif")
vals <- 1:9
class <- c("crop", "crop_mosaic", "forest", "grass", 
           "shrub", "wetland", "urban", "other", "water")
outfile <- file.path(output_dir, paste0("lu", vals, class, ".tif"))

parallel::mclapply(seq_along(vals),
                   function(x) system(paste0("gdal_calc.py -A ", infile,
                                             " --calc='((A==", x, "))*1 + (",
                                             paste0("(A==", vals[-x], ")",
                                                    collapse = " + "),
                                             ")*0' --NoDataValue=-9999",
                                             " --outfile=", outfile[x])),
                   mc.cores = length(vals), mc.preschedule = TRUE)
gdalUtils::gdalinfo(outfile[1])


## >> >> step four_clip by e ####
infile <- outfile
outfile <- gsub(".tif", "_clip.tif", infile)
outfile = gsub("\\..*", ".tif", outfile)
new_extent <- "-180 -60 180 90"

parallel::mclapply(seq_along(vals),
                   function(x) system(paste0("gdalwarp -overwrite -ot Byte",
                                             " -te ", new_extent, " ",
                                             infile[x], " ", outfile[x])),
                   mc.cores = length(vals), mc.preschedule = TRUE)

# system(paste0("gdalwarp -overwrite -ot Byte",
#               " -te ", new_extent, " ",
#               infile, " ", outfile))
gdalUtils::gdalinfo(outfile[1])



## >> >> step five_create fractional layers at 0.1 degrees (= ~10k) ####
## >> >> -- NOT WORKING -- ####
infile <- outfile
outfile <- gsub("_clip.tif", "_frac.tif", infile)
new_res <- c(0.1, 0.1)

system(paste0("gdalwarp -overwrite -r average -tr ",
              paste(new_res, collapse = " "), " ",
              infile[1], " ", outfile[1]))


## >> >> step six_reproject to Equal Earth ####
## resamplig method: near (https://support.esri.com/en/technical-article/000005606)
## see: https://gis.stackexchange.com/questions/352476/bilinear-resampling-with-gdal-leaves-holes
infile <- outfile #list.files(output_dir, pattern = "_clip", full.names = TRUE)
outfile <- gsub("_clip.tif", "_frac.tif", infile)
new_res <- c(10000,10000)
new_crs = '+proj=eqearth +ellips=WGS84 +wktext'
wgs_crs = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0" 

# parallel::mclapply(seq_along(vals),
#                    function(x) system(paste0("gdalwarp -overwrite -ot Byte -r average -tr ",
#                                              paste(new_res, collapse = " "),
#                                              " -s_srs '", wgs_crs, "'",
#                                              " -t_srs '", new_crs, "' ",
#                                              infile[x], " ", outfile[x])),
#                    mc.cores = length(vals), mc.preschedule = TRUE)

system(paste0("gdalwarp -overwrite -r average -tr ",
              paste(new_res, collapse = " "),
              " -t_srs '", new_crs, "' ",
              infile[1], " ", outfile[1]))

gdalUtils::gdalinfo(infile[1])
sort(unique(values(raster(outfile[1]))))
range(unique(values(raster(outfile[1]))), na.rm = TRUE)
file.remove(infile)

# ## >> >> step_clip: to make #cells equal between mask and layer ####
# infile <- outfile
# outfile <- sub("_ee", "_clip2", outfile)
# mask_file <- file.path(output_dir, "globalmask_10k_ee.tif")
# new_extent <- extent(raster(mask_file))
# 
# system(paste0("gdalwarp -overwrite -ot Byte -te ",
#               paste(new_extent[1], new_extent[3],
#                     new_extent[2], new_extent[4]), " ",
#               infile, " ", outfile))
# file.remove(infile)


## >> >> step seven_mask ####
## https://gitlab.unimelb.edu.au/garberj/gdalutilsaddons/-/blob/master/gdal_calc.R
infile <- outfile
outfile <- sub("_ee", "_eqar", outfile)
mask_file <- file.path(output_dir, "globalmask_10k_ee.tif")

gdalUtils::gdalinfo(mask_file)
parallel::mclapply(seq_along(vals),
                   function(x) system(paste0("gdal_calc.py -A ", infile[x], 
                                             " -B ", mask_file,
                                             " --calc='((B==1)*A) + (-9999*(B!=1))' --NoDataValue=-9999",
                                             " --outfile=", outfile[x])),
                   mc.cores = length(vals), mc.preschedule = TRUE)
gdalUtils::gdalinfo(outfile[3])
file.remove(infile)



# ## >> II. Source: Copernicus (fractional land use) ####
# ## Link @ GEE: https://developers.google.com/earth-engine/datasets/catalog/COPERNICUS_Landcover_100m_Proba-V_Global
# ## Direct download link: https://zenodo.org/communities/copernicus-land-cover/search?page=1&size=20
# ## Data processing for global layer @ GEE co-authored by MC Loor: https://code.earthengine.google.com/?scriptPath=users%2Fpayalbal%2Fglobal_layers%3Afraclu_sklu. See processing in GEE for reclassification scheme.
# 
# ## Land use classes used:
# ## 1    urban
# ## 2    crop
# ## 3    forest
# ## 4    grass
# ## 5    other
# ## 6    NA
# 
# # ... gee.py script
# # ...download by tile
# 
# ## Stiching tiles by Maria del Mar Quiroga
# ## ref: https://stackoverflow.com/a/50235578
# 
# ## step one_Get a list of all tif tiles downloaded from Google Earth Engine
# tile_files <- list.files(path = file.path(data_dir, "copernicus/global_frac/sklu_classes/tifs/"), pattern = "*.tif", full.names = TRUE)
# 
# ## step two_Build a virtual raster file stitching all tiles
# ## WARNING: This will fail if the file it is trying to write to (output.vrt) already exists
# gdalbuildvrt(gdalfile = tile_files, output.vrt = file.path(data_dir, "copernicus","lu_world.vrt"))
# 
# ## step three_Copy the virtual raster to an actual physical file
# ## WARNING: This takes ~5 minutes to run
# gdal_translate(src_dataset = file.path(data_dir, "copernicus", "lu_world.vrt"), 
#                dst_dataset = file.path(data_dir, "copernicus", "lu_world.tif"), 
#                output_Raster = FALSE,
#                options = c("BIGTIFF=YES", "COMPRESSION=LZW"))
# 
# ## step four_Save each band (i.e. land use class) as a tif file
# landuse <- file.path(file.path(data_dir, "copernicus", "lu_world.tif"))
# landuse <- brick(landuse)
# for (i in 1:nlayers(landuse)){
#   temp <- landuse[[i]]
#   writeRaster(temp, filename = file.path(data_dir, "copernicus", paste0("landuse_class", i, ".tif")))
# }
# rm(temp, landuse)
# ## landuse classes: urban, crop, forest, grass, other
# lu1 <- file.path(file.path(data_dir, "copernicus", "landuse_class1.tif"))
# lu2 <- file.path(file.path(data_dir, "copernicus", "landuse_class2.tif"))
# lu3 <- file.path(file.path(data_dir, "copernicus", "landuse_class3.tif"))
# lu4 <- file.path(file.path(data_dir, "copernicus", "landuse_class4.tif"))
# lu5 <- file.path(file.path(data_dir, "copernicus", "landuse_class5.tif"))
# landuse <- c(lu1, lu2, lu3, lu4, lu5)



## >> III. Source: LUH2 Hurtt data ####
##

hurtt_out <- file.path(data_dir, "landuse/hurtt_processing")
if(!dir.exists(output_dir)) {
  dir.create(output_dir)
}
infile <- "/home/payalb/gsdms_r_vol/tempdata/research-cifs/uom_data/gsdms_data/landuse/Hurtt/ssp5-rcp85_2015-2100.nc"
gdalUtils::gdalinfo(infile)

lu <- nc_open(infile, write=FALSE, readunlim=FALSE, 
              verbose=FALSE,auto_GMT=TRUE, 
              suppress_dimvals=FALSE)

names(lu)
names(lu$var)
names(lu$var$primf)
str(lu$var$primf)

lu$var$primf$missval
ncvar_change_missval(primf.86, ...) 
## Remember we have to open the file as writable to be able to change
## the missing value on disk!

primf.86 <- ncvar_get(lu, varid = "primf", 
                   start = NA, count = NA, 
                   verbose=FALSE,
                   signedbyte=TRUE, 
                   collapse_degen=TRUE, 
                   raw_datavals=TRUE )

typeof(primf.86); class(primf.86); str(primf.86); sum(is.na(primf.86))
dim(primf.86)

primf <- raster::brick(primf.86)

## Change missign values to NA
primf[[1]]
primf[primf == lu$var$primf$missval] <- NA
primf[[1]]

## Set extent
extent(primf) <-  c(-180, 180, -90, 90)

## Rotate rasters
primf_flip <- t(flip(primf_reproj[[1]], direction='y' ))
plot(primf_flip)

## Reproject rasters
crs(primf) <- "+proj=longlat +datum=WGS84 +no_defs"
primf_reproj <- projectRaster(primf, crs = "+proj=eqearth +ellips=WGS84 +wktext")

plot(primf[[1]])
plot(primf_reproj[[1]])

mask_file <- file.path(output_dir, "globalmask_wgs_1degree.tif")
raster(mask_file)

primf.2015 <- primf[[1]]
primf.2100 <- primf[[86]]

r <- raster(nrows = 1440, ncols = 720,
            xmn=-180, xmx=180, ymn=-90, ymx=90,
            crs = "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0", 
            resolution = c(0.25, 0.25), 
            vals=NULL)

r[] <- getValues(raster(primf.86[ , ,1]))



## Cdo operations...
system(paste0("cdo -sinfo ", infile))
system(paste0("cdo -showname ", infile))
system(paste0("cdo -ntime ", infile))
system(paste0("cdo -griddes ", infile))

outfile <- file.path(hurtt_out, "temp.nc")
system(paste0("cdo -chname, lat_bounds, bounds_lat ", infile, " ", outfile))

system(paste0("ncrename -d time_bnds,bounds_time -d lat_bounds,bounds_lat -d lon_bounds,bounds_lon -v time_bnds,bounds_time -v lat_bounds,bounds_lat -v lon_bounds,bounds_lon ", infile, " ", outfile))


## --------------------------------------------------------------------------- ##
## Processing CONTINUOUS covariate layers - Equal Area projection ####
## --------------------------------------------------------------------------- ##

## SRTM
srtm <- file.path(data_dir, "srtm/mn30_grd/srtm.adf")

## Soil
soil <- list.files(file.path(data_dir,"orders"), pattern = "*.dat", 
                   full.names = TRUE, recursive = TRUE)

## WorldClim (current) data
bio_current <- list.files(paste0(data_dir, "/bio_30s"), pattern = "bio_current*", full.names = TRUE)

## GCM data - One GCM: BCC-CSM1-1
bio_rcp45 <- list.files(file.path(data_dir, "gcm_30s"), pattern = "*bc45*", full.names = TRUE)
bio_rcp85 <- list.files(file.path(data_dir, "gcm_30s"), pattern = "*bc85*", full.names = TRUE)
# bio_rcp60 <- list.files(file.path(data_dir, "gcm_30s"), pattern = "*bc60*", full.names = TRUE)

## GCM data quartiles - TO DO...
# ## ALL GCMs
# ## Source: from regSSP scripts
# ## *** replace with gdal functions...
# rcps <- c("45", "60", "85")
# models <- c("BC", "CC", "GS", "HD", "HE", "IP", "MI", "MR", "MC", "MG", "NO")
# 
# mask_file <- file.path(data_dir, "global_mask_ee.tif")
# gcm_files <- list.files(file.path(data_dir, 'gcm_30s'), full.names = T, recursive = T)
# gcm_masked_path <- file.path(processed_dir, 'gcm_masked')
# if(!dir.exists(gcm_masked_path)){dir.create(gcm_masked_path)}
# 
# 
# ## Mask data & stack by GCM
# global_mask <- raster(mask_file)
# for(model_i in models){
#   mod_stack <- list()
#   file_mod <- files_gcm[grepl(model_i, files_gcm)]
#   for(j in 1:length(rcps)){
#     file_mod_rcp <- file_mod[grepl(rcps[j], file_mod)]
#     rcp_stack <- list()
#     for(f in 1:length(file_mod_rcp)){
#       rcp_stack[[f]] <- mask(crop(raster(file_mod_rcp[f]), global_mask), global_mask)
#     }
#     mod_stack[[j]] <- readAll(brick(rcp_stack))
#   }
#   saveRDS(mod_stack, file = paste0(gcm_reg_path, "/", model_i, ".rds"))
# }
# rm(rcp_stack, mod_stack)
# 
# ## Extract cell-wise quartiles across GCM
# quartiles <- c("q1", "q2", "q3")
# gcm_files <- list.files(gcm_masked_path, full.names = T)
# gcm_quant_path <- file.path(processed_dir, 'gcm_quant')
# if(!dir.exists(gcm_quant_path)){dir.create(gcm_quant_path)}
# 
# inds <- which(!is.na(global_mask[]))
# 
# for(j in 1:length(rcps)){
#   saveRDS(stack(), file = paste0(gcm_quant_path, "/bio", "q1_", rcps[j], ".rds"))
#   saveRDS(stack(), file = paste0(gcm_quant_path, "/bio", "q2_", rcps[j], ".rds"))
#   saveRDS(stack(), file = paste0(gcm_quant_path, "/bio", "q3_", rcps[j], ".rds"))
#   print(paste0("processing rcp", rcps[j]))
#   for(k in 1:19){
#     print(paste0("processing bioclim var: ", k))
#     bio <- stack()
#     for(i in 1:length(models)){
#       print(paste0("processing model: ", i))
#       dat <- readRDS(gcm_files[[i]])[[j]]
#       bio <- stack(bio, dat[[k]])
#     }
#     
#     print(paste0("getting quartiles..."))
#     df1 <- na.omit(as.matrix(getValues(bio)))
#     c <-rowQuantiles(df1, probs = c(0.25, 0.5, 0.75))
#     for(m in 1:3){
#       bioclim <- readRDS(file = paste0(gcm_quant_path, "/bio", quartiles[m], "_", rcps[j], ".rds"))
#       r[inds] <- c[,m]
#       names(r) <- paste0("bio", k)
#       saveRDS(readAll(stack(bioclim, r)), file = paste0(gcm_quant_path, "/bio", quartiles[m], "_", rcps[j], ".rds"))
#     }
#   }
# }
# unlink(gcm_masked_path, recursive=T)


cov_files <- c(bio_current, bio_rcp45, bio_rcp85, srtm, soil)
all(lapply(cov_files, file.exists))

## Check SRS for files
lapply(cov_files, gdalsrsinfo, as.CRS = TRUE)
## shows ERROR due to blank CRS for some files...

## >> step one_clip by e ####
infile = cov_files
outfile = file.path(output_dir, paste0(tools::file_path_sans_ext(basename(infile)), "_clip.", tools::file_ext(infile)))
outfile = gsub("\\..*", ".tif", outfile)
new_extent <- "-180 -60 180 90"

parallel::mclapply(seq_along(infile),
                   function(x) system(paste0("gdalwarp -overwrite -ot Byte",
                                             " -te ", new_extent, " ",
                                             infile[x], " ", outfile[x])),
                   mc.cores = mc.cores, mc.preschedule = TRUE)


## >> step two_reproject: Equal Earth ####
## resamplig method: bilinear  (https://support.esri.com/en/technical-article/000005606)
infile <- outfile #list.files(output_dir, pattern = "_clip", full.names = TRUE)
outfile <- gsub("_clip", "_ee", infile)
new_res <- c(10000,10000)
new_crs = '+proj=eqearth +ellips=WGS84 +wktext'
wgs_crs = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0" 

parallel::mclapply(seq_along(infile),
                   function(x) system(paste0("gdalwarp -overwrite -ot Byte -r bilinear -tr ",
                                             paste(new_res, collapse = " "),
                                             " -s_srs '", wgs_crs, "'",
                                             " -t_srs '", new_crs, "' ",
                                             infile[x], " ", outfile[x])),
                   mc.cores = mc.cores, mc.preschedule = TRUE)

## Checks 
n=62
raster(infile[n])
raster(outfile[n])

file.remove(infile)

## >> step three_clip: to make #cells equal between mask and layer ####
infile <- outfile
outfile <- sub("_ee", "_clip2", outfile)
mask_file <- file.path(output_dir, "globalmask_10k_ee.tif")
new_extent <- extent(raster(mask_file))

parallel::mclapply(seq_along(infile),
                   function(x) system(paste0("gdalwarp -overwrite -ot Byte -te ",
                                             paste(new_extent[1], new_extent[3],
                                                   new_extent[2], new_extent[4]), " ",
                                             infile[x], " ", outfile[x])),
                   mc.cores = mc.cores, mc.preschedule = TRUE)

## Checks
raster(infile[n])
raster(outfile[n])
raster(mask_file)

file.remove(infile)

## >> step four_mask ####
## https://gitlab.unimelb.edu.au/garberj/gdalutilsaddons/-/blob/master/gdal_calc.R
infile <- outfile
outfile <- sub("_clip2", "_eqar", outfile)
mask_file <- file.path(output_dir, "globalmask_10k_ee.tif")

gdalUtils::gdalinfo(mask_file)
parallel::mclapply(seq_along(infile),
                   function(x) system(paste0("gdal_calc.py -A ", infile[x], " -B ", mask_file,
                                             " --calc='((B==1)*A)+(-9999*(B!=1))' --NoDataValue=-9999",
                                             " --outfile=", outfile[x])),
                   mc.cores = mc.cores, mc.preschedule = TRUE)
gdalUtils::gdalinfo(outfile[n])

## Checks
raster(infile[n])
raster(outfile[n])
unique(values(raster(outfile[n])))
plot(raster(outfile[n]))

file.remove(infile)



## -------------------------------------------------------------------------- ##
## Find min non-NA set values across mask and covariates and sync NAs ####
## -------------------------------------------------------------------------- ##

mask_file <- file.path(output_dir, "globalmask_10k_ee.tif")
mask_file2 <- file.path(output_dir, "globalmask_10k_ee_minNA.tif")
file.copy(mask_file, mask_file2)

## Update mask based on covariate layers for syncing NAs
for(j in 1:length(infile)){
  
  print(paste0("Updating mask using file = ", j, " of ", length(infile) ,":", infile[j]))
  
  # gdalcalc(calc="((A>-9999) & (B==1))", infile = list(A=infile[j], B=mask_file2),outfile = mask_file2,
  #          NoDataValue=-9999, overwrite=TRUE)
  
  system(paste0("gdal_calc.py -A ", infile[j], " -B ", mask_file2,
                " --calc='((A>-9999) & (B==1))' --NoDataValue=-9999",
                " --outfile=", mask_file2))
}

## Update covariate layers based on updated mask for syncing NAs
for(j in 1:length(infile)){
  
  print(paste0("Masking file = ", j, " of ", length(infile) ,":", infile[j], " using updated mask"))
  
  gdalmask(infile = infile[j], mask = mask_file2, outfile = infile[j], output_Raster = FALSE, overwrite=TRUE, verbose=TRUE)
  
}

## Checks
freq(raster(mask_file))
freq(raster(mask_file2))
## this might remove (if applicable) some ones from the original mask and made them zeros. From looking at the map, it seems like it is mostly picking up on large lakes? Looks like the great lakes, Lake victoria, and Lake Baikal are some of the removed cells.
gdalinfo(mask_file, stats = TRUE)
gdalinfo(mask_file2, stats = TRUE)

lapply(infile, gdalinfo, stats = TRUE)
## check size for all
lapply(infile, gdalsrsinfo, as.CRS = TRUE)

n = 63
gdalinfo(outfile[n])
raster(outfile[n])
unique(values(raster(outfile[n])))
plot(raster(outfile[n]), main = tools::file_path_sans_ext(basename(outfile[n])))



## ------------------------------------------ ##
## Create seperate tif for SRTM variables  ####
## ----------------------------------------- ##

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




## --------------------------------------------------------- ##
## Create covariate data.tables for model and prediction ####
## -------------------------------------------------------- ##

infile <- list.files(output_dir, pattern = "_eqar(.*).tif$", full.names = TRUE)
## matching everything with _eqar and ending in .tif but also anythign in between the two using (.*)
infile <- infile[-grep('ESA_landuse_reclass_eqar.tif$', infile)]
## remove ESA_landuse_reclass_eqar.tif
lapply(infile, gdalinfo, stats = TRUE)
## check size for all


## >> Generate quadrature (background) points - TO MODIFY ####
## ***** Need to swap out rasterToPoints() and stack() when working with large rasters *****
## NOTE: Trial alternative approach using biogeoregions

mask_file <- file.path(input_dir, "globalmask_10k_ee_minNA.tif")
global_mask0 <- raster(mask_file); freq(global_mask0)
global_mask0[which(is.na(global_mask0[]))] <- 0; freq(global_mask0)
rpts <- rasterToPoints(global_mask0, spatial=TRUE)
mask_pts <- data.frame(rpts@data, X=coordinates(rpts)[,1], Y=coordinates(rpts)[,2])
mask_pts <- mask_pts[,-1]
fwrite(mask_pts, file = file.path(output_dir, "globalmask_10k_eeXY.csv"))
# mask_pts <- fread(file.path(output_dir, "globalmask_10k_eeXY.csv"))
dat <- stack(infile[-grep('bc45|bc85', infile)])
dat <- as.data.table(na.omit(cbind(mask_pts, as.matrix(dat))))


## >> Test for correlations in covariates ####
preds <- colnames(correlations(dat, thresh = 0.7, N = 200000))
preds <- sub("_current_", "", preds)
write.csv(preds, file.path(output_dir, "uncorrelated_covs.csv"), 
          row.names = FALSE)
# preds <- fread(file.path(output_dir, "uncorrelated_covs.csv"))$x

## >> Prepare background point data for model + remove correlated covariates ####
lucovs <- grep("landuse", names(dat), value = TRUE)
names(dat) <- sub("_current_", "", names(dat))
all(names(dat[,names(dat)[names(dat) %in% c(preds, lucovs, "X", "Y")], with = FALSE]) == preds)

dat <- dat[,names(dat)[names(dat) %in% c(preds, lucovs, "X", "Y")], with = FALSE] ## to ensure landuse, X, Y are retained
dim(dat)
names(dat) <- sub("_eqar", "", names(dat))
names(dat) <- sub("ESA_landuse_reclass_", "", names(dat))

## Rearrange covariate data accroding to static, dynamic and land use vars
## NOTE: Can be removed once we've figured out variable selection for models
static <- names(dat)[-grep("lu|bio", names(dat))]
dynamic <- names(dat)[grep("bio", names(dat))]
landuse <- names(dat)[grep("lu", names(dat))]
arrange_cols <- c(static, dynamic, landuse)
rm(static, dynamic, landuse)

dat <- dat[,..arrange_cols]
names(dat)
fwrite(dat, file = file.path(output_dir, "covariates_model.csv"))

dat_mod <- data.table::copy(dat)
rm(dat)

## >> Prepare prediction points data + remove correlated covariates ####
## RCP 45
dat <- stack(infile[-grep('bio_current|bc85', infile)])
dat <- as.data.table(na.omit(cbind(mask_pts, as.matrix(dat))))
names(dat) <- sub("bc45bi70", "bio", names(dat))
dat <- dat[,names(dat)[names(dat) %in% c(preds, lucovs, "X", "Y")], with = FALSE] ## retain landuse, X, Y 
dim(dat)
names(dat) <- sub("_eqar", "", names(dat))
names(dat) <- sub("ESA_landuse_reclass_", "", names(dat))
dat <- dat[,..arrange_cols]
all(names(dat) == names(dat_mod))

fwrite(dat, file = file.path(output_dir, "covariates_predict_rcp45.csv"))

dat_rcp45 <- data.table::copy(dat)
rm(dat)

## RCP 85
dat <- stack(infile[-grep('bio_current|bc45', infile)])
dat <- as.data.table(na.omit(cbind(mask_pts, as.matrix(dat))))
names(dat) <- sub("bc85bi70", "bio", names(dat))
dat <- dat[,names(dat)[names(dat) %in% c(preds, lucovs, "X", "Y")], with = FALSE] ## retain landuse, X, Y 
dim(dat)
names(dat) <- sub("_eqar", "", names(dat))
names(dat) <- sub("ESA_landuse_reclass_", "", names(dat))
dat <- dat[,..arrange_cols]
all(names(dat) == names(dat_mod))

fwrite(dat, file = file.path(output_dir, "covariates_predict_rcp85.csv"))

dat_rcp85 <- data.table::copy(dat)
rm(dat)

## Checks
dim(dat_mod); dim(dat_rcp45); dim(dat_rcp85)



