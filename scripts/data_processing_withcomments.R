## Data processing for global analyses
## NOTE: Data sources recorded in data_downnloads.R


## Set working environment ####
setwd("./land-use")

# devtools::install_github('skiptoniam/sense')
x <- c('sp', 'raster', 'rgdal', 'sense', 'tools', 'bitops', 'RCurl', 'gdalUtils', 'usethis')
lapply(x, require, character.only = TRUE)
source("./0_functions_par.R")
source("./gdal_calc.R") # by jgarber
# ## How to pull directly from jgarber's repo:
# library(RCurl)
# url <- "https://gitlab.unimelb.edu.au/garberj/gdalutilsaddons/-/blob/master/gdal_calc.R"
# url.exists(url)
# source(url)
# devtools::source_url(url)
# url <- "https://gitlab.unimelb.edu.au/garberj/gdalutilsaddons.git"
# devtools::install_git(url = url)

data_path <- "/tempdata/workdir/data" ## "/Volumes/discovery_data/gsdms_data"
data_processed <- file.path(getwd(), "processed_layers") 
dir.create(data_processed)


## WorldClim (current) data ####
## Source: https://www.worldclim.org/version1
bio_current <- list.files(paste0(data_path, "/bio_30s"), pattern = "bio_current*", full.names = TRUE)

## Create mask ####
## step one_create mask from WorldClim layer 
## crs for Worldclim layers: https://worldclim.org/data/v1.4/formats.html
global_mask <- raster(bio_current[1])
global_mask[which(!is.na(global_mask[]))] <- 1
wgs_crs = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0" 
crs(global_mask) <- wgs_crs
writeRaster(global_mask, filename = file.path(data_path, "global_mask.tif"))
plot(global_mask)

## step two_clip mask to desired extent
infile <- file.path(data_path, "global_mask.tif")
outfile <- file.path(data_path, "global_mask_crop2.tif")
reso_wgs <- res(global_mask)
e <- c(-180,180,-60,90)
sense::gdalCrop(inpath = infile, outpath = outfile, extent=e, resolution=reso_wgs, return = FALSE)

## step three_reproject mask to Equal Earth
## Ref for EE: https://proj.org/operations/projections/eqearth.html#id2
## jgarber note the '+proj=eqearth' is a new format and wont work wth gdalUtils
## instead we use '+proj=eqearth +ellips=WGS84 +wktext' which give equal earth and works here
## see: https://en.wikipedia.org/wiki/Equal_Earth_projection
reso_ee <- c(1000,1000)
new_crs = '+proj=eqearth +ellips=WGS84 +wktext'
outfile2 <- sub("_crop","_ee", outfile)
gdalwarp(outfile, outfile2, s_srs = wgs_crs, t_srs = new_crs, tr = reso_ee) #, output_Raster =TRUE, verbose=TRUE) #these extra arguments were used to see what was going on and the output


## GCM data ####
## One GCM: BCC-CSM1-1
bio_rcp45 <- list.files(file.path(data_path, "gcm_30s"), pattern = "*bc45*", full.names = TRUE)
bio_rcp60 <- list.files(file.path(data_path, "gcm_30s"), pattern = "*bc60*", full.names = TRUE)
bio_rcp85 <- list.files(file.path(data_path, "gcm_30s"), pattern = "*bc85*", full.names = TRUE)

## GCM data quartiles ####
# ## ALL GCMs
# ## Source: from regSSP scripts
# ## *** replace with gdal functions...
# rcps <- c("45", "60", "85")
# models <- c("BC", "CC", "GS", "HD", "HE", "IP", "MI", "MR", "MC", "MG", "NO")
# 
# mask_file <- file.path(data_path, "global_mask_ee.tif")
# gcm_files <- list.files(file.path(data_path, 'gcm_30s'), full.names = T, recursive = T)
# gcm_masked_path <- file.path(data_processed, 'gcm_masked')
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
# gcm_quant_path <- file.path(data_processed, 'gcm_quant')
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


## SRTM ####
srtm <- file.path(data_path, "srtm/mn30_grd/srtm.adf")

## Soil
soil <- list.files(file.path(data_path,"orders"), pattern = "*.dat", full.names = T, recursive = T)

## Landuse ####
## Source: Copernicus data (for fraclu analysis)
## Link @ GEE: https://developers.google.com/earth-engine/datasets/catalog/COPERNICUS_Landcover_100m_Proba-V_Global
## Direct download link: https://zenodo.org/communities/copernicus-land-cover/search?page=1&size=20
## Data processing for global layer @ GEE by MC Loor: https://code.earthengine.google.com/?scriptPath=users%2Fpayalbal%2Fglobal_layers%3Afraclu_mcloor_sklu. See processing in GEE for reclassification scheme.

## Land use classes used:
## 1    urban
## 2    crop
## 3    forest
## 4    grass
## 5    other
## 6    NA

# ... gee.py script
# ...download by title

## Stiching tiles by Maria del Mar Quiroga
## ref: https://stackoverflow.com/a/50235578

## step one_Get a list of all tif tiles downloaded from Google Earth Engine
tile_files <- list.files(path = file.path(data_path, "copernicus/tifs/"), pattern = "*.tif", full.names = TRUE)

## step two_Build a virtual raster file stitching all tiles
## WARNING: This will fail if the file it is trying to write to (output.vrt) already exists
gdalbuildvrt(gdalfile = tile_files, output.vrt = file.path(data_path, "copernicus","lu_world.vrt"))

## step three_Copy the virtual raster to an actual physical file
## WARNING: This takes ~5 minutes to run
gdal_translate(src_dataset = file.path(data_path, "copernicus", "lu_world.vrt"), 
               dst_dataset = file.path(data_path, "copernicus", "lu_world.tif"), 
               output_Raster = FALSE,
               options = c("BIGTIFF=YES", "COMPRESSION=LZW"))

## step four_Save each band (i.e. land use class) as a tif file
landuse <- file.path(file.path(data_path, "copernicus", "lu_world.tif"))
landuse <- brick(landuse)
for (i in 1:nlayers(landuse)){
  temp <- landuse[[i]]
  writeRaster(temp, filename = file.path(data_path, "copernicus", paste0("landuse_class", i, ".tif")))
}
rm(temp, landuse)
## landuse classes: urban, crop, forest, grass, other
lu1 <- file.path(file.path(data_path, "copernicus", "landuse_class1.tif"))
lu2 <- file.path(file.path(data_path, "copernicus", "landuse_class2.tif"))
lu3 <- file.path(file.path(data_path, "copernicus", "landuse_class3.tif"))
lu4 <- file.path(file.path(data_path, "copernicus", "landuse_class4.tif"))
lu5 <- file.path(file.path(data_path, "copernicus", "landuse_class5.tif"))
landuse <- c(lu1, lu2, lu3, lu4, lu5)

## Processing covariate layers ####
cov_files <- c(bio_current, bio_rcp45, bio_rcp60, bio_rcp85, srtm, soil, landuse)
all(lapply(cov_files, file.exists))

## step one_clip by e
infile <- cov_files
# file_in <- c(bioclim, srtm, soil, landuse)
outfile <- file.path(data_processed, paste0(tools::file_path_sans_ext(basename(infile)), "_cropped.", tools::file_ext(infile)))
reso_wgs <- res(raster(file.path(data_path, "global_mask_crop.tif")))
e <- c(-180,180,-60,90)
mapply(gdalCrop, inpath = infile, outpath = outfile, MoreArgs = list(extent=e, resolution=reso_wgs, return = FALSE)) 

## step two_reproject: Equal Earth
infile <- outfile #list.files(data_processed, pattern = "_cropped", full.names = TRUE)
outfile <- sub("_cropped", "_aligned", infile)
reso_ee <- c(1000,1000)
wgs_crs = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0" 
new_crs = '+proj=eqearth +ellips=WGS84 +wktext'
# gdalwarp(infile, outfile, s_srs = wgs_crs, t_srs = new_crs, tr = reso_ee) #, output_Raster =TRUE, verbose=TRUE)
mapply(gdalwarp, srcfile = infile, dstfile = outfile, MoreArgs = list(s_srs = wgs_crs, t_srs = new_crs, tr = reso_ee, verbose=TRUE))

## step three_mask
## https://gitlab.unimelb.edu.au/garberj/gdalutilsaddons/-/blob/master/gdal_calc.R
outfile2 <- sub("_aligned", "_masked", outfile)
outfile2 <- paste0(tools::file_path_sans_ext(outfile2), ".tif")

mask_file <- file.path(data_path, "global_mask_ee.tif")
mapply(gdalmask, infile = outfile, mask = mask_file, outfile = outfile2, MoreArgs = list(output_Raster = FALSE, overwrite=FALSE, verbose=TRUE))
# file.remove(c(infile, outfile))


## Find min non-NA set values across mask and covariates and sync NAs ####
infile <- list.files(data_processed, pattern = "_masked", full.names = TRUE)
mask_file <- file.path(data_path, "global_mask_ee.tif")
mask_file2 <- file.path(data_processed, "global_mask_ee_nona.tif")
file.copy(mask_file, mask_file2)

for(j in 1:length(infile)){
  print(paste0("processing file = ", j, " of ", length(infile) ,":", infile[j], " for updating mask"))
  # gdalmask(infile = mask_file, mask = infile[j], outfile = mask_file, output_Raster = FALSE, overwrite=TRUE, verbose=TRUE)
  gdalcalc(calc="((A>-9999) & (B==1))", infile = list(A=infile[j], B=mask_file2),outfile = mask_file2,
           NoDataValue=-9999, overwrite=TRUE)
}

for(j in 1:length(infile)){
  print(paste0("processing file = ", j, " of ", length(infile) ,":", infile[j], " using updated mask"))
  gdalmask(infile = infile[j], mask = mask_file2, outfile = infile[j], output_Raster = FALSE, overwrite=TRUE, verbose=TRUE)
}

## Check
summary(raster(mask_file))
summary(raster(mask_file2))
## jgarber comment: Im not sure what the summary file is doing, doesn't seem to count the NA's right, but if you subtract the two rasters and plot them, you will see that some of the cells in the continent have a difference of 1, meaning  that they were masked, but the no data values in the data files are now 0's in the new mask file, this should still mask them out when you run gdalmask, as it masks out every cell where the mask does not equal 1. So i think we are good to go!(see the photo below. White is 0, black is a nodatavalue)
freq(raster(mask_file))
freq(raster(mask_file2))
## jgarber comment: Ran the frequency and as you can see, it has removed some ones from the original mask and made them zeros, this should still work in the gdalmask function. From looking at the map, it seems like it is mostly picking up on large lakes? Looks like the great lakes, Lake victoria, and Lake Baikal are some of the removed cells.
gdalinfo(mask_file, stats = TRUE)
gdalinfo(mask_file2, stats = TRUE)

lapply(infile, gdalinfo, stats = TRUE)


## Create global SpatialPointsDataFrame object ####
mask_path <- file.path(data_processed, "global_mask_ee_nona.tif")
global_mask0 <- raster(mask_path)
global_mask0[which(is.na(global_mask0[]))] <- 0
rpts <- rasterToPoints(global_mask0, spatial=TRUE)
global_mask_df <- data.frame(rpts@data, X=coordinates(rpts)[,1], Y=coordinates(rpts)[,2])     
global_mask_df <- global_mask_df[,-1]
saveRDS(global_mask_df, file = file.path(data_path, "global_sp_pts_df.rds"))


## Separate out srtm variables & save as tifs ####
infile <- list.files(data_processed, pattern = "_masked", full.names = TRUE)
elevation <- raster(infile[grep("srtm", infile)])
writeRaster(elevation, filename = file.path(data_processed, "elevation_masked.tif"))
aspect <- terrain(elevation, opt = "aspect")
writeRaster(aspect, filename = file.path(data_processed, "aspect_masked.tif"))
slope <- terrain(elevation, opt = 'slope')
writeRaster(slope, filename = file.path(data_processed, "slope_masked.tif"))
roughness <- terrain(elevation, opt = "roughness")
writeRaster(roughness, filename = file.path(data_processed, "roughness_masked.tif"))


## Check for correlations, subset and save covariates ####
infile <- list.files(data_processed, pattern = "_masked.tif$", full.names = TRUE)
infile <- infile[-grep("srtm", infile)]

## Issue: some rasters show values and others don’t, but all have values when plotted
## jgarber comment: the only issue is that the aligned datafiles are integers, but the masked data files come out as float datatypes, that means the computation time is a lot longer, and I suspect that the raster package is just not bothering with calculating the statistics. I ran gdal info on the rasters and found that they had reasonable statistics. So I think its just a matter of the raster package being lazy to spend time, and our data are good to go I think.
# for (i in 1:length(infile)) {print(raster(infile[i]))}
# lapply(infile, gdalUtils::gdalinfo, stats = TRUE)

## Covariates for model versus prediction: naming convention as 'model' or 'predict'
covs_model <- infile[-grep('bc45|bc60|bc85', infile)]
# covs_predict_bc45bi70 <- infile[-grep('bio_current|bc60|bc85', infile)]
# covs_predict_bc60bi70 <- infile[-grep('bio_current|bc45|bc85', infile)]
covs_predict_bc85bi70 <- infile[-grep('bio_current|bc45|bc60', infile)]

## Create data frame for model building and prediction
covs_list <- list(covs_model = covs_model, 
                  # covs_predict_bc45bi70 = covs_predict_bc45bi70, 
                  # covs_predict_bc60bi70 = covs_predict_bc60bi70, 
                  covs_predict_bc85bi70 = covs_predict_bc85bi70)

## Save landuse covariate names
lucovs <- tools::file_path_sans_ext(basename(infile))[grep("landuse", basename(infile))]
lucovs <- sub("_masked", "", lucovs)
saveRDS(lucovs, file = file.path(data_processed, "lucovs.rds"))

## Subset data to uncorrelated variables: covs_model
global_mask_df <- readRDS(file.path(data_path, "global_sp_pts_df.rds"))

for(i in 1:length(covs_list)) {
  
  model_data <- stack(covs_list[[i]])
  model_data <- na.omit(cbind(global_mask_df, as.matrix(model_data)))
  
  ## Identify correlated variables
  if(grepl("model", names(covs_list)[[i]])){ 
    names(model_data) <- sub("_masked", "", names(model_data))
    names(model_data) <- sub("_current_", "", names(model_data))
    
    ## Test for correlations
    preds <- colnames(correlations(model_data, thresh = 0.7, N = 200000))
    saveRDS(preds, file = file.path(data_processed, paste0("preds", ".rds")))
    
  } else {
    
    if(grepl("predict", names(covs_list)[[i]])) {
      names(model_data) <- sub("_masked", "", names(model_data))
      names(model_data) <- sub("bc85bi70", "bio", names(model_data))
      
    }
  }
  
  ## Save final datasets
  model_data <- model_data[,colnames(model_data) %in% c(preds, lucovs, "X", "Y")] ## retain landuse, X, Y layers 
  saveRDS(model_data, file = file.path(data_processed, paste0(names(covs_list)[i], ".rds")))
}


mask_file <- file.path(data_path, "global_mask_ee.tif")
globalmask <- raster(mask_file)

