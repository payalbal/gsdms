## Download global Worldclim data by tiles


## Set working environment ---- ####
rm(list = ls())
gc()
# system("ps")
# system("pkill -f R")


## Load libraries
x <- c('bitops', 'RCurl', 'rstudioapi')
lapply(x, require, character.only = TRUE)
rm(x)

data_dir <- "/home/payalb/gsdms_r_vol/tempdata/research-cifs/uom_data/gsdms_data"
setwd(data_dir)


## WorldClim data - download tiles and stitch to make global layer ---- ####
## Source: https://www.worldclim.org/version1
## author: Chris Ware

## set path to download each zip archive to (will be deleted once used)
## must have .zip ext
zipdst = file.path(data_dir, 'bio_30s.zip')

## set dir to extract raster files to
rasterdst = file.path(data_dir, 'bio_30s/')
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
  # require(bitops)
  # require(RCurl)
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


## WorldClim data - Future GCMs
## Year: 2070, GCM: BCC-CSM1-1
url_bc85bi70 = "http://biogeo.ucdavis.edu/data/climate/cmip5/30s/bc85bi70.zip"
url_bc26bi70 = "http://biogeo.ucdavis.edu/data/climate/cmip5/30s/bc26bi70.zip"
download_from_url(urls = c(url_bc26bi70, url_bc85bi70), zipdst, rasterdst)

  # ## From regSSP scripts
  # ## Global variables
  # regions <- "global"
  # rcps <- c("45", "60", "85")
  # models <- c("BC", "CC", "GS", "HD", "HE", "IP", "MI", "MR", "MC", "MG", "NO")
  #
  # global_mask <- raster(bio_current[1])
  # global_mask[which(!is.na(global_mask[]))] <- 1
  # crs(global_mask) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
  # writeRaster(global_mask, filename = paste0(data_dir, "/mask_global.tif"))
  #
  #
  # ## Downlaoad GCM model predictions
  # urls_bioclim <- c()
  # for (model_i in models){
  #   for (rcp_j in rcps){
  #     fn = paste0(tolower(model_i), rcp_j, 'bi70', '.zip')
  #     thisurl = paste0('https://biogeo.ucdavis.edu/data/climate/cmip5/30s/', fn)
  #     urls_bioclim = c(urls_bioclim, thisurl)
  #   }
  # }
  # library(bitops)
  # library(RCurl)
  # test = lapply(urls_bioclim, url.exists)
  # all(unlist(test))
  #
  # zipdst = file.path(gsdms_data, 'gcm_30s.zip')
  # rasterdst = file.path(gsdms_data, 'gcm_30s/')
  # if(!dir.exists(rasterdst)) {
  #   dir.create(rasterdst)
  # } # end if !dir.exists
  # for (i in 1:length(urls_bioclim)){
  #   download_from_url(urls = urls_bioclim[i], zipdst, rasterdst)
  # }
  #
  # ## Mask data for study regions & create stacks by region and gcms
  # files_gcm <- list.files(file.path(gsdms_data, 'gcm_30s'), full.names = T, recursive = T)
  # gcm_reg_path <- file.path(data_dir, 'gcm_reg')
  # dir.create(gcm_reg_path)
  #
  # for(region in regions){
  #   for(model_i in models){
  #     reg_stack <- list()
  #     file_mod <- files_gcm[grepl(model_i, files_gcm)]
  #     for(j in 1:length(rcps)){
  #       file_mod_rcp <- file_mod[grepl(rcps[j], file_mod)]
  #       temp_stack <- list()
  #       for(f in 1:length(file_mod_rcp)){
  #         reg_mask <- readRDS(file.path(data_dir, paste0("mask_", region, ".rds")))
  #         temp_stack[[f]] <- mask(crop(raster(file_mod_rcp[f]), reg_mask), reg_mask)
  #       }
  #       reg_stack[[j]] <- readAll(brick(temp_stack))
  #     }
  #     saveRDS(reg_stack, file = paste0(gcm_reg_path, "/", region, "_", model_i, ".rds"))
  #   }
  # }
  # rm(temp_stack, reg_stack)
  #
  # ## Extract cell-wise quartiles across GCM
  # quartiles <- c("q1", "q2", "q3")
  #
  # for(region in regions){
  #   gcm <- list.files(gcm_reg_path, full.names = T, pattern = region)
  #   reg_mask <- readRDS(file.path(data_dir, paste0("mask_", region, ".rds")))
  #   r <- reg_mask
  #   inds <- which(!is.na(r[]))
  #
  #   j<-1
  #   for(j in 1:length(rcps)){
  #     saveRDS(stack(), file = paste0(data_dir, "/bio", "q1_", rcps[j], "_", region,  ".rds"))
  #     saveRDS(stack(), file = paste0(data_dir, "/bio", "q2_", rcps[j], "_", region,  ".rds"))
  #     saveRDS(stack(), file = paste0(data_dir, "/bio", "q3_", rcps[j], "_", region,  ".rds"))
  #     print(paste0("processing rcp", rcps[j]))
  #     for(k in 1:19){
  #       print(paste0("processing bioclim var: ", k))
  #       bio <- stack()
  #       for(i in 1:length(models)){
  #         print(paste0("processing model: ", i))
  #         dat <- readRDS(gcm[[i]])[[j]]
  #         bio <- stack(bio, dat[[k]])
  #       }
  #
  #       print(paste0("getting quartiles..."))
  #       df1 <- na.omit(as.matrix(getValues(bio)))
  #       c <-rowQuartiles(df1, probs = c(0.25, 0.5, 0.75))
  #       for(m in 1:3){
  #         bioclim <- readRDS(file = paste0(data_dir, "/bio", quartiles[m], "_", rcps[j], "_", region,  ".rds"))
  #         r[inds] <- c[,m]
  #         names(r) <- paste0("bio", k)
  #         saveRDS(readAll(stack(bioclim, r)), file = paste0(data_dir, "/bio", quartiles[m], "_", rcps[j], "_", region,  ".rds"))
  #       }
  #     }
  #   }
  # }
  # unlink(gcm_reg_path, recursive=T)

## by quartile data in :
## files  <- list.files(file.path(rdata_path), pattern = region, full.names = TRUE)
## bioq <- files[grepl("bioq", files)]