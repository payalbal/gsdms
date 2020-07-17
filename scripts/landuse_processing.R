## From https://github.com/kapitzas/ch2_landusemodel/blob/master/R/1_preprocessing.R
## Simon runs this for all time steps: 1990, 2000, 2006, 2012 and 2018 for Europe
## I only need it for one baseline year, e.g. 2018 for Globe

## Europe 2018 (for trial): https://land.copernicus.eu/pan-european/corine-land-cover
## Global 2015: https://zenodo.org/communities/copernicus-land-cover/search?page=1&size=20 > ProbaV_LC100_epoch2015_global_v2.0.2_discrete-classification_EPSG-4326.tif


## Set work environment ####
rm(list = ls())
library('devtools')
require(gdaltools)
library(rgdal)
library(sp)
library(raster)
require(SpaDES)

data_path <- file.path(getwd(), "data", "output")
temp_path <- file.path(getwd(), "data", "temp")
raw_path <-  file.path(getwd(), "data")

  # eu_data_path <- "/Volumes/discovery_data/gsdms_data/copernicus/pan_european"
  # global_data_path <- "/Volumes/discovery_data/gsdms_data/copernicus/global"


  # ## Download CORINE data ####
  # system("wget https://land.copernicus.eu/land-files/7ac95361f9ac3cecdf37785bc183ff02dd765a16.zip -O data/raw/corine_clc2018.zip")
  # unzip(file.path(raw_path, "corine_clc2018.zip"), exdir = file.path(raw_path))
  # unzip(file.path(raw_path, "clc2018_clc2018_v2018_20_raster100m.zip"), list = TRUE)
  # unzip(
  #   file.path(raw_path, "clc2018_clc2018_v2018_20_raster100m.zip"),
  #   exdir = file.path(raw_path)
  # )
  # file.remove(file.path(raw_path, "corine_clc2018.zip"))


legend <- read.csv(file.path(raw_path, "clc2018_clc2018_v2018_20_raster100m", "CLC2018_CLC2018_V2018_20.txt"), header = FALSE)
urb <- legend[1:10, 1]
agr <- legend[11:22, 1]
forest <- legend[c(23:25, 27:29), 1]
grass <- legend[26, 1]
other <- legend[30:38, 1]
water <- legend[c(39:44), 1]


## Prepapre mask ####
#Europe mask
lu_list <- list.files(raw_path, recursive = TRUE, pattern = "tif$", full.names = TRUE)[-1]
lu <- raster(lu_list[2])

infile <- file.path(temp_path, "mask.tif")
writeRaster(lu[[1]], filename = infile, driver = "GTiff", overwrite = TRUE)
outfile <- file.path(temp_path, "mask_repr.tif")
reproj_ras(infile, outfile, crs = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs", res = c(0.00833333, 0.00833333), ext = c(-32, 46, 26, 72), method = "near")
r <- raster(outfile)
r <- trim(r)

r[r %in% water] <- NA
r[r == 999] <- NA
r[!is.na(r[])] <- 1
saveRDS(readAll(r), file = file.path(data_path, "mask_eur.rds"))

#Countries
shp <- readOGR("~/OneDrive - The University of Melbourne/PhD - Large Files/PhD - Raw Data/Global/ne_10m_admin_0_countries/ne_10m_admin_0_countries.shp")

shp <- readOGR("~/OneDrive - The University of Melbourne/PhD - Large Files/PhD - Raw Data/Global/ne_50m_admin_0_sovereignty/ne_50m_admin_0_sovereignty.shp")
mask <- readRDS(file = file.path(".", "data", "mask_eur.rds"))
europe <- mask
europe[!is.na(europe[])] <- NA
shp <- shp[shp@data$CONTINENT == "Europe",]
l <- nrow(shp@data)

for(i in 1:l){
  print(i)
  m <- raster::extract(mask, shp[i,], cellnumbers = T)
  if(!is.null(m[[1]])){
    europe[m[[1]][,1]] <- i
  }
}
europe[which(europe[]%in%c(11,3))] <- NA #get rid of Russia and Ukraine pixels
europe[which(europe[]%in%c(1))] <- 27 #get rid of Russia and Ukraine pixels
europe <- mask(europe, mask)
mask <- mask(mask, europe)
saveRDS(readAll(europe), file.path(".", "data", "mask_eur.rds"))

## Land use data processing ####
mask <- readRDS(file.path(data_path, "mask_eur.rds"))

lu_list <- list.files(raw_path, recursive = TRUE, pattern = "tif$", full.names = TRUE)[-1]
lu <- stack(lu_list)

lu_stack <- list()

logfile <- file.path(data_path, paste0("corine_log.txt"))
writeLines(c(""), logfile)
j <- 1
for(j in 1:nlayers(lu)){
  cat(paste0("Splitting time step ", j, "\n"), file = logfile, append = TRUE)
  lu_split <- splitRaster(lu[[j]], nx = 5, ny = 5, path = temp_path)
  aggr_layers <- list()
  #cl <- makeCluster(length(lu_split))
  #registerDoParallel(cl)
  cat(paste0("Starting reclassification and aggregation \n"), file = logfile, append = TRUE)
  aggr_layers <- list()
  #aggr_layers <- foreach(i = 1:length(lu_split), .packages = c("raster")) %do% {
  
  for(i in 1:length(fls)){
    cat(paste0("Reclassification split ", i,   "\n"), file = logfile, append = TRUE)
    l <- raster(list.files(temp_path, pattern = paste0('tile', i, '.grd$'), full.name = TRUE))
    l[l %in% urb] <- 1
    l[l %in% agr] <- 2
    l[l %in% grass] <- 3
    l[l %in% forest] <- 4
    l[l %in% other] <- 5
    l[l %in% water] <- NA
    l[l == 999] <- NA
    
    if(!all(is.na(l[]))){
      cat(paste0("Layerization split ", i,   "\n"), file = logfile, append = TRUE)
      layers <- stack(layerize(l))
      layers <- mask(layers, l)
      cat(paste0("Aggregation split ", i,   "\n"), file = logfile, append = TRUE)
      aggr <- aggregate(layers, 10, fun = mean)
      cat(paste0("Reprojection split ", i,   "\n"), file = logfile, append = TRUE)
      aggr_layers[[i]] <- stack(projectRaster(aggr, mask, method = "ngb"))
      nulls <- which(sapply(aggr_layers, FUN = function(x) {is.null(x)}))
    }
    file.remove(list.files(dirname(tempdir()), recursive = TRUE, full.names = TRUE))
  }
  #stopCluster(cl)
  if(length(nulls) > 0){
    aggr_layers <- aggr_layers[-nulls]
  }
  lu_merged <- tile_merge(aggr_layers)
  lu_merged <- mask(lu_merged, mask)
  lu_stack[[j]] <- readAll(lu_merged)
}

saveRDS(lu_stack, file = file.path(data_path, "lu_eur.rds"))

test <- calc(lu_stack[[1]], sum)
plot(test)
summary(test[])
lu_stack <- readRDS(file.path(data_path, paste0("lu_eur.rds")))
for(j in 1:length(lu_stack)){
  stack_out <- lu_stack[[j]]
  for(i in 1:nlayers(stack_out)){
    writeRaster(stack_out[[i]], file = file.path(temp_path, paste0("lu",j,"_", i, "_eur.tif")), format = "GTiff")
  }
}