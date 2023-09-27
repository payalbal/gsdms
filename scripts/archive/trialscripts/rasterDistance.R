## Calculate distance to feature and create a raster based on mask
## Author: Roozbeh Valavi (Feb 2019)
##  

rasterDistance <- function(feature, rastermask){
  require(raster)
  require(sf)
  require(progress)
  p <- st_as_sf(rasterToPoints(rastermask, spatial = TRUE))
  p$indx <- st_nearest_feature(st_geometry(p), st_geometry(feature))
  pb <- progress::progress_bar$new(format = " Progress [:bar] :percent in :elapsed",
                                   total=nrow(p), clear=FALSE, width=75) # add progress bar
  for(i in 1:nrow(p)){
    p$dist[i] <- st_distance(st_geometry(p[i,]), st_geometry(feature[p$indx[i],]))
    pb$tick() # update progress bar
  }
  output <- raster::rasterize(p, rastermask, field = "dist")
  return(output)
}


## Trials using NASA road dataset: 
#
library(data.table)
library(sp)
library(raster)
library(rgdal)
library(rgeos)
library(sf)
library(dplyr)
library(progress)
# 
## load mask
# global_mask <- raster("data/processed/global_mask100.tif")
#
## test on a subset of the data
# plot(global_mask) # to draw a ploygon by hand
# pol <- drawPoly()
# crs(pol) <- "+proj=longlat +datum=WGS84 +no_defs" # specify projection
# 
# sub_mask <- crop(global_mask, pol) %>%
#   mask(pol)
# plot(sub_mask)
# writeRaster(sub_mask, "sub_mask", format = "GTiff", options="COMPRESS=LZW", overwrite = TRUE)
# 
## convert the raster mask to SpatialPointsDataFrame
# sub_points <- st_as_sf(rasterToPoints(sub_mask, spatial = TRUE))
# 
## load road data
## data: Global Roads Open Access Data Set, Version 1 (gROADSv1)
## source: http://sedac.ciesin.columbia.edu/data/set/groads-global-roads-open-access-v1
# 
# roads <- sf::st_read(dsn = "data/gROADS_v1.gdb", layer = "Global_Roads")
# 
## 
# sub_road <- st_intersection(roads, st_as_sf(pol))
# plot(st_geometry(sub_road), add = TRUE)
# 
## subset the road data by bla bla, e.g.
# sub_road <- roads %>% 
#   filter(FCLASS == (0:5) & (EXS == c(0,1)))
# plot(st_geometry(sub_road))
# st_write(sub_road, "sub_road.shp")
# 
## for each feature (geometry) in x the index of the nearest feature (geometry) in y 
# sub_points$indx <-  st_nearest_feature(st_geometry(sub_points), st_geometry(sub_road))
# 
## find distnce of raster points to nearest feature
# for(i in 1:nrow(sub_points)){
#   sub_points$dist[i] <- st_distance(st_geometry(sub_points[i,]), st_geometry(sub_road[sub_points$indx[i],]))
#   print(i)
# }
#
## alternately (but discard for now...)
## mapply(function(x, y, z){st_distance(x[z,], y[z,])}, sub_points, sub_road, sub_points$indx)
# 
# plot(raster::rasterize(sub_points, sub_mask, field = "dist"))
# plot(st_geometry(sub_road), add = TRUE)
# 
## Calculate system time and compare methods
# startime <- Sys.time()
# a <- rasterDistance(sub_road, sub_mask)
# Sys.time() - startime
# plot(a)
# 
# startime <- Sys.time()
# aa <- st_distance(sub_points, sub_road)
# sub_points$dist <- apply(aa, 1, min)
# rasterize(sub_points, sub_mask, field = "dist")
# Sys.time() - startime
# 
# 
## Notes. https://gis.stackexchange.com/questions/233443/finding-distance-between-raster-pixels-and-line-features-in-r
## for sp object: dd = gDistance(as(roads, "Spatial"), as(global_mask,"SpatialPoints"), byid=TRUE)
## for sf object: dist.roads <- st_distance(roads, as(global_mask,"SpatialPoints"), by_element = TRUE)
## converting sf to sp: st_as_sf(sp_object)
## converting sp to sf: as (sf_object, "Spatial)
##
## All 'st_functionname' belong to the sf package
##
## Plotting sf ojects: plot(st_geometry(sf_object)) OR plot(sf_object$share/geometry/geom)





