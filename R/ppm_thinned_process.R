#' @name thinned_process
#' @description This is how we setup the thinned process for the point process 
#' models. Basically we are going to assume a "plug-in" thinned process (we don't 
#' estimate this), using the observed effort for all species (or taxa). This 
#' will give us the probability that any given cell has been sampled (within the 
#' region of interest. A key point is that the thinned process needs to be generated 
#' at the same resolution of the point process models. As the probability is given
#' at the resolution (area of the cells). The output of this function will be a
#' raster which matches the input conditions/data. 
#' We can then use this layer as an offset in the model. It has a few cavets, such as 
#' it assumes that sampling is representative of the effort across the region of interest.
#' It also has a slight computation problem which is when a cell is not sampled 
#' at all, where the log(0) results in numerical issues. We pass a small value 
#' epsilon that essentially makes the probability of cells having not been 
#' sampled very close to zero. In reality this could have implication on the 
#' model and care needs to be taken if weird results start to emerge.
#' @param locs The locations of all points in the region, id's are not needed (ala species).
#' @param window
#' @param res
#' @param extent
#' @param epsilon 
#' 
#' 
library(ppmData)
locs <- snails

thinned_process <- function(locs, window, res=NULL, extent=NULL, epsilon = 1e-9){
  
  
  
  
}

arear <- raster::area(r)
lambda <- rasterize(temp.ala, arear, fun = function(x,...){length(x)})
l1 <- crop(lambda, extent(-180, 0, -90, 90)) ## convert to equal earth: extent(global_mask)
l2 <- crop(lambda, extent(0, 180, -90, 90)) ## convert to equal earth: extent(global_mask)
extent(l1) <- c(180, 360, -90, 90) ## convert to equal earth: extent(global_mask)
lm <- merge(l1, l2)
lm[is.na(lm)]<-0
lm1 <- lm
lambda_mask <- mask(lm,global_mask)
lambda_offset <- 1-exp(-1 * lambda_mask)
lambda_offset[lambda_offset == 0] <- 1e-6
area_rast_mask <- mask(raster::area(lm), global_mask)
## If offset as just area: spp ~ blah + foo + offset(log(area_rast))
## IF offset as area + bias: spp ~ blah + foo + offset(log(off.area_rast)-log(off.bias))

writeRaster(lambda_offset,
            filename = file.path(output_dir, "effort_offset_lambda_0360.tif"),
            format = "GTiff",overwrite = TRUE)
saveRDS(lambda_offset, file.path(output_dir, "effort_offset_lambda_0360.rds"))

writeRaster(area_rast_mask,
            filename = file.path(output_dir, "area_offset_0360.tif"),
            format = "GTiff", overwrite = TRUE)
saveRDS(area_rast_mask,filename = file.path(output_dir, "area_offset_0360.rds"))
