## Update mask based on NAs in covariate stack

align.maskNA <- function(raster_stack, region_mask) {
  print(paste0("# NAs in input mask: ", summary(region_mask)[6]))
  
  for (i in names(raster_stack)){
    if (!sum(is.na(region_mask@data@values)) == sum(is.na(raster_stack[[i]]@data@values))) {
      nona <- which(!is.na(values(region_mask))) # non-na values in mas
      nas <- which(is.na(values(raster_stack[[i]])[nona])) # na values in covariate
      xys <- xyFromCell(region_mask, nona)
      xys <- xys[nas,]
      values(region_mask)[cellFromXY(region_mask, xys)] <- NA
    }
  }
  
  new_mask <- region_mask
  
  print(paste0("# NAs in output mask: ", summary(new_mask)[6]))
  return(new_mask)
}
