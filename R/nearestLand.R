## Author: Nick Goulding
## Notes: Not suitable for 


nearestLand <- function (points, raster, max_distance) {
  # get nearest non_na cells (within a maximum distance) to a set of points
  # points can be anything extract accepts as the y argument
  # max_distance is in the map units if raster is projected
  # or metres otherwise
  
  # function to find nearest of a set of neighbours or return NA
  nearest <- function (lis, raster) {
    neighbours <- matrix(lis[[1]], ncol = 2)
    point <- lis[[2]]
    # neighbours is a two column matrix giving cell numbers and values
    land <- !is.na(neighbours[, 2])
    if (!any(land)) {
      # if there is no land, give up and return NA
      return (c(NA, NA))
    } else{
      # otherwise get the land cell coordinates
      coords <- xyFromCell(raster, neighbours[land, 1])
      
      if (nrow(coords) == 1) {
        # if there's only one, return it
        return (coords[1, ])
      }
      
      # otherwise calculate distances
      dists <- sqrt((coords[, 1] - point[1]) ^ 2 +
                      (coords[, 2] - point[2]) ^ 2)
      
      # and return the coordinates of the closest
      return (coords[which.min(dists), ])
    }
  }
  
  # extract cell values within max_distance of the points
  neighbour_list <- extract(raster, points,
                            buffer = max_distance,
                            cellnumbers = TRUE)
  
  # add the original point in there too
  neighbour_list <- lapply(1:nrow(points),
                           function(i) {
                             list(neighbours = neighbour_list[[i]],
                                  point = as.numeric(points[i, ]))
                           })
  
  return (t(sapply(neighbour_list, nearest, raster)))
}
  
  
  
  ## OTHER
  ## ----------------------------------------------------------------------
  notMissingIdx <- function(raster) {
    # return an index for the non-missing cells in raster
    which(!is.na(getValues(raster)))
  }
  
  missingIdx <- function(raster) {
    # return an index for the missing cells in raster
    which(is.na(getValues(raster)))
  }
  