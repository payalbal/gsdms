# quick and dirty way to generate random points for weighted poisson regression

#' @title random_background_weights
#' @description Quick way to generate psuedo-random points using the terria 
#' package. Weights are taken as the the area(window)/npoints. The default unit 
#' is km^2, but other units can be used such as meters squared "m" or hectars 
#' "ha".
#' @param window A SpatRaster from terra package which will represent the extent
#'  and resolution of the point process model. 
#' @param npoints The number of background points to use. The default is 10000, 
#' but it is generally recommended to use more when fitting a ppm.
#' @param unit The scale of the area weights, default is kilometers squared "km"
#' but meters squared "m" or hectars "ha" can be used.  
#' @examples 
#' library(ppmData)
#' library(terra)
#' path <- system.file("extdata", package = "ppmData")
#' lst <- list.files(path=path,pattern='*.tif',full.names = TRUE)
#' window <- terra::rast(lst[1])
#' res <- random_background_points(window,10000)
#' plot(window)
#' points(res[,1:2],pch=".")


random_background_points <- function(window = NULL,
                                                 npoints = 10000,
                                                 unit = "km"){
  
  if(is.null(window)) stop("This function requires a window (terra raster) to work.")
  if(class(window)[1]!="SpatRaster") stop("'window' needs to be a 'SpatRaster' from the 'terra' package.")
  
  background_sites <- terra::spatSample(x = window,
                           size = npoints,
                           na.rm = TRUE,
                           as.df = TRUE,
                           xy = TRUE)
   
   window_area <- terra::expanse(window, unit = unit)
   
   bck_wts <- rep(window_area/npoints,npoints)
   
   res <- data.frame(background_sites[,-3],weights=bck_wts)
   
   return(res)
  
}
