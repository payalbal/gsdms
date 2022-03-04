#' Here we are going to set up the ppm data structure for model fitting.
#' We are going to use the ppmData package and hope that it doesn't break 
#' with very large rasters. We are going to fit these models on a species by 
#' species basis. We might need to generate a single set of quadrature points 
#' (will have a think on this).

ppm_setup <- function(presences, window, covariates){
  require(ppmData)
  quad <- ppmData(presences = presences, window = window, covariates = covariates,
            
          
            )
  
}
  