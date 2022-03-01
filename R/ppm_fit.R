## function for fitting ppms to GBIF species records
## A function for fitting ppms using the weighted Poisson approximation approach

ppm.fit <- function(x, y, weights, family, method=, control=list()){
  
  if(method=="glm")
    ft <- glm.fit(x = x ,y = y/wts, weights = wts, family = poisson())
  if(method=="gam")
    ft <- 
  if(method=="glm.net")
    ft <- 
  else 
    stop("'method' not known choose from 'glm', 'glm.net' or 'gam'")    
  
  
  
  
}

