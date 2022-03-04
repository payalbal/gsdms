#' A function for fitting ppms to GBIF species records
#' @description A function for fitting PPMs using the weighted Poisson approximation approach or gibbs (maxent approach)
#' @author Skipton Woolley
#' @param x The design matrix which represents 
#' @param y The presence and absences (backgrounds sites) for the species.
#' @param wts The quadrature weights for the presences and background points.
#' @param method Which model type to 


## 
form_spp <- y/wts ~ 1 + blah + foo
form_bias <- ~ -1 + bias1 + bias2


ppm.fit <- function(x, y, wts, method=c("glm","gam","lasso","ridge","gibbs"),lambda=NULL, control=list()){
  
  # lambda will be a vector of 
  
  if(method=="glm")
    ft <- glm.fit(x = x ,y = y/wts, weights = wts, family = poisson())
  if(method=="gam"){
    ft <- mgcv::gam(formula = ,data = cbind(x,y),weights = wts,family = poisson() )
  }
  if(method=="lasso")
    ft <- glmnet::glmnet(x=y, y=y/wts, weights = wts, family = "poisson", alpha = 1) #lasso
  if(method=="ridge")
    ft <- glmnet::glmnet(x=y, y=y/wts, weights = wts, family = "poisson", alpha = 0) #lasso
  if(method=="gibbs")
    ft <- glm.fit(x = x, y = y/wts, weights = wts, family = binomial())
  else 
    stop("'method' not known choose from 'glm', 'glm.net' or 'gam'")  
  
}

