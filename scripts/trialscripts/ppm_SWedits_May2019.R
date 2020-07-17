SKIP"S EDITS TO ppms.R

#'@title PPM models
#'
#'@description Fit PPMs using GBIF data for IUCN listed species
#'GitHub: https://github.com/payalbal/gsdms/blob/master/R/6_ppms.R
#'
#'
#'

## skip's comments you can load multiple packages using lapply
x <- c('sp', 'raster', 'spatstat.data', 'nlme', 'rpart', 'spatstat', 'ppmlasso', 'parallel', 'doParallel', 'doMC')
lapply(x, require, character.only = TRUE)

## Copy required files to an 'output' folder
if(!dir.exists("./output")) {
  dir.create("./output")
}

## skip's comments: I've put all the data in a a 'data' directory, makes it a little easier. 
# current.folder <- "./data/"
# ppm_files <- paste0(current.folder, c("covariates.rds", "covariates_predict.rds", "bio1.bil", "bio1.hdr", "2019-05-14_gbif_iucnsp.csv"))
## Load covariates
covariates <- readRDS("./data/covariates.rds")
covariates_predict <- readRDS("./data/covariates_predict.rds")

## Load mask
global_mask <- raster("./data/bio1.bil") 
global_mask[which(!is.na(global_mask[]))] <- 1
global_mask[which(is.na(global_mask[]))] <- 0
plot(global_mask)

## Quadrature (background) points & associated covariate data
rpts <- raster::rasterToPoints(global_mask, spatial=TRUE)
backxyz <- data.frame(rpts@data, X=coordinates(rpts)[,1],
                      Y=coordinates(rpts)[,2])     
backxyz <- backxyz[,-1]
backxyz <- na.omit(cbind(backxyz, as.matrix(covariates)))

## skip's comments: its also bad practice to keep on redefining the same thing over and over
## skip's comments: make a new object
backxyz200k <- backxyz[sample(nrow(backxyz), 200000), ] ## skip's comments: I've made it 200000
backxyz200k$Pres <- rep(0, dim(backxyz200k)[1])

# ## Checks
summary(covariates)
summary(backxyz200k)
plot(global_mask, legend = FALSE)
plot(rasterFromXYZ(backxyz200k[,1:3]), col = "black", add = TRUE, legend=FALSE)


## Prediction points 
## skip's comments: Just use all the data.
predxyz <- data.frame(rpts@data, X=coordinates(rpts)[,1],
                      Y=coordinates(rpts)[,2])     
predxyz <- predxyz[,-1]
predxyz <- na.omit(cbind(predxyz, as.matrix(covariates_predict)))
# predxyz <- predxyz[sample(nrow(predxyz), 500000), ]

## Load filtered biodiversity data
gbif <- read.csv("./data/2019-05-14_gbif_iucnsp.csv", header = TRUE)
gbif_sp <- as.vector(unique(gbif$species)) 
spdat <- gbif


## Define global ppm variables
# the pred error you where getting bas because there was different names between the fitted data and the predict data you had:
## "ppmxyz$X" not "X" so the model was just returning the fitted values from the ppmlasso model.
ppm_terms <- names(backxyz)[1:length(names(backxyz))-1]
RCPs <- c(26, 85)

## Specify covariates with interactions
interaction_terms <- c("X", "Y", "srtm")

# Fix n.fits = 20 and max.lambda = 100
lambdaseq <- round(sort(exp(seq(-10, log(100 + 1e-05), length.out = 20)), 
                        decreasing = TRUE),5)



## Initiate log files
## skip's comments: 
## calculate the area of the globe and then work out the weights based on the total area divided by number of points
ar <- raster::area(global_mask)
new_mask <- raster("./data/bio1.bil") 
ar <- mask(ar,new_mask)
totarea <- cellStats(ar,'sum')*1000 ### skip's comments in meters^2
area_offset <- extract(ar, backxyz200k[,c('X','Y')], small = TRUE, fun = mean, na.rm = TRUE)*1000 ## in meters
bkgrd_wts <- c(totarea/area_offset)


## this was meant to calculate the species weights - but it's not quite right
# species_names <- levels(factor(spdat$species))
# spwts <- list()
# for(i in seq_along(species_names)){
#       print(i)
#       spxy <- spdat[spdat$species %in% species_names[i], c(4,3)]
#       names(spxy) <- c("X", "Y")
#       cellNo <- cellFromXY(ar,spxy)
#       cellNoCounts <- table(cellNo)
#       tmp_cell_area <- extract(ar, spxy, fun = mean, na.rm = TRUE)*1000
#       tmp_dat <- data.frame(area=tmp_cell_area,cell_number=cellNo)
#       tmp_wts <- ((tmp_dat$area*1000)/cellNoCounts[as.character(cellNo)])/totarea
#       spwts[[i]] <- tmp_wts
# }

## There is no hard and fast rule about how many covariates to fit to data, and it will change depending on the data and the amount of information in each observation and how they vary with the covariates, how the covariates are correlated ect... But to start I'd only fit sqrt(n_observations) covariates. So if you have 20 occurrences that's 4-5 covariates (including the polynomials) so that's only two variables! You might have to identify the most important variable and start from there.

## let's do a PCA on the data to work out which are the most variable coefs.
Xoriginal=t(as.matrix(backxyz200k))

# Center the data so that the mean of each row is 0
rowmns <- rowMeans(Xoriginal)
X <-  Xoriginal - matrix(rep(rowmns, dim(Xoriginal)[2]), nrow=dim(Xoriginal)[1])

# Calculate P
A <- X %*% t(X)
E <- eigen(A,TRUE)
P <- t(E$vectors)

dimnames(P) <- list(colnames(backxyz200k),paste0("PCA",1:ncol(P)))
df <- as.data.frame(t(P[,1:5]))
df$row.names<-rownames(df)
library(reshape2)
library(ggplot2)
long.df<-melt(df,id=c("row.names"))
pca_plot <- ggplot(long.df,aes(x=row.names,y=variable,fill=value))+
  geom_raster()+
  scale_fill_viridis_c()+
  theme_minimal()

## not surprising that elevation, slope and aspect are all correlated (choose one).
pca_plot


fit_ppms_apply <- function(i, spdat, bkdat, bkwts, interaction_terms, ppm_terms, n.fits=50){
   
  
  species_names <- levels(factor(spdat$species))
  cat('Fitting a ppm to', species_names[i],'\nThis is the', i,'^th model of',length(species_names),'\n')
  logfile <- paste0("./output/ppm_log_",gsub(" ","_",species_names[i]),"_",gsub("-", "", Sys.Date()), ".txt")
  writeLines(c(""), logfile)
  spxy <- spdat[spdat$species %in% species_names[i], c(4,3)]
  names(spxy) <- c("X", "Y")
  spxyz <- extract(covariates, spxy, buffer = 1000000,
                   small = TRUE, 
                   fun = mean, na.rm = TRUE)
  
  ## see notes below... 
  spxyz <- cbind(spxy, spxyz)
  ## Remove remaining NA values & add check for < 20 records
  spxyz <- na.omit(spxyz)
  if (!(dim(spxyz)[1] < 20)) { 
    spxyz$Pres <- rep(1, dim(spxyz)[1])
    
    ## Specify weights for PPM data
    ppmxyz <- rbind(spxyz, bkdat)
    ppmxyz$wt <- NA
    ppmxyz[ppmxyz$Pres == 0,]$wt <- bkwts
    ppmxyz[ppmxyz$Pres == 1,]$wt <- 1e-6 #spwts[[i]]#
    
    ## work out how many covariates
    maxk <- 5 + length(ppm_terms[!(ppm_terms %in% interaction_terms)])*2
    nk <- ceiling(sqrt(nrow(spxy)))
    
    ## Specify ppm formula - independent with quadratic terms
    ## if low number of observations just use the covariates interaction terms for the spatial variables
    if(nk <= 5){   ppmform <- formula(paste0(" ~ poly(", paste0(interaction_terms, collapse = ", "),
                                             ", degree = 2, raw = FALSE)",collapse =""))
    } else  {
      extra_covar <- ceiling((nk - 5)/2)
      if(extra_covar>8) extra_covar <- 8
      ppmform <- formula(paste0(paste0(" ~ poly(", paste0(interaction_terms, collapse = ", "),
                                       ", degree = 2, raw = FALSE)"), 
                                paste0(" + poly(", ppm_terms[!(ppm_terms %in% interaction_terms)][1:extra_covar],
                                       ", degree = 2, raw = FALSE)", 
                                       collapse = ""), collapse =""))
    }
    ## Fit ppm & save output
    cat(paste("\nFitting ppm model for",i , " @ ", Sys.time(), "\n"), 
        file = logfile, append = T)
    cat(paste("   # original records for",i , " = ", dim(spxy)[1], "\n"), 
        file = logfile, append = T)
    cat(paste("   # extracted records for",i , " = ", dim(spxyz)[1], "\n"), 
        file = logfile, append = T)
    ## data was being standardised twice!
    mod <- try(ppmlasso(formula = ppmform, data = ppmxyz, n.fits = n.fits, criterion = "bic",standardise = FALSE), silent=TRUE)
    gc()
    return(mod)
  } else {
    return(NULL)
  }
  
}

spdat <- gbif
bkdat <- backxyz200k
bkwts <- bkgrd_wts
spp <- levels(factor(spdat$species))
mc.cores <- 1
seq_along(spp)
model_list <- parallel::mclapply(1:5, fit_ppms_apply, spdat, #spwts,
                           bkdat, bkwts, interaction_terms, ppm_terms,
                           n.fits=100, mc.cores = mc.cores)

predict_ppms_apply <- function(i, models_list, newdata, bkdat, RCPs = c(26, 85)){
  
  cat('Predicting ', i,'^th model\n')
  
  if(class(models_list)== "try-error") return(NULL)
  
  predmu <- list()
  for (rcp in RCPs) {
    if (rcp == 26) {
      preddat <- newdata[which(names(newdata) %in% 
                                 names(newdata)[-grep('85', names(newdata))])]
      names(preddat) <- names(bkdat)[1:(length(names(bkdat))-2)]
      predmu$rcp26 <- predict.ppmlasso(model_list[[i]], newdata = preddat)
      # predmu <- 1-exp(-predmu) ## gives reative probabilities
      
    } else {
      preddat <- newdata[which(names(newdata) %in% 
                                 names(newdata)[-grep('26', names(newdata))])]
      names(preddat) <- names(bkdat)[1:(length(names(bkdat))-2)]
      predmu$rcp85 <- predict.ppmlasso(model_list[[i]], newdata = preddat)
      # predmu <- 1-exp(-predmu) ## gives reative probabilities
    }
    rm(preddat)
  }
  return(predmu)
}

newdata <- predxyz
prediction_list <- parallel::mclapply(seq_along(model_list), predict_ppms_apply,
                                      model_list, newdata, bkdat, RCPs, mc.cores = mc.cores)

# Prediction for sp1
preds_sp1 <- rasterFromXYZ(cbind(predxyz[,1:2],prediction_list[[1]]))
plot(preds_sp1[[1]],main=spp[1])
points(spdat[spdat$species==spp[1],c(4,3)],pch=16,cex=.5)

# Prediction for sp2
preds_sp2 <- rasterFromXYZ(cbind(predxyz[,1:2],prediction_list[[2]]))
plot(preds_sp2[[1]],main=spp[2])
points(spdat[spdat$species==spp[2],c(4,3)],pch=16,cex=.5)

# Prediction for sp3
preds_sp3 <- rasterFromXYZ(cbind(predxyz[,1:2],prediction_list[[3]]))
plot(preds_sp3[[1]],main=spp[3])
points(spdat[spdat$species==spp[3],c(4,3)],pch=16,cex=.5)

# Ect ect ect