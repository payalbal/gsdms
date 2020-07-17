#'@title PPM models
#'
#'@description Fit PPMs using GBIF data for IUCN listed species
#'GitHub: https://github.com/payalbal/gsdms/blob/master/R/6_ppms.R
#'
#'


## Set working environment
# system("ps")
# system("pkill -f R")

library(pacman)
p_load(sp, raster, spatstat.data, nlme, rpart, spatstat, ppmlasso, 
       parallel, doParallel, doMC)


## Copy required files to an 'output' folder
# if(!dir.exists("./output")) {
#   dir.create("./output")
# }
# current.folder <- "provide path to downloaded files"
# ppm_files <- paste0(current.folder, c("covariates.rds", "covariates_predict.rds", "bio1.bil", "bio1.hdr", "2019-05-14_gbif_iucnsp.csv"))
# file.copy(ppm_files, "./output/")


## Load covariates
covariates <- readRDS("./data/covariates.rds")
covariates_predict <- readRDS("./data/covariates_predict.rds")


## Load mask
global_mask <- raster("./data/bio1.bil") 
  ## global_mask <- raster::raster("/Volumes/discovery_data/data/raw/climate/wc10/bio1.bil")
global_mask[which(!is.na(global_mask[]))] <- 1
global_mask[which(is.na(global_mask[]))] <- 0


## Quadrature (background) points & associated covariate data
# random sample of 100000 non-NA covraiate values from global grid
rpts <- raster::rasterToPoints(global_mask, spatial=TRUE)
backxyz <- data.frame(rpts@data, X=coordinates(rpts)[,1],
                      Y=coordinates(rpts)[,2])     
backxyz <- backxyz[,-1]
backxyz <- na.omit(cbind(backxyz, as.matrix(covariates)))
backxyz <- backxyz[sample(nrow(backxyz), 100000), ]
backxyz$Pres <- rep(0, dim(backxyz)[1])

    # ## Checks
    # summary(covariates)
    # summary(backxyz)
    # plot(global_mask, legend = FALSE)
    # plot(rasterFromXYZ(backxyz[,1:3]), col = "black", add = TRUE)


## Prediction points
# random sample of 100000 non-NA covraiate values from global grid
predxyz <- data.frame(rpts@data, X=coordinates(rpts)[,1],
                      Y=coordinates(rpts)[,2])     
predxyz <- predxyz[,-1]
predxyz <- na.omit(cbind(predxyz, as.matrix(covariates_predict)))
predxyz <- predxyz[sample(nrow(predxyz), 500000), ]


## Load filtered biodiversity data
gbif <- read.table("./data/2019-05-14_gbif_iucnsp.csv", sep = ",", 
                   header = TRUE)
gbif_sp <- as.vector(unique(gbif$species)) 


## Define global ppm variables
ppmvars <- paste0("ppmxyz$", names(backxyz)[1:length(names(backxyz))-1])
rcps <- c(26, 85)


## Specify covariates with interactions
interact_cov <- paste0("ppmxyz$", c("X", "Y", "srtm"))


# Fix n.fits = 20 and max.lambda = 100
lambdaseq <- round(sort(exp(seq(-10, log(100 + 1e-05), length.out = 20)), 
                        decreasing = TRUE),5)



## Initiate log files
logfile <- paste0("./output/ppm_log_", gsub("-", "", Sys.Date()), ".txt")
writeLines(c(""), logfile)
skipfile <- "./output/skipspecies.txt"
writeLines(c(""), skipfile)


## Start foreach loop for each species in data
ppmout <- list()

# registerDoMC(30)
ppmout <- foreach (i = gbif_sp, #subset by: gbif_sp[c(1,17,24,25)]
                   .export = c("covariates", "backxyz", "predxyz", 
                               "rcps", "ppmvars", "lambdaseq", "ppmout"),
                   .packages = c("spatstat.data", "nlme", "rpart", "spatstat", 
                                 "ppmlasso")) %do% {
                                   
   ## Presence points per species & associated covariate data  
   spxy <- gbif[gbif$species == i, c(4,3)]
   names(spxy) <- c("X", "Y")
   spxyz <- extract(covariates, spxy, buffer = 1000000, small = TRUE, 
                    fun = mean, na.rm = TRUE)
   ## buffer method to avoid NA covariate values for occurrence records on the edge of mask
   ## see notes below... 
   spxyz <- cbind(spxy, spxyz)
   ## Remove remaining NA values & add check for < 20 records
   spxyz <- na.omit(spxyz)
   if (!(dim(spxyz)[1] < 20)) {
     spxyz$Pres <- rep(1, dim(spxyz)[1])
     
         # ## Check by plotting
         # plot(global_mask, legend = FALSE)
         # points(spxy, col = "red", cex = 0.5, pch = 20)
     
     
     ## Specify weights for PPM data
     ppmxyz <- rbind(spxyz, backxyz)
     ppmxyz$wt <- NA
     ppmxyz[ppmxyz$Pres == 0,]$wt <- 1
     ppmxyz[ppmxyz$Pres == 1,]$wt <- 0.1
     
     
     ## Specify ppm formula - independent with quadratic terms
     ppmform <- formula(paste0(paste0(" ~ poly(", paste0(interact_cov, collapse = ", "),
                                      ", degree = 2, raw = FALSE)"), 
                               paste0(" + poly(", ppmvars[!(ppmvars %in% interact_cov)],
                                      ", degree = 2, raw = FALSE)", 
                                      collapse = ""), collapse =""))
     
     
     ## Fit ppm & save output
     cat(paste("\nFitting ppm model for",i , " @ ", Sys.time(), "\n"), 
         file = logfile, append = T)
     cat(paste("   # original records for",i , " = ", dim(spxy)[1], "\n"), 
         file = logfile, append = T)
     cat(paste("   # extracted records for",i , " = ", dim(spxyz)[1], "\n"), 
         file = logfile, append = T)
     mod <- ppmlasso(formula = ppmform, data = ppmxyz, lamb = lambdaseq, 
                     criterion = "bic")

     
     ## Predict & save output
     predmu <- list()
     for (rcp in rcps) {
       if (rcp == 26) {
         preddat <- predxyz[which(names(predxyz) %in% 
                                    names(predxyz)[-grep('85', names(predxyz))])]
         names(preddat) <- names(ppmxyz)[1:(length(names(ppmxyz))-2)]
         predmu$rcp26 <- predict.ppmlasso(mod, newdata = preddat)
         # predmu <- 1-exp(-predmu) ## gives reative probabilities

       } else {
         preddat <- predxyz[which(names(predxyz) %in% 
                                    names(predxyz)[-grep('26', names(predxyz))])]
         names(preddat) <- names(ppmxyz)[1:(length(names(ppmxyz))-2)]
         predmu$rcp85 <- predict.ppmlasso(mod, newdata = preddat)
         # predmu <- 1-exp(-predmu) ## gives reative probabilities
       }
       
       rm(preddat)
     }

     out <- list(species = i, lambda = mod$lambda, 
                 beta = mod$beta, predictions = predmu)
     return(out)
     
     rm(mod, predmu)

   } else {
     warning(i , " skipped because it has < 20 occurrence records. \n")
     cat(paste("INSUFFICIENT RECORDS:", i , "skipped because it has < 20 occurrence records", "\n"), 
         file = logfile, append = T)
     cat(paste("   # original records for", i , " = ", dim(spxy)[1], "\n"), 
         file = logfile, append = T)
     cat(paste("   # extracted records for", i , " = ", dim(spxyz)[1], "\n\n"), 
         file = logfile, append = T)
     cat(paste(i , "\n"), 
         file = skipfile, append = T)
     rm(spxy, spxyz)
     out <- list(species = i, model = "NULL. INSUFFICIENT RECORDS")
     return(out)
   }
  }

## Save outputs to disk
saveRDS(ppmout, file = paste0("./output/ppm_deg2.rds"))
cat(paste("\nTOTAL NUMBER OF PPMS = " , length(ppmout) , "\n"), 
    file = logfile, append = T)
temp <- sapply(ppmout, "[", 1)
temp <- unlist(temp, use.names = FALSE)
cat(paste("\nPPMS fitted for:\n " , paste0(temp, collapse = ", "), "\n"), 
    file = logfile, append = T)


## Clear workspace
rm(list=ls())
gc()




## CHECK 1 : Count occurrence points falling off the mask (edge) - Run on boab
# out <-data.frame(matrix(NA, length(gbif_sp), 5))
# colnames(out) <- c("# records", "extract", "extract.by.buffer", "extract.by.bilinear", "NAcovariates")
# for (i in 1:length(gbif_sp)) {
#   spxy <- gbif[gbif$species == gbif_sp[i], c(4,3)]
#   out[i,1] <- dim(spxy)[1]
#   spxyz1 <- extract(covariates, spxy)
#   out[i,2] <- sum(rowSums(is.na(spxyz1))!=0)
#   spxyz2 <- extract(covariates, spxy, buffer = 1000000, small = TRUE, fun = mean, na.rm = TRUE)
#   out[i,3] <- sum(rowSums(is.na(spxyz2))!=0)
#   spxyz3 <- extract(covariates, spxy, method = 'bilinear')
#   out[i,4] <- sum(rowSums(is.na(spxyz3))!=0)
#   out[i,5] <- paste(dimnames(spxyz1)[[2]][colSums(is.na(spxyz1)) != 0], collapse = ", ")
# }
# out$species <- gbif_sp
# out
# saveRDS(out, file="NNAcounts.rds")

## Note: extracting by buffer gets rid of more NAs than by using bilinear method, 
##  because buffer considers larger area around the point

readRDS(".output/NAcounts.rds")

