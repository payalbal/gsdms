## Fitting PPMs to GBIF data

## Set up work environment ####
rm(list = ls())
gc()
Sys.setenv(TZ='Australia/Melbourne')

# system("ps")
# system("pkill -f R")

x <- c('data.table', 'sp', 'raster', 
       'sense', 'tools', 'bitops', 'RCurl', 
       'rgdal', 'gdalUtils', 'usethis', 'rgeos',
       'plyr', 'dismo',
       'spatstat.data', 'nlme', 'rpart', 'spatstat', 'ppmlasso', 
       'parallel', 'doParallel', 'doMC', 'future.apply', 'future')

lapply(x, require, character.only = TRUE)
rm(x)

## Functions
source("/home/payalb/gsdms_r_vol/tempdata/workdir/gsdms/scripts/0_functions.R")

## Server paths
data_dir <- "/home/payalb/gsdms_r_vol/tempdata/research-cifs/uom_data/gsdms_data"
input_dir <- file.path(data_dir, "processed_data_10k")
output_dir <- file.path(data_dir, "ppm_outputs_10k")
if(!dir.exists(output_dir)) {
  dir.create(output_dir)
}


## Specify global parameters here ####
bkpts <- 10000 # number of background points
n.fits <- 20


## Master log file ####
log_dir <- paste0(output_dir, gsub("-", "", format(Sys.time(), "%F_%H%M")))
dir.create(log_dir)

job_start <- Sys.time()
masterlog <- paste0(log_dir, "/ppm_run_bkp_",bkpts, "_fits_", n.fits, ".txt")
writeLines(c(""), masterlog)
cat(paste0(">> Job start = ", job_start, "\n"), file = masterlog, append = TRUE)
cat(paste0("---------------------------------------------"), file = masterlog, append = TRUE)

datestamp <- gsub("-", "", format(Sys.time(), "%F")) ## used later to extract results


## Scenario labels
ssps <- c("ssp1","ssp5") #paste0("ssp", 1:3)
rcps <- c("45", "85") #c("45", "60", "85")
quartiles <- "q2" #c("q2", "q1", "q3")
scens <- sort(apply(expand.grid(quartiles, ssps), 1, paste0, collapse="_"))
scens_rcps <- sort(apply(expand.grid(quartiles, rcps), 1, paste0, collapse="_"))




## Biodiversity data ####
  # ## Section to be moved to gbif_processing script later
  # gbif_wgs <- fread(file.path(data_dir, "gbif/2019-05-14_gbif_iucnsp.csv"))
  # gbif_wgs[, .N, species]
  # gbif_wgs <- as.data.frame(gbif_wgs)
  # 
  # ## Project biodiv data points in equal area projection - TO MODIFY
  # ## Need to swap out SpatialPointsDataFrame() when working with full GBIF dataset
  # wgs_crs  <-  "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
  # eqarea_crs <- "+proj=eqearth +datum=WGS84 +ellps=WGS84 +wktext" ## check conversion; looks finne in QGIS
  # spdf <- SpatialPointsDataFrame(coords = gbif_wgs[c("decimallongitude", "decimallatitude")],
  #                                data =gbif_wgs, proj4string = CRS(wgs_crs))
  # spdf <- spTransform(spdf, CRSobj = CRS(eqarea_crs))
  # 
  # gbif <- as.data.table(spdf)
  # gbif <- gbif[,c(6,7,2,1,5)]
  # names(gbif)[1:2] <- c("X", "Y")
  # head(gbif)
  # fwrite(gbif, file = file.path(output_dir, "2019-05-14_gbif_iucnsp_EA.csv"))

occdat <- fread(file.path(input_dir, "2019-05-14_gbif_iucnsp_EA.csv"))[,.(X, Y, species)]




## Offset/Samplig bias (using all records) - SKIP TO PLUG IN CODE HERE ####
## Section to be moved to data_processing script later
## ISSUE 1: Uses up too much memory for global analysis causing session/MRC instance crash
## ISSUE 2: raster::area() only works with WGS projection
## NOTES: 
## Trial terra package here
## Alternarte approach: https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0079168

  # ## Option 1:
  # ## Create empty raster
  # global_mask <- raster(file.path(input_dir, "globalmask_10k_ee_minNA.tif"))
  # r <- global_mask
  # r[] <- 0
  # 
  # ## Count the number of records per cell
  # tab <- table(cellFromXY(r, occdat[,1:2])) ## sample occ_data when using full dataset
  # r[as.numeric(names(tab))] <- tab
  # r <- mask(r,global_mask)
  # 
  # ## Make zeros a very small number otherwise issues with log(0).
  # r[r[]==0] <- 1e-6
  # arear <- raster::area(r)
  # 
  # ## Calculate the probability that a cell has been sampled while accounting for area differences in lat/lon
  # off.bias <- (-log(1-exp(-r*arear)) - log( arear))
  # names(off.bias) <- "off.bias"

  
  # ## Option 2:
  # temp.ala<-SpatialPoints(occdat[,c(1, 2)], proj4string=crs(global_mask))
  # 
  # ## Craete empty raster to catch data
  # global_mask <- raster(file.path(input_dir, "globalmask_10k_ee_minNA.tif"))
  # r <- global_mask
  # r[] <- 0
  # 
  # ## Calculate area to adjust for
  # ## QUESTION: Difference between lambda_offset and area_offset?
  # arear <- raster::area(r)
  # lambda <- rasterize(temp.ala, arear, fun = function(x,...){length(x)})
  # l1 <- crop(lambda, extent(-180, 0, -90, 90)) ## convert to equal earth: extent(global_mask)
  # l2 <- crop(lambda, extent(0, 180, -90, 90)) ## convert to equal earth: extent(global_mask)
  # extent(l1) <- c(180, 360, -90, 90) ## convert to equal earth: extent(global_mask)
  # lm <- merge(l1, l2)
  # lm[is.na(lm)]<-0
  # lm1 <- lm
  # lambda_mask <- mask(lm,global_mask)
  # lambda_offset <- 1-exp(-1 * lambda_mask)
  # lambda_offset[lambda_offset == 0] <- 1e-6
  # area_rast_mask <- mask(raster::area(lm), global_mask)
  # ## If offset as just area: spp ~ blah + foo + offset(log(area_rast))
  # ## IF offset as area + bias: spp ~ blah + foo + offset(log(off.area_rast)-log(off.bias))
  # 
  # writeRaster(lambda_offset,
  #             filename = file.path(output_dir, "effort_offset_lambda_0360.tif"),
  #             format = "GTiff",overwrite = TRUE)
  # saveRDS(lambda_offset, file.path(output_dir, "effort_offset_lambda_0360.rds"))
  # 
  # writeRaster(area_rast_mask,
  #             filename = file.path(output_dir, "area_offset_0360.tif"),
  #             format = "GTiff", overwrite = TRUE)
  # saveRDS(area_rast_mask,filename = file.path(output_dir, "area_offset_0360.rds"))
  # 
  # ## Add offset to covariate stacks for model and prediction
  # off.bias <- ...



## Covariate data ####
## >> Background points ####
## See data_processing for background point generation & if an alternate approach is to be used
backxyz <- fread(file.path(input_dir, "covariates_model.csv"))

## Add offset
backxyz <- cbind(bacxyz, off.bias)

## Sample points
backxyz200k <- backxyz[sample(nrow(backxyz), 200000), ] ## ~1/5th : to thik about
backxyz200k$Pres <- rep(0, dim(backxyz200k)[1])

## Checks
summary(backxyz200k)
global.mask <- raster(file.path(input_dir, "globalmask_10k_ee_minNA.tif"))
plot(global.mask, legend = FALSE)
# plot(rasterFromXYZ(backxyz[,1:3]), col = "grey", add = TRUE, legend=FALSE)
plot(rasterFromXYZ(backxyz200k[,1:3]), col = "red", add = TRUE, legend=FALSE)
points(occdat[,1:2], col = "yellow", pch = 20, cex = 0.3)


## >> Prediction points (according to scenarios) ####
# predxyz <- fread(file.path(input_dir, "covariates_predict_rcp45.csv"))
predxyz <- fread(file.path(input_dir, "covariates_predict_rcp85.csv"))

## Add offset
predxyz <- cbind(covs_predict, off.bias)



## >> Check that NA syncs across all covariate layers ####
## use gdal_calc function if NAs need to be synced






## Define model parameters ####
## >> Specify covariates with interactions ###
interaction_terms <- c("X", "Y")

## >> Specify factor variable ####
factor_terms <- names(backxyz)[grep("lu", names(backxyz))]


## >> Define continuous ppm variables to be used in the model ####
## NOTE: X, Y & landuse are called in the ppm formula directly, hence removed here
ppm_terms <- names(backxyz)[!grepl("off.bias", names(backxyz))]
ppm_terms <- ppm_terms[!(ppm_terms %in% interaction_terms)]
ppm_terms <- ppm_terms[!(ppm_terms %in% factor_terms)]


# ## Fix max.lambda = 100 - SKIP, IS THIS NEEDED? ####
# lambdaseq <- round(sort(exp(seq(-10, log(max.lambda + 1e-05), length.out = n.fits)), 
#                         decreasing = TRUE),5)


## >> Estimate weights for background points - SKIP TO FIX/MODIFY ####
## See ppmlasso::ppmdat 
## To replace with Dirichlet tessalation
ppmlasso_weights <- function (sp.xy, quad.xy, coord = c("X", "Y")){
  sp.col = c(which(names(sp.xy) == coord[1]), which(names(sp.xy) ==
                                                      coord[2]))
  quad.col = c(which(names(quad.xy) == coord[1]), which(names(quad.xy)
                                                        == coord[2]))
  X.inc = sort(unique(quad.xy[, quad.col[1]]))[2] -
    sort(unique(quad.xy[, quad.col[1]]))[1]
  Y.inc = sort(unique(quad.xy[, quad.col[2]]))[2] -
    sort(unique(quad.xy[, quad.col[2]]))[1]
  quad.0X = min(quad.xy[, quad.col[1]]) - floor(min(quad.xy[,
                                                            quad.col[1]])/X.inc) * X.inc
  quad.0Y = min(quad.xy[, quad.col[2]]) - floor(min(quad.xy[,
                                                            quad.col[2]])/Y.inc) * Y.inc
  X = c(sp.xy[, quad.col[1]], quad.xy[, quad.col[1]])
  Y = c(sp.xy[, quad.col[2]], quad.xy[, quad.col[2]])
  round.X = round((X - quad.0X)/X.inc) * X.inc
  round.Y = round((Y - quad.0Y)/Y.inc) * Y.inc
  round.id = paste(round.X, round.Y)
  round.table = table(round.id)
  wt = X.inc * Y.inc/as.numeric(round.table[match(round.id,
                                                  names(round.table))])
}


## >> Estimate species weights - SKIP TO REVIEW ####
## (delete/modify/keep/add not to work on this later?)
# spdat <- gbif
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


## >> Initialise log file for comparing model runtimes ####
timelog <- paste0(log_dir, "/ppm_timelog_",Sys.Date(), ".txt")
writeLines(c(""), timelog)
cat(paste(c("i", "species_name", "pr_pts", "back_pts",	"ppmfit_0folds_min",	
            "ppmfit_1fold_avg_min",	"ppmfit_5folds_min",	"ppmeval_train_sec",	
            "ppmeval_test_sec", "distribution", "\n"), collapse = ","), 
    file = timelog, append = T)


## >> Define model function ####
fit_ppms_apply <- function(i, spdat, bkdat, interaction_terms, 
                           ppm_terms, species_names, mask_path, 
                           n.fits=50, min.obs = 60, 
                           modeval = TRUE, 
                           output_list = FALSE) {
  
  ## >> >> Initialise log file for species ####
  cat("Fitting a ppm to", species_names[i],"\nThis is the", i,"^th model of",length(species_names),"\n")
  specieslog <- paste0(log_dir, "/ppm_log_",gsub(" ","_",species_names[i]),"_",gsub("-", "", Sys.Date()), ".txt")
  writeLines(c(""), specieslog)
  cat(paste("Fitting a ppm to", species_names[i],"\nThis is the", i,"^th model of",length(species_names),"\n"),
      file = specieslog, append = T)
  cat(paste("-------------------------------------------\n"),
      file = specieslog, append = T)
  ppm_dat_start <- Sys.time()
  
  ## >> >> Get species specific data ####
  spxy <- spdat[spdat$species %in% species_names[i], c(1,2)]
  names(spxy) <- c("X", "Y")
  
  # # ## Check if points fall outside the mask - TO REVIEW
  # ## Decide what to do here: discard or move points
  # ## This check/fix will be moved to the gbif processing script later
  # plot(global_mask)
  # points(spxy)
  # 
  # 
  # ## Move species occurrence points falling off the mask to nearest 'land' cells
  # ## First find points which fall in NA areas on/off the raster
  # global_mask <- raster(mask_path)
  # vals <- extract(global_mask, spxy)
  # outside_mask <- is.na(vals)
  # if(sum(outside_mask) > 0){
  #   outside_pts <- spxy[outside_mask, ]
  #   ## find the nearest land within 5 decimal degrees of these
  #   land <- nearestLand(outside_pts, global_mask, 1000000)
  #   ## replace points falling in NA with new points on nearest land
  #   spxy[outside_mask, ] <- land
  #   ## count how many were moved
  #   sum(!is.na(land[, 1]))
  # }
  # ## Checks
  # # nrow(unique(outside_pts))
  # # nrow(unique(land))
  # ## NOTE: We lose data because number of unique locations is reduced.
  
  
  ## >> >> Extract covariates for species presence points - TO FIX ####
  ## Too long to run raster funcs, alternative gdal function to extract points from grid
  
  ## Extract covariate values excluding land use 
  cov_mod <- rasterFromXYZ(bkdat)
  r <- cov_mod[[-grep("lu", names(cov_mod))]]
  spxy_nolu <- extract(r, spxy,
                       fun = mean, na.rm = TRUE)
  # ## NOTE: Extracting by buffer takes too long to run ...
  # spxyz_nolu <- extract(r, spxy, buffer = 1000000, small = TRUE, 
  #                       fun = mean, na.rm = TRUE)
  
  ## Extract landuse values: Take the non-NA value at shortest distance from point
  r = cov_mod[[grep("lu", names(cov_mod))]] ## landuse raster
  spxy_lu <- extract(r, spxy, method='simple', na.rm=TRUE, factor = TRUE) 
  ## NOTE: extract(r, spxy, method = 'simple') gives NAs
  ## WARNING: Even with the na.rm argumennt, we may still get NAs.
  ## Alternate approach (below) but it takes too long to run; no good, abort!
  # spxy_lu <- apply(X = spxy, MARGIN = 1, FUN = function(X) r@data@values[which.min(replace(distanceFromPoints(r, X), is.na(r), NA))])
  
  spxyz <- cbind(spxy, spxy_lu, spxy_nolu)
  spxyz <- na.omit(spxyz)
  
  ## >> >> Return error if NAs found in data ####
  if (all(is.na(spxyz))) {
    stop("NA values not allowed in model data.")
  }
  
  ## >> >> Return NULL for species if number of observations < 20 (i.e. do not fit model) ####
  nk <- nrow(spxyz)
  
  if (nk < 20) {
    return(NULL)
  } else {
    
    ## >> >> Build model ####
    ## Add 'Pres' column and calculate weights for PPM data
    ppmxyz <- rbind(cbind(spxyz, Pres = 1),
                    cbind(bkdat, Pres = 0))
    wts <- ppmlasso_weights(sp.xy = spxyz, quad.xy = bkdat)
    ppmxyz <- cbind(ppmxyz, wt = wts)
    
    ## Log run time for data preparationn
    ppm_dat_end <- Sys.time() - ppm_dat_start
    cat(paste("\nRun time for data preparation: ", ppm_dat_end, attr(ppm_dat_end, "units"), "\n"), 
        file = specieslog, append = T)
    
    ## !! Specify PPM formula based on number of observatios - TO FIX ####
    ## If number of observations < 20, do not fit a model.
    ## If number of observations >= 20, fit a model accordig to the followig rule:
    ##  1. If number of observations <= min.obs (default as 60), 
    ##    use only the spatial covariates i.e. lat and long and offset. This gives 6 
    ##    terms: X, Y, X^2, Y^2, X:Y (linear, quadratic with interaction), offset
    ##  2. If number of observatons > min.obs, 
    ##    use the one-in-ten rule (https://en.wikipedia.org/wiki/One_in_ten_rule)
    ##    where, an additional covariate is added for every 10 additonal observations.
    ##    Because we fit poly with degree 2, adding a covariate will add two 
    ##    terms x and x^2 to the model. Therefore, to operationalize this rule, 
    ##    we add 1 raw covariate for every 20 additional observations. 
    ## NOTE: if area offset considered, use "+ offset(log(off.area)-log(off.bias))"
    
    ## *** ISSUE: There is inherent bias in how covariates are added i.e.
    ##  which covariates are included in the model as they are selected 
    ##  based on their order in the dataframe.
    
    if(nk <= min.obs){
      
      ppmform <- formula(paste0(" ~ poly(", paste0(interaction_terms, collapse = ", "),
                                ", degree = 2, raw = FALSE) + offset(log(off.bias))",collapse =""))
    } else {
      
      ## Fit covariates as per 1-in-10 rule
      extra_covar <- ceiling((nk - min.obs)/20)
      extra_covar <- ifelse(extra_covar > length(ppm_terms), length(ppm_terms), extra_covar) 
      
      ## !! TO FIX LU VARIABLE - each LU class showing as categorical ####
      ppmform <- formula(paste0(" ~ poly(", paste0(interaction_terms, collapse = ", "),
                                       ", degree = 2, raw = FALSE) + ", 
                                       paste0("factor(", factor_terms, collapse = " + ", ")"),  
                                       paste0(" + poly(", ppm_terms[1:extra_covar],
                                              ", degree = 2, raw = FALSE)",
                                              collapse = ""), " + offset(log(off.bias))" ,
                                       collapse =""))
    }
  }
}

    
    ## Fit ppm & save output
    cat(paste("\nFitting ppm model for species ",i , ": ", spp[i], " @ ", Sys.time(), "\n"), 
        file = specieslog, append = T)
    cat(paste("   # presence points for species (original) = ", 
              dim(spxy)[1], "\n"), file = specieslog, append = T)
    cat(paste("   # presence points for species (extracted) = ", 
              dim(spxyz)[1], "\n\n"), file = specieslog, append = T)
    cat(paste("   # total data points for species (presence + background) = ", 
              dim(ppmxyz)[1], "\n\n"), file = specieslog, append = T)
    
    if(!modeval){ ## if model_eval == FALSE
      ## Fit model
      ppm_mod_start <- Sys.time()
      
      mod <- tryCatch(ppmlasso(formula = ppmform, data = ppmxyz, n.fits = n.fits, 
                               criterion = "bic", standardise = FALSE),
                      error = function(e){ 
                        cat(paste("\nModel ",i," for species ", spp[i], " has errors. \n"),
                            file = specieslog, 
                            append = T)
                        return(NA) 
                      })
      
      if(is.na(mod)){
        cat(paste("\nERROR: Output for model ",i , " for species ", spp[i], " is NA. \n"), 
            file = specieslog, append = T)
      }
      
      # ## Print warnings to screen (for sequence runs only)
      # cat('Warnings for ', species_names[i],':\n')
      # warnings()
      
      # ## Capture messages and errors in file (for sequence runs only)
      # sink(specieslog, type = "message", append = TRUE, split = FALSE)
      # try(warnings())
      # ## reset message sink and close the file connection
      # sink(type="message")
      # close(specieslog)
      
      ## Log run time for model run
      ppm_mod_end <- Sys.time() - ppm_mod_start
      cat(paste(">> Run time for ppm model fitting: ", ppm_mod_end, 
                attr(ppm_mod_end, "units"), "\n"), 
          file = specieslog, append = T)
      
      ## Log time for comparison
      time_vec <- c(dim(spxy)[1], dim(bkdat)[1], ppm_mod_end, rep(NA, 5))
      cat(paste(c(i, trimws(gsub(" ", "_", spp[i])), time_vec, "\n"), collapse = ","),
          file = timelog, append = T)
      
      if (output_list == TRUE) {
        return(mod)
      } else {
        saveRDS(mod, paste0(output_dir, "/ppm_out_",gsub(" ","_",species_names[i]),"_",gsub("-", "", Sys.Date()), ".rds"))
      }
      
    } else { ## if model_eval == TRUE
      
      ppm_mod_start <- Sys.time()
      
      ## Run a basic k-fold for a ppmlasso
      ## ... Brendan might want to check this
      nk <- k <- 5 # number of folds
      kfold_train <- list()
      kfold_test <- list()
      rm(folds)
      
      ## Sample without replacement 
      createFolds <- function(x,k){
        n <- nrow(x)
        x$folds <- rep(1:k,length.out = n)[sample(n,n)]
        x
      }
      
      folds <- plyr::ddply(ppmxyz,.(Pres),createFolds,k = nk)
      for(ii in 1:nk){
        kfold_test[[ii]]<-folds[folds$folds==ii,] ## dim(folds)[1]/k rows in kfolds_test[[ii]]
        kfold_train[[ii]]<-folds[folds$folds!=ii,] ## all rows not in kfolds_test[[ii]]
        # print(colSums(folds[folds$folds==ii,],na.rm = T))
      }
      
      ppmCV <- list()
      ppmCV_dats <- list()
      ppmCV_end_kmodels <- c() ## to log time per run
      
      ## Fit model
      ## fit k models with a random (without replacement) basic K-fold
      for (ii in seq_len(k)){
        ppmCV_dats[[ii]] <- kfold_train[[ii]]
        ppmCV[[ii]] <- tryCatch(ppmlasso(formula = ppmform, 
                                         data = ppmCV_dats[[ii]], 
                                         n.fits = n.fits, 
                                         criterion = "bic", 
                                         standardise = FALSE),
                                error = function(e){ 
                                  cat(paste("\n Model",ii,"/",k," for species ", i , ": ", spp[i], " has errors. \n"),
                                      file = specieslog, 
                                      append = T)
                                  return(NA) 
                                })
        
        ppmCV_end <- Sys.time()
        
        if(is.na(ppmCV[[ii]])){
          cat(paste("\nERROR: Output for model ",ii,"/",k," for species ", i , ": ", spp[i], " is NA. \n"), 
              file = specieslog, append = T)
          stop("\nERROR: Output for model ",ii,"/",k," for species ", i , ": ", spp[i], " is NA \n")
        } else {
          cat(paste("\nModel",ii,"/",k, " done @ ", ppmCV_end, "\n"),
              file = specieslog, append = T)
        }
        ppmCV_end_kmodels[[i]] <- ppmCV_end
      }
      
      ## Log run time for k models run
      ppm_mod_end <- Sys.time() - ppm_mod_start
      cat(paste(">> Run time for all", k," models for species ",
                i , ": ", ppm_mod_end, 
                attr(ppm_mod_end, "units"), "\n"), 
          file = specieslog, append = T)
      
      ## Model evaluation with training data
      ppm_traineval_start <- Sys.time()
      
      cat(paste("\n\nModel evaluation with training data: ", Sys.time(), "\n"), 
          file = specieslog, append = T)
      
      cell_area <- prod(res(global_mask))
      train.preds <- lapply(1:k,function(x)predict(ppmCV[[x]],
                                                   newdata=ppmCV_dats[[x]])*cell_area)
      model.evaluations <- lapply(1:k,function(x)dismo::evaluate(train.preds[[x]][ppmCV_dats[[x]]$Pres==1],
                                                                 train.preds[[x]][ppmCV_dats[[x]]$Pres==0]))
      
      # mean AUC from k-folds the model
      train.meanAUC <- mean(sapply(model.evaluations,function(x)x@auc))
      
      # ## ROC curves
      # par(mfrow=c(3,2))
      # sapply(1:k,function(ii)plot(model.evaluations[[ii]],"ROC"))
      # 
      # ## TPR plots
      # par(mfrow=c(3,2))
      # sapply(1:k,function(ii)plot(model.evaluations[[ii]],"TPR"))
      # 
      # ## Density plots
      # par(mfrow=c(3,2))
      # sapply(1:k,function(ii)density(model.evaluations[[ii]]))
      
      ## Log run time for evaluation using training data
      ppm_traineval_end <- Sys.time() - ppm_traineval_start
      cat(paste(">> Run time for model evaluation with training data: ", 
                ppm_traineval_end, attr(ppm_traineval_end, "units"), "\n"), 
          file = specieslog, append = T)
      
      
      ## Model evaluations with test data
      ppm_testeval_start <- Sys.time()
      
      cat(paste("\n\nModel evaluation with test data: ", Sys.time(), "\n"), 
          file = specieslog, append = T)
      
      test.preds <- lapply(1:k,function(x)predict(ppmCV[[x]],
                                                  newdata=kfold_test[[x]])*cell_area)
      test.evaluations <- lapply(1:k,function(x)dismo::evaluate(test.preds[[x]][kfold_test[[x]]$Pres==1],
                                                                test.preds[[x]][kfold_test[[x]]$Pres==0]))
      
      
      ## Mean AUC from k-folds the model
      test.meanAUC <- mean(sapply(test.evaluations,function(x)x@auc))
      
      # ## ROC curves
      # par(mfrow=c(3,2))
      # sapply(1:k,function(ii)plot(test.evaluations[[ii]],"ROC"))
      # 
      # ## TPR plots
      # par(mfrow=c(3,2))
      # sapply(1:k,function(ii)plot(test.evaluations[[ii]],"TPR"))
      # 
      # ## Density plots
      # par(mfrow=c(3,2))
      # sapply(1:k,function(ii)density(model.evaluations[[ii]]))
      
      ## Log run time for evaluation using test data
      ppm_testeval_end <- Sys.time() - ppm_testeval_start
      cat(paste(">> Run time for model evaluation with test data: ", 
                ppm_testeval_end, attr(ppm_testeval_end, "units"), "\n"), 
          file = specieslog, append = T)
      
      ## Save output
      mod <- list(ppmCV = ppmCV, 
                  train_eval = model.evaluations, 
                  train_AUC = train.meanAUC, 
                  test_eval = test.evaluations,
                  test_AUC = test.meanAUC)
      
      ## Log time for comparison
      time_vec <- c(dim(spxy)[1], dim(bkdat)[1], NA, mean(ppmCV_end_kmodels), 
                    ppm_mod_end, ppm_traineval_end, ppm_testeval_end, NA)
      cat(paste(c(i, trimws(gsub(" ", "_", spp[i])), time_vec, "\n"), collapse = ","),
          file = timelog, append = T)
      
      gc()
      
      if (output_list == TRUE) {
        return(mod)
      } else {
        saveRDS(mod, paste0(output_dir, "/ppm_out_",gsub(" ","_",species_names[i]),"_",gsub("-", "", Sys.Date()), ".rds"))
        out <- list(likelihood = sapply(mod$ppmCV, `[[`, 8), train_AUC = mod$train_AUC, test_AUC = mod$test_AUC)
        return(out)
      }
    }
  }
}






## III. Fit models & log info ####
ppm_start <- Sys.time()
cat(paste0("\n\nModel runs start = ", ppm_start, "\n"), file = masterlog, append = TRUE)

trialN <- 10

spdat <- as.data.frame(occdat)
bkdat <- backxyzK
modeval <- TRUE
spp <- unique(spdat$species)[1:trialN]
mc.cores <- trialN
seq_along(spp)
ppm_models <- list()

plan(multiprocess, workers = mc.cores)
options(future.globals.maxSize = +Inf) ## CAUTION: Set this to a value, e.g. availablecores-1?/RAM-10?
ppm_models <- future.apply::future_lapply()(1:length(spp), fit_ppms_apply, spdat,
                                            bkdat, interaction_terms, ppm_terms,
                                            species_names = spp,
                                            mask_path = file.path(input_dir, "aligned_mask_aus.rds"),
                                            n.fits = n.fits, min.obs = 60,
                                            modeval = modeval, output_list = FALSE)
# ppm_models <- parallel::mclapply(1:length(spp), fit_ppms_apply, spdat,
#                                  bkdat, interaction_terms, ppm_terms,
#                                  species_names = spp,
#                                  mask_path = file.path(input_dir, "aligned_mask_aus.rds"),
#                                  n.fits = n.fits, min.obs = 60, mc.cores = mc.cores,
#                                  modeval = modeval, output_list = FALSE)

# ppm_models <- lapply(1:length(spp), fit_ppms_apply, spdat,
#                      bkdat, interaction_terms, ppm_terms,
#                      species_names = spp,
#                      mask_path = file.path(input_dir, "aligned_mask_aus.rds"),
#                      n.fits=10, min.obs = 60, modeval = TRUE, output_list = TRUE)
# names(ppm_models) <- tolower(gsub(" ", "_", spp))
# saveRDS(ppm_models, paste0(output_dir, "/ppm_out_", gsub("-", "", Sys.Date()), ".rds"))

# ppm_models <- fit_ppms_apply(1, spdat, bkdat, interaction_terms,
#                              ppm_terms, species_names = spp,
#                              mask_path = file.path(input_dir, "aligned_mask_aus.rds"),
#                              n.fits=5, min.obs = 60, modeval = TRUE, output_list = TRUE)
# names(ppm_models) <- tolower(gsub(" ", "_", spp))
# saveRDS(ppm_models, paste0(output_dir, "/ppm_out_", gsub("-", "", Sys.Date()), ".rds"))

ppm_runtime <- Sys.time()-ppm_start
cat(paste0("\n\nfit_ppms_apply() run time for ", length(spp), " species: ", 
           ppm_runtime, " ", attr(ppm_runtime, "units"), "\n\n"))
cat(paste0("Model runs stop = ", Sys.time(), "\n"), file = masterlog, append = TRUE)
cat(paste0(">> Model runtime = ", ppm_runtime, " ",
           attr(ppm_runtime, "units"), "\n"), file = masterlog, append = TRUE)

cat(paste0("\n\n---------------------------------------------\n"), 
    file = masterlog, append = TRUE)
cat(paste0("Run details: \n", 
           "modeval = ", modeval, "\n",
           "# species = ", length(spp), "\n",
           "# background points = ", bkpts, "\n",
           "# n.fits = ", n.fits, "\n"),
    file = masterlog, append = TRUE)
cat(paste0("\n---------------------------------------------\n\n"), 
    file = masterlog, append = TRUE)

cat("Session details: \n", file = masterlog, append = TRUE)
cat(paste0("Date:\n"), file = masterlog, append = TRUE)
sink(file = masterlog, append = TRUE)
Sys.Date()
sink(NULL)
cat("\n\nSession Info:\n", file = masterlog, append = TRUE)
sink(file = masterlog, append = TRUE)
sessionInfo()
sink(NULL)
gc()

## Load ppm_models from rds files (when output_list = FALSE) ####
## Species ppm files
ppm_files <- list.files(output_dir, pattern = "ppm_out", full.names = TRUE)

## Extract species names from files
ppm_sp <- regmatches(ppm_files,regexec(paste0("ppm_out_(.*?)_", datestamp), ppm_files))
ppm_sp <- sapply(ppm_sp, `[[`, 2)

## Create list of results
ppm_out <- list()
for(i in 1:length(ppm_files)){
  ppm_out[[i]] <- readRDS(ppm_files[i])
}
names(ppm_out) <- ppm_sp

## Check model fits/predictions
names(ppm_out[[1]])
ppm_test_auc <- sapply(ppm_out, `[[`, 5)

## Check species distributions
global_mask <- readRDS(file.path(input_dir, "aligned_mask_aus.rds"))
par(mfrow = c(2,5))
for (i in 1:length(spp)){
  ## Define species specific data & remove NA (if any)
  spxy <- spdat[spdat$species %in% spp[i], c(1,2)]
  names(spxy) <- c("X", "Y")
  
  ## Check if points fall outside the mask
  plot(global_mask, axes = FALSE, legend = FALSE, box = FALSE, main = spp[i])
  points(spxy)
}

## Check for NULL models - WHERE ARE MY NULL MODELS!! (not being stored in output - why?)
length(ppm_out[!sapply(ppm_out,is.null)])
length(ppm_out[sapply(ppm_out,is.null)])
which(sapply(ppm_out,is.null))


## IV. Catch errors in models & save outputs ####
error_list <- list()
model_list <- list()
n <- 1
m <- 1
errorfile <- paste0(log_dir, "/errorfile_", gsub("-", "", Sys.Date()), ".txt")

for(i in 1:length(ppm_models)){
  if(!class(ppm_models[[i]])[1] == "try-error") {
    model_list[[n]] <- ppm_models[[i]]
    n <- n+1
  }else{
    print(paste0("Model ",i, " for '", spp[i], "' has errors"))
    cat(paste(i, ",", spp[i], "\n"),
        file = errorfile, append = T)
    error_list[[m]] <- ppm_models[[i]]
    m <- m+1
  }
}

saveRDS(model_list, file = paste0(output_dir, "/modlist_",  gsub("-", "", Sys.Date()), ".rds"))
saveRDS(error_list, file = paste0(output_dir, "/errlist_",  gsub("-", "", Sys.Date()), ".rds"))



## V. Predict and save output ####
## Define prediction function ####
predict_ppms_apply <- function(i, models_list, newdata, bkdat, RCPs = c(45, 85)){
  
  cat("Predicting ", i,"^th model\n")
  
  if(class(models_list)== "try-error") { ## redundant because error models are removed
    return(NULL)
  } else {
    
    predmu <- list()
    
    ## Current predictions
    preddat <- newdata[which(names(newdata) %in% 
                               names(newdata)[-grep("45|85", names(newdata))])]
    predmu$current <- predict.ppmlasso(models_list[[i]], newdata = preddat)
    
    ## Future predictions
    for (rcp in RCPs) {
      if (rcp == 45) {
        preddat <- newdata[which(names(newdata) %in% 
                                   names(newdata)[-grep("85|bio", names(newdata))])]
        names(preddat) <- names(bkdat)[1:(length(names(bkdat))-1)]
        predmu$rcp26 <- predict.ppmlasso(models_list[[i]], newdata = preddat)
        # predmu <- 1-exp(-predmu) ## gives reative probabilities
        
      } else {
        preddat <- newdata[which(names(newdata) %in% 
                                   names(newdata)[-grep("45|bio", names(newdata))])]
        names(preddat) <- names(bkdat)[1:(length(names(bkdat))-1)]
        predmu$rcp85 <- predict.ppmlasso(models_list[[i]], newdata = preddat)
        # predmu <- 1-exp(-predmu) ## gives reative probabilities
      }
      rm(preddat)
    }
    return(predmu)
  }
}

# ## Bioregions layer - TO USE FOR CLIPPING...
# bioreg <- readRDS(file.path(input_dir, "bioregions_aus.rds")) 

newdata <- predxyz
bkdat <- backxyzK
RCPs <- c(45, 85)
now <- Sys.time()
prediction_list <- parallel::mclapply(seq_along(model_list), predict_ppms_apply,
                                      model_list, newdata, bkdat, RCPs, mc.cores = mc.cores)
pred_time <- Sys.time()-now
saveRDS(prediction_list, file = paste0(output_dir, "/predlist_",  gsub("-", "", Sys.Date()), ".rds"))



## ---- MOVE TO NEW SCRIPT ---- ####
## Predict and save output ####
newdata <- predxyz
bkdat <- backxyz200k
RCPs <- c(45, 85)
now <- Sys.time()
prediction_list <- parallel::mclapply(seq_along(model_list), predict_ppms_apply,
                                      model_list, newdata, bkdat, RCPs, mc.cores = mc.cores)
pred_time <- Sys.time()-now
saveRDS(prediction_list, file = paste0("./output/predlist_",  gsub("-", "", Sys.Date()), ".rds"))



## Locate errors and rerun analysis for species with errors ###
##  To be automated if error problem is not solved by species grouping
##  At the moment, it appears that error might be when species data is spatially restricted.
error_species <- read.table("./output/errorfile_1_20190725.txt", header = FALSE, sep = ",")
colnames(error_species) <- c("index", "species")
error_index <- error_species$index

names(model_list) <- tolower(gsub(" ","_", levels(factor(spdat$species))[-error_index]))
names(prediction_list) <- tolower(gsub(" ","_", levels(factor(spdat$species))[-error_index]))

spdat <- gbif
bkdat <- backxyz200k
bkwts <- bkgrd_wts
spp <- levels(factor(spdat$species))[error_index]
seq_along(spp)
mc.cores <- 1
error_models <- parallel::mclapply(1:length(spp), fit_ppms_apply, spdat, #spwts,
                                   bkdat, bkwts, interaction_terms, ppm_terms,
                                   species_names = spp, n.fits=100, min.obs = 50, mc.cores = mc.cores)
names(error_models) <- tolower(gsub(" ","_", levels(factor(spdat$species))[error_index]))

newdata <- predxyz
bkdat <- backxyz200k
RCPs <- c(26, 85)
error_pred <- parallel::mclapply(seq_along(error_models), predict_ppms_apply,
                                 error_models, newdata, bkdat, RCPs, mc.cores = mc.cores)
names(error_pred) <- tolower(gsub(" ","_", levels(factor(spdat$species))[error_index]))

errorfile <- paste0("./output/errorfile_2_", gsub("-", "", Sys.Date()), ".txt")
n <- length(model_list)+1
m <- length(errlist)+1
error_list <- list()
for (i in 1:length(error_models)){
  if(!class(error_models[[i]])[1] == "try-error") {
    model_list[[n]] <- error_models[[i]]
    n <- n+1
  }else{
    print(paste0("Model ",i, " for '", spp[i], "' has errors"))
    cat(paste(i, ",", spp[i], "\n"),
        file = errorfile, append = T)
    error_list[[m]] <- error_models[[i]]
    m <- m+1
  }
}
names(model_list)[((length(model_list)-length(error_models))+1):length(model_list)] <- names(error_models)

n <- length(prediction_list) + 1
for (i in 1:length(error_pred)) {
  prediction_list[n] <- error_pred[i]  
  n <- n + 1
}
names(prediction_list)[((length(prediction_list)-length(error_pred))+1):length(prediction_list)] <- names(error_pred)


## Catch remianing errors ####
n <- 1
m <- 1
catch_errors(seq_along(error_models), ppm_models = model_list, species_names = spp, errorfile = errorfile)





## ------ EXTRAS ----------------

## Testing for corerelations in covariates ####
## There is no hard and fast rule about how many covariates to fit to data, and it will change depending on the data and the amount of information in each observation and how they vary with the covariates, how the covariates are correlated ect... But to start I'd only fit sqrt(n_observations) covariates. So if you have 20 occurrences that's 4-5 covariates (including the polynomials) so that's only two variables! You might have to identify the most important variable and start from there.

## see slide 14 & 26 in: http://www.bo.astro.it/~school/school09/Presentations/Bertinoro09_Jasper_Wall_3.pdf

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
long.df<-reshape2::melt(df,id=c("row.names"))
pca_plot <- ggplot2::ggplot(long.df,aes(x=row.names,y=variable,fill=value))+
  geom_raster()+
  scale_fill_viridis_c()+
  theme_minimal()

## not surprising that elevation, slope and aspect are all correlated (choose one).
pca_plot


##ALTERNATIVE - Specify ppm formula based on number of observations for a species
nk <- ceiling(sqrt(nrow(spxy)))
## Specify ppm formula - independent with quadratic terms
## if low number of observations, (i.e. =< 25), just use the covariates interaction terms for the spatial variables
if(nk <= 5){ ## because linear and quad for X ad Y gives 5 terms
  ppmform <- formula(paste0(" ~ poly(", paste0(interaction_terms, collapse = ", "),
                            ", degree = 2, raw = FALSE)",collapse =""))
  ## including land-use change for species with low number of observations
  # paste0(" ~ poly(", paste0(interaction_terms, collapse = ", "),
  #        ", degree = 2, raw = FALSE) + poly(landuse, degree = 2, raw = FALSE)",collapse ="")
  
} else  {
  ## if number of observations, (i.e. > 25), fit independent with linear and quadratic terms startign with landuse
  extra_covar <- ceiling((nk - 5)/2) ## fit additional covariates based on 
  if(extra_covar > 10) extra_covar <- 10 ## if nrow(spxy) > 630 and nk > 25, only then all parameters are used!! ,...fix
  ppmform <- formula(paste0(paste0(" ~ poly(", paste0(interaction_terms, collapse = ", "),
                                   ", degree = 2, raw = FALSE)"), 
                            paste0(" + poly(", ppm_terms[!(ppm_terms %in% interaction_terms)][1:extra_covar],
                                   ", degree = 2, raw = FALSE)", 
                                   collapse = ""), collapse =""))
}

## Explore relationship between a method to specify cut-off and 3obs
n.obs <- seq(20,1000, 10)
cutoff <- ceiling(sqrt(n.obs))
add.params <- ceiling((sqrt(n.obs) - 5)/2)
max.params <- length(names(cov.mod))
obs.cutoff <- n.obs[which(add.params > max.params)[1]]
add.params[which(add.params > max.params)] <- max.params
plot(n.obs,add.params, type ="l")
abline(v=obs.cutoff)
text(obs.cutoff - 20, 4, paste0("max obs for params cut off: ", obs.cutoff), srt=90)


## 
## >> Quadrature (background) points ####
cov.mod <- stack(raster(infile))
names(cov.mod) <- c("bio1", "srtm")

global.mask <- raster(mask_file)
global.mask[which(is.na(global.mask[]))] <- 0
rpts <- rasterToPoints(global.mask, spatial=TRUE)
backxyz <- data.frame(rpts@data, X=coordinates(rpts)[,1], Y=coordinates(rpts)[,2])     
## if NAs in mask
backxyz <- backxyz[,-1]
backxyz <- na.omit(cbind(backxyz, as.matrix(cov.mod)))
## if 0s in mask
# backxyz <- backxyz[backxyz$globalmask_ee_minNA == 1,]
# backxyz <- backxyz[,-1]
backxyz200k <- backxyz[sample(nrow(backxyz), 200000), ]
backxyz200k$Pres <- rep(0, dim(backxyz200k)[1])