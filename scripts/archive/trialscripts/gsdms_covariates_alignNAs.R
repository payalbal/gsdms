## CHECK AND CORRECT FOR NAS IN DATA


## Set working environment
library(sp)
library(raster)
library(viridis)
library(viridisLite)

data_raw <- "/Volumes/discovery_data/data/raw"
data_gsdms <- "/Volumes/discovery_data/data/gsdms_data"


## Count occurrence points falling off the mask (edge) - Run on boab
# NAcounts <-data.frame(matrix(NA, length(gbif_sp), 5))
# colnames(NAcounts) <- c("# records", "extract", "extract.by.buffer", "extract.by.bilinear", "NAcovariates")
# for (i in 1:length(gbif_sp)) {
#   spxy <- gbif[gbif$species == gbif_sp[i], c(4,3)]
#   NAcounts[i,1] <- dim(spxy)[1]
#   spxyz1 <- extract(covariates, spxy)
#   NAcounts[i,2] <- sum(rowSums(is.na(spxyz1))!=0)
#   spxyz2 <- extract(covariates, spxy, buffer = 1000000, small = TRUE, fun = mean, na.rm = TRUE)
#   NAcounts[i,3] <- sum(rowSums(is.na(spxyz2))!=0)
#   spxyz3 <- extract(covariates, spxy, method = 'bilinear')
#   NAcounts[i,4] <- sum(rowSums(is.na(spxyz3))!=0)
#   NAcounts[i,5] <- paste(dimnames(spxyz1)[[2]][colSums(is.na(spxyz1)) != 0], collapse = ", ")
# }
# NAcounts$species <- gbif_sp
# NAcounts
# saveRDS(NAcounts, file="NAcounts.rds")

## Note: extracting by buffer gets rid of more NAs than by using bilinear method, 
##  because buffer considers larger area around the point

readRDS("./output/NAcounts.rds")


## Load global_mask and covariates
global_mask <- raster(".data/bio1.bil") # 10min resolution ~ 345km2 (source: climate/wc10/bio1.bil")
global_mask[which(!is.na(global_mask[]))] <- 1
covariates <- readRDS("./data/covariates_all.rds")


## Find and plot NAs in soil variables: "bulkdens", "pawc", "soilcarb", "totaln"
plot(global_mask)
plot(covariates[["bulkdens"]], col = viridis(100))

  ## Index of non-NA values in global_mask
  nona <- which(!is.na(values(global_mask)))

    ## CHECK: Compare indices (nas) across all soil variables 
    ##    - indicates these are same across soil variables
    nabulkdens <- which(is.na(values(covariates[["bulkdens"]])[nona]))
    napawc <- which(is.na(values(covariates[["pawc"]])[nona]))
    nasoilcarb <- which(is.na(values(covariates[["soilcarb"]])[nona]))
    natotaln <- which(is.na(values(covariates[["totaln"]])[nona]))
    paste0("# NAs in ", c('bulkdens', 'pawc', 'soilcarb', 'totaln'), " = ", 
           c(length(nabulkdens), length(napawc), length(nasoilcarb), length(natotaln)))
    sum(nabulkdens == napawc)
    sum(napawc == nasoilcarb)
    sum(nasoilcarb == natotaln)
  
  ## Index of NA values in soil variables where global_mask is non-NA
  nas <- which(is.na(values(covariates[["bulkdens"]])[nona]))
  
  ## Plot 'nas' on map
  xys <- xyFromCell(global_mask, nona)
  head(xys)
  xys <- xys[nas,]
  dim(xys)
  plot(global_mask)
  points(xys, pch = 21, cex = .2)

  ## Remove these values from global_mask
  values(global_mask)[cellFromXY(global_mask, xys)] <- NA
  
  ## CHECK: NA values have been removed from global_mask
  nona <- which(!is.na(values(global_mask)))
  nas <- which(is.na(values(covariates[["bulkdens"]])[nona]))
  xys <- xyFromCell(global_mask, nona)
  xys <- xys[nas,]
  head(xys)
  paste0("No NAs found = ", dim(xys)[1] == 0)
  plot(global_mask)
  points(xys, pch = 21, cex = .2)
  
  
## Find and plot NAs in srtm variables: "srtm", "slope", "roughness", "aspect"
plot(global_mask)
plot(covariates[["srtm"]])
plot(covariates[["slope"]], col = viridis(100))
  
  ## Index of non-NA values in global_mask
  nona <- which(!is.na(values(global_mask)))
  
    ## CHECK: Compare indices (nas) across soil variables 
    ##    - indicates these are same across soil variables
    nasrtm <- which(is.na(values(covariates[["srtm"]])[nona]))
    naslope <- which(is.na(values(covariates[["slope"]])[nona]))
    naroughness <- which(is.na(values(covariates[["roughness"]])[nona]))
    naaspect <- which(is.na(values(covariates[["aspect"]])[nona]))
    paste0("# NAs in ", c('srtm', 'slope', 'roughness', 'aspect'), " = ", 
           c(length(nasrtm), length(naslope), length(naroughness), length(naaspect)))
    sum(nasrtm == naslope)
    sum(naslope == naroughness)
    sum(naroughness == naaspect)
  
  ## Index of NA values in slope/roughnness/aspect variables where global_mask is non-NA
  nas <- which(is.na(values(covariates[["slope"]])[nona]))
  
  ## Plot nas on map
  xys <- xyFromCell(global_mask, nona)
  head(xys)
  xys <- xys[nas,]
  dim(xys)
  plot(global_mask)
  points(xys, pch = 21, cex = 1)
  
  ## Remove these values from global_mask
  values(global_mask)[cellFromXY(global_mask, xys)] <- NA
  
  ## CHECK: NA values have been removed from global_mask
  nona <- which(!is.na(values(global_mask)))
  nas <- which(is.na(values(covariates[["slope"]])[nona]))
  xys <- xyFromCell(global_mask, nona)
  xys <- xys[nas,]
  head(xys)
  paste0("No NAs found = ", dim(xys)[1] == 0)
  plot(global_mask)
  points(xys, pch = 21, cex = .2)
  

## Save updated global_mask
covariates_new <- raster::mask(covariates, mask = global_mask)
summary(covariates)
summary(covariates_new)





