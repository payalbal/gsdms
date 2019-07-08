### PPM plots for IUCN species

gbif <- read.csv("./data/2019-05-14_gbif_iucnsp.csv", header = TRUE)
spp <- levels(factor(gbif$species))
error_sp <- read.csv("./output/errorfile_20190528.txt", header = FALSE, sep = ",")
colnames(error_sp) <- c("sp_idx", "sp_name")
error_idx <- error_sp$sp_idx
pred_spp <- spp[!spp %in% trimws(as.vector(error_sp$sp_name))]


## Load covariate data
covariates <- readRDS("./data/covariates.rds")
covariates_predict <- readRDS("./data/covariates_predict.rds")
raw_mask <- raster("./data/bio1.bil") 

## Mask
raw_mask[which(!is.na(raw_mask[]))] <- 1
source("./R/align.maskNA.R")
global_mask <- align.maskNA(covariates, raw_mask)
covariates <- mask(covariates, global_mask)
covariates_predict <- mask(covariates_predict, global_mask)
# covariates_predict_rcp26 <- covariates_predict[[names(covariates_predict)[grep('26', names(covariates_predict), invert = TRUE)]]]
# covariates_predict_rcp85 <- covariates_predict[[names(covariates_predict)[grep('26', names(covariates_predict), invert = TRUE)]]]

## Prediction grid
global_mask0 <- global_mask
global_mask0[which(is.na(global_mask0[]))] <- 0
rpts <- rasterToPoints(global_mask0, spatial=TRUE)
predxyz <- data.frame(rpts@data, X=coordinates(rpts)[,1], Y=coordinates(rpts)[,2])     
predxyz <- predxyz[,-1]
predxyz <- na.omit(cbind(predxyz, as.matrix(covariates_predict)))


predlist <- readRDS("~/Dropbox/discovery_trade/analyses/gsdms/output/predlist_20190528.rds")
prediction_list <- predlist

mu_rcp26 <- matrix(NA, length(prediction_list[[1]]$rcp26), length(pred_spp))
colnames(mu_rcp26) <- pred_spp
mu_rcp85 <- matrix(NA, length(prediction_list[[1]]$rcp26), length(pred_spp))
colnames(mu_rcp85) <- pred_spp

for (i in 1:length(pred_spp)) {
  mu_rcp26[,i] <- prediction_list[[i]]$rcp26
  mu_rcp85[,i] <- prediction_list[[i]]$rcp85
}

## Species based indices
colSums(mu_rcp26)...

## Site based indices
site_richness26 <- rowSums(mu_rcp26)
pred_rcp26 <- rasterFromXYZ(cbind(predxyz[,1:2],site_richness26))
site_richness85 <- rowSums(mu_rcp85)
pred_rcp85 <- rasterFromXYZ(cbind(predxyz[,1:2],site_richness85))

## Plot richess for RCP26 and RCP85
plot(pred_rcp26, main = "@rcp26", ext = extent(global_mask), col = viridisLite::viridis(10), breaks=seq(minValue(pred_rcp26), maxValue(pred_rcp26), length.out = 50), axes=F, box=F, legend = FALSE)
plot(pred_rcp85, main = "@rcp85", ext = extent(global_mask), col = viridisLite::viridis(10), breaks=seq(minValue(pred_rcp85), maxValue(pred_rcp85), length.out = 50), axes=F, box=F, legend = FALSE)

## copy legend from
# plot(pred_rcp26, main = "@rcp26", ext = extent(global_mask), col = viridisLite::viridis(10), breaks=seq(minValue(pred_rcp26), maxValue(pred_rcp26), length.out = 11), axes=F, box=F)

## PLot differece in RCP scenarios
temp <- overlay(pred_rcp26, pred_rcp85, fun=function(x,y) as.logical(round(x,2)==round(y,2)))
plot(overlay(pred_rcp26, pred_rcp85, fun=function(x,y) as.logical(round(x,2)==round(y,2))), col= c("salmon", "steelblue"), box=FALSE, axes=FALSE, legend = FALSE)
legend(30, -35, legend = c("change", "no change"), fill=c("salmon", "steelblue"), bty = "n", cex = 1)
