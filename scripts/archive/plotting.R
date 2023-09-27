### PPM plots for IUCN species

prediction_list <- readRDS("~/Dropbox/discovery_trade/analyses/gsdms/output/predlist_20190708.rds")


## Summed
mu_current <- matrix(NA, length(prediction_list[[1]]$rcp26), length(prediction_list))
colnames(mu_current) <- names(prediction_list)
mu_rcp26 <- matrix(NA, length(prediction_list[[1]]$rcp26), length(prediction_list))
colnames(mu_rcp26) <- names(prediction_list)
mu_rcp85 <- matrix(NA, length(prediction_list[[1]]$rcp26), length(prediction_list))
colnames(mu_rcp85) <- names(prediction_list)

for (i in 1:length(prediction_list)) {
  mu_current[,i] <- prediction_list[[i]]$current
  mu_rcp26[,i] <- prediction_list[[i]]$rcp26
  mu_rcp85[,i] <- prediction_list[[i]]$rcp85
}

## Site based indices
richness_current <- rowSums(mu_current)
pred_current <- rasterFromXYZ(cbind(predxyz[,1:2],richness_current))
richness_rcp26 <- rowSums(mu_rcp26)
pred_rcp26 <- rasterFromXYZ(cbind(predxyz[,1:2],richness_rcp26))
richness_rcp85 <- rowSums(mu_rcp85)
pred_rcp85 <- rasterFromXYZ(cbind(predxyz[,1:2],richness_rcp85))

## Plot richess for RCP26 and RCP85
plot(pred_rcp26, main = "@rcp26", ext = extent(global_mask), col = viridisLite::viridis(10), breaks=seq(minValue(pred_rcp26), maxValue(pred_rcp26), length.out = 50), axes=F, box=F, legend = FALSE)
plot(pred_rcp85, main = "@rcp85", ext = extent(global_mask), col = viridisLite::viridis(10), breaks=seq(minValue(pred_rcp85), maxValue(pred_rcp85), length.out = 50), axes=F, box=F, legend = FALSE)

## copy legend from
# plot(pred_rcp26, main = "@rcp26", ext = extent(global_mask), col = viridisLite::viridis(10), breaks=seq(minValue(pred_rcp26), maxValue(pred_rcp26), length.out = 11), axes=F, box=F)

## PLot differece in RCP scenarios
plot(overlay(pred_current, pred_rcp26, fun=function(x,y) as.logical(round(x,2)==round(y,2))), col= c("salmon", "steelblue"), box=FALSE, axes=FALSE, legend = FALSE, main = "Changed expected under RCP 26")
legend(30, -35, legend = c("change", "no change"), fill=c("salmon", "steelblue"), bty = "n", cex = 1)

plot(overlay(pred_current, pred_rcp85, fun=function(x,y) as.logical(round(x,2)==round(y,2))), col= c("salmon", "steelblue"), box=FALSE, axes=FALSE, legend = FALSE, main = "Changed expected under RCP 85")
legend(30, -35, legend = c("change", "no change"), fill=c("salmon", "steelblue"), bty = "n", cex = 1)


## By species

species <- names(prediction_list)

for (i in 1:length(prediction_list)){
  
  png(file = paste0("./output/plot_",tolower(gsub(" ","",species[i])), ".png"),
      width=1500, height=600, res = 130, pointsize = 12)
  par(mfrow = c(1,2))
  par(mar=c(3,3,5,7))
  
  preds_sp <- rasterFromXYZ(cbind(predxyz[,1:2],prediction_list[[i]]))
  plot(preds_sp$rcp26, main = "rcp26", ext = extent(global_mask))
  points(spdat[spdat$species == species[i], c(4,3)], pch=16, cex=.3)
  plot(preds_sp$rcp85, main = "rcp85", ext = extent(global_mask))
  points(spdat[spdat$species == species[i], c(4,3)], pch=16, cex=.3)
  mtext(species[i], side = 3, line = -2, outer = TRUE, cex = 1.5, font = 2)
  rm(preds_sp)
  dev.off()
}




