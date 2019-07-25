## SPECIES WITH ERROR MODELS

pacman::p_load(sp, raster)

raw_mask <- raster("./data/bio1.bil") 
raw_mask[which(!is.na(raw_mask[]))] <- 1
gbif <- read.csv("./data/2019-05-14_gbif_iucnsp.csv", header = TRUE)

spdat <- gbif
spp <- levels(factor(spdat$species))

## Species with error in ppm models
err_spp <- spp[c(11,12,18)]
err_spp <- tolower(gsub(" ", "", err_spp))
nangerdama <- spdat[which(spdat$species == spp[11]),]
dim(nangerdama)
notomysaquilo <- spdat[which(spdat$species == spp[12]),]
dim(notomysaquilo)
phascogalepirata <- spdat[which(spdat$species == spp[18]),]
dim(phascogalepirata)

plot(raw_mask, legend = FALSE)
points(nangerdama[,c(4,3)], col = "black", cex = 0.5)
points(notomysaquilo[,c(4,3)], col = "red", cex = 0.5)
points(phascogalepirata[,c(4,3)], col = "blue", cex = 0.5)        


## Plot all
for (i in spp){
  dat <- spdat[which(spdat$species == i),]
  plot(raw_mask, col = "grey", legend = FALSE)
  points(dat[,c(4,3)], col = "blue", cex = 0.3)
  title(paste(i, "\n", "#records = ", dim(dat)[1]))
}


## Restricted
