

## -----------------------------------------------------------------------------  
## EXTRAS

## FIXING LAMBDAS vs. N.FITS
## see: https://web.stanford.edu/~hastie/glmnet/glmnet_alpha.html
## Notes: Defining lambda sequence externally 
##  (using skip’s code to fix n.fits =20 and max.lambda = 100) takes ~3 times 
##  less time than specifying n.fits = 100 to generate lambdas internally.
##  lambdaseq ~6 mins for Panthera leo (~8 on boab) & 
##  ~16 mins for 100 n.fits (~25 on boab)

# mod.by.nfits <- ppmlasso(formula = ppmform, data = ppmxyz, 
#                          n.fits = 20)

## CLIP OBSERVATIONS
## Method 1: by IUCN range maps - XX
## Method 2: 200 k m buffer around observation point (convex hull)
## Method 3: by Biogeoghraphic zones (Olson et al. 2001)

## Method 1
dat <- rgdal::readOGR(dsn=paste0(data_path, "raw/biodiversity/IUCN/rangemaps/TERRESTRIAL_MAMMALS"), 
                      layer = "TERRESTRIAL_MAMMALS") ## as downloaded from IUCN
source("./R/filter.iucn.R")
iucn_range <- filter.iucn(dat, filter_fields = c("binomial"),
                          species_list = FALSE, write_file = FALSE)
iucn_range <- iucn_range[which(iucn_range$binomial %in% gbif_sp),]
as.matrix(table(as.vector(iucn_range$binomial)))
sprange <- sf::st_as_sf(iucn_range[which(iucn_range$binomial == i),])
spxy_spatial <- sf::st_as_sf(spxy, coords = c("X", "Y"), crs = 4326)
sp_inrange <- st_intersection(spxy_spatial, sprange) # or: spxy_spatial[sprange,]
dim(sp_inrange)[1]/dim(spxy_spatial)[1]


# plot(sprange, axes = TRUE)
# plot(global_mask)
# plot(sprange, axes = FALSE, add = TRUE, col = "black")
plot(sprange, axes = TRUE)

plot(global_mask)
plot(sprange, axes=FALSE, add = TRUE)

# plot(spxy_spatial, col = "red", axes = TRUE)
plot(spxy_spatial, col = "red", axes = FALSE, add = TRUE, pch = 16, cex = 0.3)
plot(spxy_spatial[sprange,], col = "yellow", axes = FALSE, add = TRUE, pch = 20)


length(iucn_range[which(iucn_range$binomial == i),]@polygons)
str(iucn_range[which(iucn_range$binomial == i),]@polygons)

length(iucn_range[which(iucn_range$binomial == "Pipistrellus maderensis"),]@polygons)
str(iucn_range[which(iucn_range$binomial == "Pipistrellus maderensis"),]@polygons)


## ALTERNATE METHOD TO GENERATE spxy and env.grid using ppmdat()
scales <- c(0.5, 1, 2, 4, 8, 16)
## Find resolution for analysis
out <- findres(scales, formula = ppmform, spxy = spxyz, env.grid = backxyz)
## Generate weights
ppmxyz2 <- ppmdat(spxy = spxyz, sp.scale = 1, back.xy = backxyz, coord = c("X", "Y"))
## code from Skip on generating weights


## PLOTTING
quadraster <- rasterFromXYZ(quad)
plot(quadraster[[1]], col = "black", maxpixels = 5e100)
length(which(!is.na(values(quadraster))))

quad_spatial <- rasterToPoints(quadraster[[1]])
quad_spatial <- SpatialPointsDataFrame(coords = quad[complete.cases(quad[,1:3]),1:2], 
                                       data = quad[complete.cases(quad[,1:3]),1:3],
                                       proj4string = crs(global_mask))
points(quad_spatial, cex = 0.1, pch = 16)
crs(quadraster) <- crs(global_mask)
points(quad_spatial, pch = 16, col= "red")
points(spxy, pch = 20, cex= 0.4)


## sf commands from Roozbeh's tutorial
dat <- gbif[1:500, 2:4]
data_sf <- st_as_sf(dat, coords = c("decimallongitude", "decimallatitude"), crs = 4326)
plot(data_sf, axes = TRUE)
plot(data_sf["species" == "Panthera pardus"], axes = TRUE)
plot(st_geometry(data_sf), axes = TRUE)
as(data_sf, "Spatial")

library(rworldmap)
map <- getMap(resolution = "low")
# aus <- map[which(map$NAME == "Australia")]
aus <- st_as_sf(aus)
map <- st_as_sf(map)
plot(st_geometry(map))