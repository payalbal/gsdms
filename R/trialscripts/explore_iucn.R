## Opening IUCN mammals range shapefile 

library(sp)
library(maptools) # to work with shapefiles
library(rgdal)
library(raster)

## Load data
data_path <- "/Volumes/payal_umelb/data/"
data <- readOGR(dsn=paste0(data_path, "raw/biodiversity/IUCN/rangemaps/TERRESTRIAL_MAMMALS"), layer = "TERRESTRIAL_MAMMALS")
n_species <- length(unique(data@data$binomial)) # number of species
names(data)
# data_spatial <- readShapeSpatial(file.path(dat_wd,'/TERRESTRIAL_MAMMALS'))

## Filter data
## By threatened status: include VU, EN, CR; exclude EW, EX and DD (because spatial data won't be be availabel for these)spatial info might be sparse); exclude NT, LC
unique(data$code)
dat <- data[which(data$code == "VU" | data$code == "EN" | data$code == "CR"),] #2581

## Remove subspecies  
sum(is.na(dat$binomial))
sum(sapply(strsplit(as.character(dat$binomial), " "), length) < 2)
sum(sapply(strsplit(as.character(dat$binomial), " "), length) == 2)
sum(sapply(strsplit(as.character(dat$binomial), " "), length) > 2)

dat <- dat[which(sapply(strsplit(as.character(dat$binomial), " "), length) == 2), ] #2578

## Remove species classified as marine
length(unique(dat[which(dat$marine == "t"),]$binomial))
dim(dat[which(dat$marine == "t"),])

dat<- dat[which(dat$marine == "f"), ] #2555
  
length(unique(dat[which(dat$freshwater == "t"),]$binomial))
dim(dat[which(dat$freshwater == "t"),])

dat <- dat[which(dat$freshwater == "f"),] #2480
length(unique(dat$binomial))

iucn_species <- unique(dat$binomial)

## Additional filters
  ##Only keep native species
dim(dat[which(dat$origin == 1),])
dat <- dat[which(dat$origin == 1),] #2378
  ## Remove extinct
dim(dat[which(dat$presence == 4 | dat$presence == 5 | dat$presence == 6),])
dat <- dat[which(dat$presence == 4 | dat$presence == 5 | dat$presence == 6),]
iucn_species2 <- unique(dat$binomial)


## Save species list
# write.table(iucn_species, "./output/iucn_species.txt",  row.names = FALSE, col.names = FALSE)
write.table(iucn_species2, "./output/iucn_species.txt",  row.names = FALSE, col.names = FALSE)


## Shapefile per species - SAVED ON DISK
## Split One Shapefile into Many (Useful for IUCN Spatial Data)
## http://rfunctions.blogspot.com.au/2013/03/spatial-analysis-split-one-shapefile.html

# ##Loop through species and save separate shapefile for each species
# species <- unique(data@data$binomial)
# for (i in 1:length(species)) {
#   tmp <- data[data$binomial == species[i], ] 
#   writeOGR(tmp, dsn=paste0(data_path, 'raw/biodiversity/IUCN/rangemaps/by_species', sep=""), species[i], driver="ESRI Shapefile",
#            overwrite_layer=TRUE)
# }











