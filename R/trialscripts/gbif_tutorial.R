# GBIF DATA
# Tutorial: https://ropensci.org/tutorials/rgbif_tutorial/

library(rgbif)

# setwd("~/Users/payalb/Dropbox/Discovery/Data/")



## ---------------------------------------------------------------------------------------------------
# List of countries in GBIF dataset
isocodes
dim(isocodes)

# AUSTRALIA DATA
oz_code <- isocodes[grep("Australia", isocodes$name), "code"] # get australia code
occ_count(country=oz_code) # count records

occ_count(country=oz_code, georeferenced=TRUE) # occurrence records for Australia with groreferenced data

# to see the levels used for the different fields
ozdat <- occ_search(country=oz_code,  hasCoordinate = TRUE, basisOfRecord = "OBSERVATION", limit = 20, start = 0, 
                    fields = "all", return = "all")
str(ozdat)
names(ozdat)
names(ozdat$data)

# extract relevant columns and store data
ozdat <- occ_search(country=oz_code,  hasCoordinate = TRUE, basisOfRecord = "OBSERVATION", limit = 20, start = 0, 
                    fields = c('name','decimalLatitude','decimalLongitude'), return = "data")

names(ozdat)


# GBIF plotting
sort(unique(ggplot2::map_data("world")$region))   # find country
gbifmap(ozdat, region = "Australia")


## ---------------------------------------------------------------------------------------------------
# GLOBAL DATA
occ_count(georeferenced=TRUE) # to see the number of spatially refernced records

globdat <- occ_search(hasCoordinate = TRUE, 
                      basisOfRecord = "OBSERVATION", 
                      limit = 20, start = 0, 
                      fields = c('name','decimalLatitude','decimalLongitude'), return = "data")
names(globdat)

# for better speed use 'occ_data':
globdat <- occ_data(hasCoordinate = TRUE, 
                    basisOfRecord = "HUMAN_OBSERVATION", 
                    limit = 20, start = 0, 
                    fields = c('name','decimalLatitude','decimalLongitude'), return = "data") 

# relevant values basisOfRecord include HUMAN_OBSERVATION, LITERATURE, MACHINE_OBSERVATION, OBSERVATION
# only accepts one value for basisOfRecord at a time; to keep in mind for later


## ---------------------------------------------------------------------------------------------------
## SPECIES DATA
  # ignoring data issues, not checked for if all data is there...

# Frilled neck lizrd (Chlamydosaurus kingii)
occ_search(scientificName = "Chlamydosaurus kingii", hasCoordinate = TRUE)

dat1 <- occ_search(scientificName = "Chlamydosaurus kingii", hasCoordinate = TRUE, 
                    limit = 500, start = 0, 
                    fields = c('name','decimalLatitude','decimalLongitude'), return = "data") 
  # occ_data did not work here
dat2 <- occ_search(scientificName = "Chlamydosaurus kingii", hasCoordinate = TRUE, 
                   limit = 500, start = 501, 
                   fields = c('name','decimalLatitude','decimalLongitude'), return = "data") 
dat3 <- occ_search(scientificName = "Chlamydosaurus kingii", hasCoordinate = TRUE, 
                   limit = 500, start = 1001, 
                   fields = c('name','decimalLatitude','decimalLongitude'), return = "data") 
dat4 <- occ_search(scientificName = "Chlamydosaurus kingii", hasCoordinate = TRUE, 
                   limit = 500, start = 1501, 
                   fields = c('name','decimalLatitude','decimalLongitude'), return = "data") 

frank <- rbind(dat1, dat2, dat3, dat4)

gbifmap(frank, region = "Australia")
# does TD occur on the mainland ?
# how can you specify extent for Tasmania only ?

# OR
library(sp)
library(raster)
library(dismo)
library(maptools)
library(rgeos)

P4S <- CRS("+proj=longlat +datum=WGS84")  # the right projection and datum for the map
myclip <- extent(112,155,-45,-10)    # create a mask to remove Macquarie and Lord Howe islands
aust_bound<-crop(readShapeLines("boundaries/australia_shape/ausborder_polyline.shp", verbose=TRUE, proj4string=P4S),myclip)
state_bound<-crop(readShapeLines("boundaries/australia_shape/state_boundaries.shp", verbose=TRUE, proj4string=P4S),myclip)
plot(aust_bound)
lines(state_bound)
points(tddat[,2:3], col = "red")

# on global map
P4S <- CRS("+proj=longlat +datum=WGS84")  # for country file on gadm: http://www.gadm.org/country
glob_bound <- readShapeLines("boundaries/gadm28_shp/gadm28.shp", verbose=TRUE, proj4string=P4S)
plot(glob_bound)
points(tddat[,2:3], col = "red")


























