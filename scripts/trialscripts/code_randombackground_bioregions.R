## Generate random background points based on specified region
## Author: Alejandra Moran-Ordonezlibrary(raster)
library(dismo)
library(rgdal)
library(maptools)


IB
RA<-readShapeSpatial(paste(wd, "gridsARC/IBRA7_subregions.shp", sep="")) ##  read POLY shapefile with all BIOREGIONS
IBRA.grid<-mask.clim## creates an empty raster using another raster as 'model'/mask (in my case my reference raster was 'mask.clim')
IBRA.grid[]<- NA ## it changes all IBRA.grid values to NA
IBRA.grid<-rasterize(IBRA, IBRA.grid, field="OBJECTID") ## the grid takes the values of the ID IBRA region (bioregion)
#   writeRaster(IBRA.grid,paste0(wd,"gridsARC/IBRA_grid_maskclim.tif"),overwrite=TRUE) 

### Now we need to know the spatial relationships between bioregions (who is neighbour of who)
### I created a table using the "Polygon Neighbours tool"  in ARCGIS
### ARCGIS> Analysis tools> Proximity> Polygon Neigbours.
## Column "src_OBJECT" corresponds to the ID values for each polygon in the shapefile and "nbr_OBJECT" gives the unique ID values for their neighbour polygons
### I?ve attached the table I used for Australia for your reference

IBRA.neigh<-read.csv(paste(wd,"gridsARC/IBRA_neighbours.csv", sep="")) 
IBRA.neigh<-IBRA.neigh[, 1:2] ## extracts the two first columns of the table
colnames(IBRA.neigh)<-c("OBJECTID", "NEIGHID") ## Change the names of the columns


### Now, I have my data set 'mammals.AU.record', Colums 8 and 9 are the 'long' 'lat' coordinates of my species records
mammalsAU.record$IBRA<-extract(IBRA.grid,mammalsAU.record[,8:9]) ## Creates a column in my data set indicating the corresponding IBRA regionfor each species presence record 


for(i in 1:length(sp.list)) ## species.list is a vector with the names of all the species you are working with
{
  ## From IBRA.neight, subset the IBRA regions with presence of species i
  presence.poly <- IBRA.neigh[IBRA.neigh$OBJECTID %in% unique(mammalsAU.record[mammalsAU.record$species==sp.list[i],'IBRA']),] ## subset tregion with presence of species i 
  ## vector with all unique values for the neighbour IBRA regions of those with presence data.
  neigh.poly <-as.vector(unique(presence.poly$NEIGHID))
  ## vector with unique values of IBRA regions with presence data
  presence.poly <-as.vector(unique(presence.poly$OBJECTID)) 
  
  all.poly<-unique(c(presence.poly, neigh.poly)) ## All IBRA regions= presence + neighbours of interest to generate background points
  
  z <-IBRA.grid %in% all.poly ## z is a raster with a subset of IBRA regions of interest for species X
  
  g<- IBRA.grid %in% presence.poly #if we ONLY want to consider as background the IBRA regions with presence recordss
  
  ## reclassifies Z raster to a binary map 1 background 0 area of no interest for the species
  m<-c(0,NA, 1, 1)
  rclmat <- matrix(m, ncol=2, byrow=TRUE)
  z <- reclassify(z, rclmat, include.lowest= T) 
 
  ## reclassifies Z raster to a binary map 1 background 0 area of no interest for the species
  g <- reclassify(g, rclmat, include.lowest= T) 
  
  ## this creates a 5000 random points in the IBRA regions selected area for each species (z)
  back<-randomPoints(z, 5000,p=mammalsAU.record[mammalsAU.record$species==sp.list[i],c("long_NEW","lat_NEW")],excludep=T)

  back<-data.frame(rep("background", 5000),back)
  names(back)[1:3] <- c("Background", "X", "Y") ## this creates a table with long and lats for the 5000 random points that you are going to use in MaxEnt. You can use this lat long data to extract the environmental covariates values at those points
  
  cc<-extract(pred.ext, back[,2:3])   ## extract  predictor values on background locations (pred.ext is a rasterStack of environmental predictors)
  back.sp<-cbind(back[,],cc)
  back.sp<-na.omit(back.sp) ## remove potential NA in the data frame
  write.table(back.sp, paste0(wd,"Extremes_analyses/mammals_models/background/IBRAtest/", sp.testname [i],"_ext_IBRAneigh.csv"), 
              sep=",", col.names= T, row.names=F) ## this saves the background dataframe  of each species into a csv file
  ### You could also repeat lines 52-61 using as reference the raster 'g' so you create random background points in the bioregion where there are species occurrence data
  rm(g, z, cc, back, back.sp); gc()
  }

