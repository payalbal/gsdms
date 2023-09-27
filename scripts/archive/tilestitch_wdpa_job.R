## WDPA protected area data - rasterize & stich tiles

rm(list = ls())
gc()
# system("ps")
# system("pkill -f R")


x <- c('rgdal', 'gdalUtils')
lapply(x, require, character.only = TRUE)
rm(x)

data_dir <- "/home/payalb/gsdms_r_vol/tempdata/research-cifs/uom_data/gsdms_data"



## Rasterise ####

## >> Get a list of all shp tiles downloaded
shpfiles <- list.files("/home/payalb/gsdms_r_vol/tempdata/research-cifs/uom_data/gsdms_data/protectedareas",
                       pattern="polygons", full.names = TRUE, recursive = TRUE)
shpfiles <- grep(".shp$", shpfiles, value = TRUE)


for (infile in shpfiles){
  ## >> Specify output filename
  outfile <- gsub(".shp$", ".tif", infile)

  ## >> Rasterize @ .0025 ~ 250m
  system(paste0("gdal_rasterize -at -burn 1 -ot Byte -tr .0025 .0025 ",
                infile, " ",
                outfile))

  ## >> Get tile number
  x <- as.integer(sub(".*?shp_.*?(\\d+).*", "\\1", outfile))

  ## >> Rename file
  file.rename(outfile, gsub(".tif$", paste0(x, ".tif"), outfile))

}



## Tile stitching - gdal method ####
## https://stackoverflow.com/a/50235578

## >> Get a list of all tif tiles
tiffiles <-  list.files("/home/payalb/gsdms_r_vol/tempdata/research-cifs/uom_data/gsdms_data/protectedareas", 
                        pattern=".tif$", full.names = TRUE, recursive = TRUE)
  
## >> Build a virtual raster file stitching all tiles
## WARNING: This will fail if the file it is trying to write to (output.vrt) already exists
gdalbuildvrt(gdalfile = tiffiles, output.vrt = file.path(data_dir, "protectedareas/wdpa_temp.vrt"))

## >>  Copy the virtual raster to an actual physical file
## WARNING: This takes ~5 minutes to run
gdal_translate(src_dataset = file.path(data_dir, "protectedareas/wdpa_temp.vrt"),
               dst_dataset = file.path(data_dir, "protectedareas/wdpa.tif"),
               output_Raster = TRUE,
               options = c("BIGTIFF=YES", "COMPRESSION=LZW"))


