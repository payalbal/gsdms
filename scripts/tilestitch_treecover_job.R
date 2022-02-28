## Hansen treecover data - tile stitching

rm(list = ls())
gc()
# system("ps")
# system("pkill -f R")


x <- c('rgdal', 'gdalUtils')
lapply(x, require, character.only = TRUE)
rm(x)

data_dir <- "/home/payalb/gsdms_r_vol/tempdata/research-cifs/uom_data/gsdms_data"


## Tile stitching - gdal method
## https://stackoverflow.com/a/50235578
## Get a list of all tif tiles downloaded
tiffiles <- list.files("/home/payalb/gsdms_r_vol/tempdata/research-cifs/uom_data/gsdms_data/Hansen_treecover/earthenginepartners-hansen/GFC-2020-v1.8/", pattern="*.tif", full.names = TRUE)

## >> Build a virtual raster file stitching all tiles
## WARNING: This will fail if the file it is trying to write to (output.vrt) already exists
gdalbuildvrt(gdalfile = tiffiles, output.vrt = file.path(data_dir, "Hansen_treecover/hansen_treetemp.vrt"))

## >> Copy the virtual raster to an actual physical file
## WARNING: This takes ~5 minutes to run
gdal_translate(src_dataset = file.path(data_dir, "Hansen_treecover/hansen_treetemp.vrt"),
               dst_dataset = file.path(data_dir, "Hansen_treecover/hansen_tree.tif"),
               output_Raster = TRUE,
               options = c("BIGTIFF=YES", "COMPRESSION=LZW"))




# ## Tile stitching - raster method ####
# ## loop over each variable, collect all tiles, and mosaic.
# f = list.files("/home/payalb/gsdms_r_vol/tempdata/research-cifs/uom_data/gsdms_data/Hansen_treecover/", 
#                pattern="*.tif", full.names = TRUE)
# f = lapply(f, raster)
# 
# ## specify a function for any overlapping cells
# f$fun = mean
# mos = do.call(mosaic, f)
# 
# this_dst = paste0(path, 'file.tif') # or whatever
# writeRaster(mos, this_dst)