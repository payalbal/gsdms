
rm(list = ls())
gc()

x <- c('rgdal', 'tools', 'bitops', 'gdalUtils', 'usethis')
lapply(x, require, character.only = TRUE)
rm(x)

proj.res.km <- 10
proj_res <- proj.res.km*1000

gsdms_dir <- "/tempdata/research-cifs/6300-payalb/uom_data/gsdms_data"
output_dir <- file.path(gsdms_dir, "outputs", sprintf("layers_%sk", proj.res.km))
dst_folder <- file.path(gsdms_dir, "outputs", "distance_layers_wgs")

roadsgdb <- file.path(gsdms_dir, "groads/groads-v1-global-gdb/gROADS_v1.gdb")
roads <- readOGR(dsn=roadsgdb,layer="Global_Roads")
roadsub <- roads[roads$FCLASS == (0:2),]
roadsub <- roadsub[, 4]

writeOGR(roadsub, dsn = file.path(gsdms_dir, "groads"), layer = "groads", driver="ESRI Shapefile")
system(paste0("gdal_rasterize -at -burn 1 -ot Byte -tr .0025 .0025 ",
              file.path(gsdms_dir, "groads", "groads.shp "),
              file.path(gsdms_dir, "groads", "groads.tif")))
infile <- file.path(gsdms_dir, "groads", "groads.tif")
outfile <- file.path(dst_folder, "dst_roads.tif")
system(sprintf("gdal_proximity.py %s %s -values %d -of GTiff -distunits GEO -nodata %d -co compress=LZW -co BIGTIFF=YES", infile, outfile, 1, -9999))


