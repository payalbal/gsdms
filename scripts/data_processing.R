## Covariate data processing for global analyses
## NOTE: Data sources recorded in data_downloads.R


## ------------------------------------------------------------------------- ##
## Set working environment ####
## ------------------------------------------------------------------------- ##

rm(list = ls())
gc()
# system("ps")
# system("pkill -f R")


## Load libraries
# install.packages("remotes")
# remotes::install_github("skiptoniam/geostuff")
x <- c('data.table', 'sp', 'raster', 'terra',
       'rgdal', 'geostuff', 'tools', 'bitops', 
       'ncdf4', 'processNC',
       'RCurl', 'gdalUtils', 'usethis',
       'parallel', 'doMC',
       'reticulate')
lapply(x, require, character.only = TRUE)
rm(x)


## >> Specify global variables ####
proj.res.km <- 10
mc.cores = future::availableCores()-2

luvars <- c("primf", "primn", "secdf", "secdn", "urban", "pastr", "range", "crop" )

## Folder paths
data_dir <- "/home/payalb/gsdms_r_vol/tempdata/research-cifs/uom_data/gsdms_data"
output_dir <- file.path(data_dir, "outputs", sprintf("layers_%sk", proj.res.km))
if(!dir.exists(output_dir)) {
  dir.create(output_dir)
}


## Functions
source("/home/payalb/gsdms_r_vol/tempdata/workdir/gsdms/scripts/0_functions.R")
source("/home/payalb/gsdms_r_vol/tempdata/workdir/landuse_projects/land-use/gdal_calc.R") # by jgarber

  # ## To pull directly from jgarber's repo:
  # library(RCurl)
  # url <- "https://gitlab.unimelb.edu.au/garberj/gdalutilsaddons/-/blob/master/gdal_calc.R"
  # url.exists(url)
  # source(url)
  # devtools::source_url(url)
  # url <- "https://gitlab.unimelb.edu.au/garberj/gdalutilsaddons.git"
  # devtools::install_git(url = url)




## ------------------------------------------------------------------------- ##
## Create mask ####
## ------------------------------------------------------------------------- ##

## step one_create mask from WorldClim layer 
## crs for Worldclim layers: https://worldclim.org/data/v1.4/formats.html
bio_current <- list.files(paste0(data_dir, "/bio_30s"), pattern = "bio_current*", full.names = TRUE)
global_mask <- raster(bio_current[1])
global_mask[which(!is.na(global_mask[]))] <- 1
wgs_crs = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0" 
crs(global_mask) <- wgs_crs
writeRaster(global_mask, filename = file.path(output_dir, "globalmask_wgs_1degree.tif"))


## step two_clip mask to desired extent
## Note: Clipping has to be done before the next step for proj and reso
infile <- file.path(output_dir, "globalmask_wgs_1degree.tif")
outfile <- file.path(output_dir, "globalmask_clip.tif")
new_extent <- "-180 -60 180 90"
system(paste0("gdalwarp -overwrite -ot Byte",
              " -te ", new_extent, " ",
              infile, " ", outfile))
raster(infile)
raster(outfile)


## step three_reproject mask to Equal Earth and desired resolution
## Ref for EE: https://proj.org/operations/projections/eqearth.html#id2
## jgarber note the '+proj=eqearth' is a new format and wont work wth gdalUtils
## instead we use '+proj=eqearth +ellips=WGS84 +wktext' which give equal earth and works here
## see: https://en.wikipedia.org/wiki/Equal_Earth_projection
infile <- outfile
outfile <- gsub("_clip", "_ee", infile)
new_res <- c(proj.res.km*1000, proj.res.km*1000)
# new_res <- c(1000,1000) ## 1km res
new_crs = '+proj=eqearth +ellips=WGS84 +wktext'
system(paste0("gdalwarp -overwrite -ot Byte -tr ",
              paste(new_res, collapse = " "),
              " -t_srs '", new_crs, "' ",
              infile, " ", outfile))
# gdalwarp(outfile, outfile2, s_srs = wgs_crs, t_srs = new_crs, tr = reso_ee) #, output_Raster =TRUE, verbose=TRUE) #these extra arguments were used to see what was going on and the output
raster(infile)
raster(outfile)
unique(values(raster(outfile)))
plot(raster(outfile))

file.remove(infile)


## step four_change NoData values to -9999
## ref: https://gis.stackexchange.com/questions/298230/change-no-data-value-geotif-file-with-qgis-or-gdal
change_values = "/home/payalb/gsdms_r_vol/tempdata/workdir/gsdms/scripts/change_values.sh"
system(paste0("bash ", change_values, " ", outfile, " -9999"))

gdalUtils::gdalinfo(outfile)
gdalUtils::gdalinfo(file.path(output_dir, "globalmask_ee_edt.tif"))
unique(values(raster(file.path(output_dir, "globalmask_ee_edt.tif"))))
plot(raster(file.path(output_dir, "globalmask_ee_edt.tif")))

file.rename(file.path(output_dir, "globalmask_ee_edt.tif"), 
            file.path(output_dir, sprintf("globalmask_%sk_ee.tif", proj.res.km)))

file.remove(file.path(output_dir, "globalmask_ee_temp.tif"))
file.remove(outfile)





## ------------------------------------------------------------------------- ##
## Prepare Hurtt layers ####
## ------------------------------------------------------------------------- ## 
## https://luh.umd.edu/
## 2015 - 2100 (86 time steps)

hurtt_out <- file.path(data_dir, "outputs", "hurtt")
if(!dir.exists(hurtt_out)) {
  dir.create(hurtt_out)
}
# file.remove(list.files(hurtt_out, full.names = TRUE))


## Explore Hurtt input files
infiles = list.files(file.path(data_dir, "landuse/Hurtt"), 
                     pattern = ".nc$", full.names = TRUE)

gdalUtils::gdalinfo(infiles[1])
lu <- ncdf4::nc_open(infiles[1], write=FALSE, readunlim=FALSE, 
                     verbose=FALSE,auto_GMT=TRUE, 
                     suppress_dimvals=FALSE)
names(lu)
names(lu$var)

temp <- ncvar_get(temp, varid = "primf", 
                  start = NA, count = NA, 
                  verbose=FALSE,
                  signedbyte=TRUE, 
                  collapse_degen=TRUE, 
                  raw_datavals=TRUE )
temp <- plot(raster(temp2))


## Specify time steps and index
time.all <- seq(2015, 2100, 1)
time.steps <- seq(2015, 2070, 5)
time.idx <- which(t.all %in% t.steps)


## >> step one_Extract layers by ssp, land use and year ####
system("bash /tempdata/workdir/gsdms/scripts/extract_nc.sh")

  ## Check files
  length(list.files(hurtt_out))
  x <- list.files(hurtt_out, full.names = TRUE)
  raster(x)
  plot(raster(x))


## >> step two_Rename files by year of time slice ####
for (i in 1:length(t.idx)){
  x <- list.files(hurtt_out, full.names = TRUE, pattern = paste0("t", t.idx[i], ".nc$"))
  y <- gsub(paste0("t", t.idx[i], ".nc$"), paste0("t", t.steps[i], ".nc"), x)
  file.rename(x, y)
}


## >> step three_Aggregate crop layers ####
system("bash /tempdata/workdir/gsdms/scripts/max_nc.sh")

  
## >> step four_Reproject to Equal Eath ####
infile <- list.files(hurtt_out, 
                     pattern = paste0(luvars, collapse = "|"), 
                     full.names = TRUE)
length(infile) == length(luvars)*12*3 
  ## number of variables *  number of time steps * number of scenarios
gdalUtils::gdalinfo(infile[1])
raster(infile[1])

outfile <- file.path(output_dir, paste0("hurtt_", basename(tools::file_path_sans_ext(infile)), ".tif"))
new_res <- c(proj.res.km * 1000, proj.res.km * 1000)
new_crs = '+proj=eqearth +ellips=WGS84 +wktext'
wgs_crs = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0" 

parallel::mclapply(seq_along(infile),
                   function(x) system(paste0("gdalwarp -overwrite -ot Byte -r bilinear -tr ",
                                             paste(new_res, collapse = " "),
                                             " -s_srs '", wgs_crs, "'",
                                             " -t_srs '", new_crs, "' ",
                                             infile[x], " ", outfile[x])),
                   mc.cores = mc.cores, mc.preschedule = TRUE)

## Checks 
n=103
raster(infile[n])
raster(outfile[n])
raster(file.path(output_dir, sprintf("globalmask_%sk_ee.tif", proj.res.km)))
  ## Check extent with maskfile


## >> step five_Mask layers ####
## https://gitlab.unimelb.edu.au/garberj/gdalutilsaddons/-/blob/master/gdal_calc.R
infile <- outfile
outfile <- sub(".tif", "_eqar.tif", outfile)
mask_file <- file.path(output_dir, sprintf("globalmask_%sk_ee.tif", proj.res.km))

gdalUtils::gdalinfo(mask_file)
parallel::mclapply(seq_along(infile),
                   function(x) system(paste0("gdal_calc.py -A ", infile[x], " -B ", mask_file,
                                             " --calc='((B==1)*A)+(-9999*(B!=1))' --NoDataValue=-9999",
                                             " --outfile=", outfile[x])),
                   mc.cores = mc.cores, mc.preschedule = TRUE)

## Checks
gdalUtils::gdalinfo(outfile[n])
raster(infile[n])
raster(outfile[n])
plot(raster(outfile[n]))

file.remove(infile)



## ------------------------------------------------------------------------- ##
## Prepare Worldlim GCM data quartile layers - TO FIX. SKIP. USE CHELSA. ####
## ------------------------------------------------------------------------- ##
# ## ALL GCMs
# ## Source: from regSSP scripts
# ## *** replace with gdal functions...
# rcps <- c("45", "60", "85")
# models <- c("BC", "CC", "GS", "HD", "HE", "IP", "MI", "MR", "MC", "MG", "NO")
# 
# mask_file <- file.path(data_dir, "global_mask_ee.tif")
# gcm_files <- list.files(file.path(data_dir, 'gcm_30s'), full.names = T, recursive = T)
# gcm_masked_path <- file.path(processed_dir, 'gcm_masked')
# if(!dir.exists(gcm_masked_path)){dir.create(gcm_masked_path)}
# 
# 
# ## Mask data & stack by GCM
# global_mask <- raster(mask_file)
# for(model_i in models){
#   mod_stack <- list()
#   file_mod <- files_gcm[grepl(model_i, files_gcm)]
#   for(j in 1:length(rcps)){
#     file_mod_rcp <- file_mod[grepl(rcps[j], file_mod)]
#     rcp_stack <- list()
#     for(f in 1:length(file_mod_rcp)){
#       rcp_stack[[f]] <- mask(crop(raster(file_mod_rcp[f]), global_mask), global_mask)
#     }
#     mod_stack[[j]] <- readAll(brick(rcp_stack))
#   }
#   saveRDS(mod_stack, file = paste0(gcm_reg_path, "/", model_i, ".rds"))
# }
# rm(rcp_stack, mod_stack)
# 
# ## Extract cell-wise quartiles across GCM
# quartiles <- c("q1", "q2", "q3")
# gcm_files <- list.files(gcm_masked_path, full.names = T)
# gcm_quant_path <- file.path(processed_dir, 'gcm_quant')
# if(!dir.exists(gcm_quant_path)){dir.create(gcm_quant_path)}
# 
# inds <- which(!is.na(global_mask[]))
# 
# for(j in 1:length(rcps)){
#   saveRDS(stack(), file = paste0(gcm_quant_path, "/bio", "q1_", rcps[j], ".rds"))
#   saveRDS(stack(), file = paste0(gcm_quant_path, "/bio", "q2_", rcps[j], ".rds"))
#   saveRDS(stack(), file = paste0(gcm_quant_path, "/bio", "q3_", rcps[j], ".rds"))
#   print(paste0("processing rcp", rcps[j]))
#   for(k in 1:19){
#     print(paste0("processing bioclim var: ", k))
#     bio <- stack()
#     for(i in 1:length(models)){
#       print(paste0("processing model: ", i))
#       dat <- readRDS(gcm_files[[i]])[[j]]
#       bio <- stack(bio, dat[[k]])
#     }
#     
#     print(paste0("getting quartiles..."))
#     df1 <- na.omit(as.matrix(getValues(bio)))
#     c <-rowQuantiles(df1, probs = c(0.25, 0.5, 0.75))
#     for(m in 1:3){
#       bioclim <- readRDS(file = paste0(gcm_quant_path, "/bio", quartiles[m], "_", rcps[j], ".rds"))
#       r[inds] <- c[,m]
#       names(r) <- paste0("bio", k)
#       saveRDS(readAll(stack(bioclim, r)), file = paste0(gcm_quant_path, "/bio", quartiles[m], "_", rcps[j], ".rds"))
#     }
#   }
# }
# unlink(gcm_masked_path, recursive=T)




## ------------------------------------------------------------------------- ##
## Prepare distance layers ####
## ------------------------------------------------------------------------- ##
dst_folder <- file.path(data_dir, "distance_layers")
if(!dir.exists(dst_folder)) {
  dir.create(dst_folder)
}

## >> Builtup ####
## Distance is calculated only from most densely populated areas, i.e. -values 100 
infile <- file.path(data_dir, "builtup", "GHS_BUILT_LDS2014_GLOBE_R2018A_54009_1K_V2_0.tif")
gdalinfo(infile)

## Convert Projection
t_crs <- "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"
outfile <- gsub(".tif$", "_wgs.tif", infile)
system(sprintf("gdalwarp -overwrite -ot Byte -r bilinear -s_srs '%s' -t_srs '%s' %s %s",
               t_crs, wgs_crs, infile, outfile))
gdalinfo(outfile)
## Calculate distance layer
infile <- outfile
outfile <- file.path(dst_folder, "dst_builtup.tif")

system(sprintf("gdal_proximity.py %s %s -values %d -of GTiff -distunits GEO -nodata %d -co compress=LZW -co BIGTIFF=YES", infile, outfile, 100, -9999))
gdalinfo(outfile)


## >> Roads ####
## >> >> The input file geodatabase
## https://gis.stackexchange.com/questions/151613/reading-feature-class-in-file-geodatabase-using-r
roadsgdb <- file.path(data_dir, "groads/groads-v1-global-gdb/gROADS_v1.gdb")

## >> >> List all feature classes in a file geodatabase
subset(ogrDrivers(), grepl("GDB", name))
ogrListLayers(roadsgdb)

## >> >> Read the feature class
roads <- readOGR(dsn=roadsgdb,layer="Global_Roads")
summary(roads)
unique(roads$FCLASS)
unique(roads$EXS)

## >> >> Subset the road attributes (see gROADSv1-ReadMe.txt) 
## 1. FClass = Functional class (1=Highway, 2=Primary, 3=Secondary, 4=Tertiary, 5=Local/ Urban, 6=Trail, 7=Private, 0=Unspecified)
## 2. EXS = Existence Category (1=Definite, 2=Doubtful, 0=Unspecified)
roadsub <- roads[(roads$FCLASS == (0:5) & roads$EXS == c(0,1)),]
unique(roadsub$FCLASS)
unique(roadsub$EXS)

## >> Drop columns
roadsub <- roadsub[, 4]

## >> >> Write shapefile
writeOGR(roadsub, dsn = file.path(data_dir, "groads"), layer = "groads", driver="ESRI Shapefile")

## >> >> Rasterise @ .0025 ~ 250m
system(paste0("gdal_rasterize -at -burn 1 -ot Byte -tr .0025 .0025 ",
              file.path(data_dir, "groads", "groads.shp "),
              file.path(data_dir, "groads", "groads.tif")))

rm(roads, roadsub, roadsgdb)

## >> >> Caluclate distance raster
infile <- file.path(data_dir, "groads", "groads.tif")
gdalinfo(infile)

outfile <- file.path(dst_folder, "dst_roads.tif")

system(sprintf("gdal_proximity.py %s %s -values %d -of GTiff -distunits GEO -nodata %d -co compress=LZW -co BIGTIFF=YES", infile, outfile, 1, -9999))
gdalinfo(outfile)


## >> Population ####
## Distance is calculated from all non 0 values
infile <- file.path(data_dir, "population/SEDAC/gpw-v4-population-density-rev11_2000_30_sec_tif", "gpw_v4_population_density_rev11_2000_30_sec.tif")
gdalinfo(infile, stats = TRUE)

outfile <- file.path(dst_folder, "dst_pop.tif")

system(sprintf("gdal_proximity.py %s %s -of GTiff -distunits GEO -nodata %s -co compress=LZW -co BIGTIFF=YES", infile, outfile, -9999))
gdalinfo(outfile)



## >> Water bodies ####
## Source: https://www.worldwildlife.org/publications/global-lakes-and-wetlands-database-lakes-and-wetlands-grid-level-3
infile <- file.path(data_dir, "lakesrivers/GLWD/glwd_3")
gdalinfo(infile)


## >> >> Lakes & rivers ####
## >> >> Extract rivers and lakes from data
outfile <- file.path(data_dir, "lakesrivers/GLWD", "lakesrivers.tif")
system(paste0("gdal_calc.py -A ", infile,
              " --calc='((A==1) + (A==2) + (A==3) + (A==4) + (A==5))*1 + ( + (A==6) + (A==7) + (A==8) + (A==9) + (A==10) + (A==11) + (A==12))*0' --NoDataValue=-9999",
              " --outfile=", outfile))
gdalUtils::gdalinfo(outfile)

## >> >> Caluclate distance raster
infile <- outfile
outfile <- file.path(dst_folder, "dst_lakesrivers.tif")
system(sprintf("gdal_proximity.py %s %s -values %d -of GTiff -distunits GEO -nodata %d -co compress=LZW -co BIGTIFF=YES", infile, outfile, 1, -9999))
gdalinfo(outfile)

## >> >> Specify projection
infile <- outfile
outfile <- gsub(".tif$", "_wgs.tif", infile)
wgs_crs = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0" 
system(sprintf("gdalwarp -overwrite -ot Byte -t_srs '%s' %s %s",
               wgs_crs, infile, outfile))


## >> >> Wetlands ####
## >> >> Extract wetlands from data
infile <- file.path(data_dir, "lakesrivers/GLWD/glwd_3")
outfile <- file.path(data_dir, "lakesrivers/GLWD", "wetlands.tif")
system(paste0("gdal_calc.py -A ", infile,
              " --calc='((A==1) + (A==2) + (A==3) + (A==4) + (A==5))*0 + ( + (A==6) + (A==7) + (A==8) + (A==9) + (A==10) + (A==11) + (A==12))*1' --NoDataValue=-9999",
              " --outfile=", outfile))
gdalUtils::gdalinfo(outfile)

## >> >> Caluclate distance raster
infile <- outfile
outfile <- file.path(dst_folder, "dst_wetlands.tif")
system(sprintf("gdal_proximity.py %s %s -values %d -of GTiff -distunits GEO -nodata %d -co compress=LZW -co BIGTIFF=YES", infile, outfile, 1, -9999))
gdalinfo(outfile)

## >> >> Specify projection
infile <- outfile
outfile <- gsub(".tif$", "_wgs.tif", infile)
wgs_crs = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0" 
system(sprintf("gdalwarp -overwrite -ot Byte -t_srs '%s' %s %s",
               wgs_crs, infile, outfile))


## >> Protected areas ####
## Distance is calculated from all non 0 values (only 1 and 0 values in this layer)
infile <- file.path(data_dir, "protectedareas/wdpa.tif")
gdalinfo(infile)

outfile <- file.path(dst_folder, "dst_wdpa.tif")

system(sprintf("gdal_proximity.py %s %s -values %d -of GTiff -distunits GEO -nodata %s -co compress=LZW -co BIGTIFF=YES", infile, outfile, 1, -9999))
gdalinfo(outfile)




## ------------------------------------------------------------------------- ##
## Process environmental layers - Equal Area projection ####
## ------------------------------------------------------------------------- ##

## >> Load filenames & check/specify CRS for each ####
## SRTM
srtm <- file.path(data_dir, "srtm/mn30_grd/srtm.adf")
gdalinfo(srtm)
gdalsrsinfo(srtm)

## Soil
soil <- list.files(file.path(data_dir,"soil"), pattern = "*.dat", 
                   full.names = TRUE, recursive = TRUE)
lapply(soil, gdalsrsinfo, as.CRS = TRUE)
gdalinfo(soil[1])

## >> >> Specify projection

wgs_crs = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0" 
mclapply(seq_along(soil),
         function(x) system(sprintf("gdalwarp -overwrite -ot Byte -t_srs '%s' %s %s",
                                    wgs_crs, soil[x], gsub(".dat$", "_wgs.tif", soil)[x])),
         mc.cores = length(soil), mc.preschedule = TRUE)
## OR
# mclapply(seq_along(soil),
#          function(x) system(sprintf("gdal_translate -ot Byte -a_srs '%s' %s %s",
#                                     wgs_crs, soil[x], gsub(".dat$", "_wgs2.tif", soil)[x])),
#          mc.cores = length(soil), mc.preschedule = TRUE)
soil <- list.files(file.path(data_dir,"soil"), pattern = "_wgs.tif", 
                   full.names = TRUE, recursive = TRUE)
lapply(soil, gdalsrsinfo)

  # ## WorldClim: current data
  # bio_current <- list.files(paste0(data_dir, "/bio_30s"), pattern = "bio_current*", full.names = TRUE)
  # 
  # 
  # ## Worldclim: future GCM data - One GCM: BCC-CSM1-1
  # bio_rcp45 <- list.files(file.path(data_dir, "gcm_30s"), pattern = "*bc45*", full.names = TRUE)
  # bio_rcp85 <- list.files(file.path(data_dir, "gcm_30s"), pattern = "*bc85*", full.names = TRUE)
  # # bio_rcp60 <- list.files(file.path(data_dir, "gcm_30s"), pattern = "*bc60*", full.names = TRUE)


## CHELSA current
biocurrent <- list.files(file.path(data_dir, "CHELSA/envicloud/chelsa/chelsa_V2/GLOBAL/climatologies/1981-2010/bio"),
                         full.names = TRUE, recursive = TRUE)
gdalinfo(biocurrent[1])
lapply(biocurrent, gdalsrsinfo)

## >> >> Rename files
file.rename(biocurrent, gsub("_V.2.1.tif$", ".tif", biocurrent))
biocurrent <- list.files(file.path(data_dir, "CHELSA/envicloud/chelsa/chelsa_V2/GLOBAL/climatologies/1981-2010/bio"),
                         full.names = TRUE, recursive = TRUE)

kgvars <- grep("CHELSA_kg", biocurrent, value = TRUE)
biocurrent <- grep("CHELSA_bio", biocurrent, value = TRUE)


## CHELSA future
# biofuture <- list.files("/home/payalb/gsdms_r_vol/tempdata/research-cifs/proj-2200_nature_futures21-1128.4.411/global/cmip6/bio/mean", full.names = TRUE, recursive = TRUE)
biofuture <- list.files(file.path(data_dir, "CHELSA/future/mean"), full.names = TRUE, recursive = TRUE)
gdalinfo(biofuture[1])

## >> >> Rename files
# file.rename(biofuture, gsub("_V.2.1.tif$", ".tif", biofuture))


## Hansen treecover
tree <- list.files(file.path(data_dir, "Hansen_treecover"), 
                   pattern = ".tif$", full.names = TRUE)
gdalinfo(tree)
gdalsrsinfo(tree)


## Distance layers
dst_layers <- list.files(file.path(data_dir, "distance_layers"),
                         full.names = TRUE)
dst_layers <- dst_layers[-c(3,8)]
gdalinfo(dst_layers[4])
lapply(dst_layers, gdalsrsinfo)



## >> Combine filenames ####
cov_files <- c(biocurrent, biofuture, 
               kgvars, srtm, soil, tree, dst_layers)
all(lapply(cov_files, file.exists))



## >> step one_Clip by e ####
infiles = cov_files
outfiles = file.path(output_dir, paste0(tools::file_path_sans_ext(basename(infiles)), "_clip.tif"))
# outfiles = gsub("\\..*", ".tif", outfiles)
new_extent <- "-180 -60 180 90"

parallel::mclapply(seq_along(infiles),
                   function(x) system(paste0("gdalwarp -overwrite -ot Byte",
                                             " -te ", new_extent, " ",
                                             infiles[x], " ", outfiles[x])),
                   mc.cores = mc.cores, mc.preschedule = TRUE)

length(list.files(output_dir, pattern = "_clip"))
length(infiles)

# ## ...
# outfiles = outfiles[!tools::file_path_sans_ext(basename(infiles)) %in% gsub("_clip", "", tools::file_path_sans_ext(basename(list.files(output_dir, pattern = "_clip"))))]
# infiles = infiles[!tools::file_path_sans_ext(basename(infiles)) %in% gsub("_clip", "", tools::file_path_sans_ext(basename(list.files(output_dir, pattern = "_clip"))))]



## >> step two_Reproject to Equal Earth ####
## resamplig method: bilinear  (https://support.esri.com/en/technical-article/000005606)
x <- outfiles
infiles <- outfiles[26:37]
lapply(infiles, gdalsrsinfo)

outfiles <- gsub("_clip", "_temp", infiles)
parallel::mclapply(seq_along(infiles),
                   function(x) system(sprintf("gdal_translate -ot Byte -a_srs 'EPSG:4326' %s %s",
                                              infiles[x], outfiles[x])),
                   mc.cores = mc.cores, mc.preschedule = TRUE)
  ## this step is nneeded to resolve the error: ERROR 1: Too many points (xx out of xx) failed to transform, unable to compute output bounds. Warning 1: Unable to compute source region for output window x,x,x,x, skipping.


infiles <- outfiles
outfiles <- gsub("_clip", "_ee", infiles)
new_res <- c(proj.res.km*1000, proj.res.km*1000)
wgs_crs = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0" 
new_crs = '+proj=eqearth +ellips=WGS84 +wktext'
parallel::mclapply(seq_along(infiles),
                   function(x) system(sprintf("gdalwarp -overwrite -ot Byte -r bilinear -tr %s -s_srs '%s' -t_srs '%s' %s %s",
                                      paste(new_res, collapse = " "), 
                                      wgs_crs, new_crs, 
                                      infiles[x], outfiles[x])),
                   mc.cores = mc.cores, mc.preschedule = TRUE)
file.remove(infiles)

list.files(output_dir, pattern = "_ee")
lapply(outfiles, gdalsrsinfo)
# system(sprintf("gdalwarp -overwrite -ot Byte -r bilinear -tr %s -s_srs '%s' -t_srs '%s' srcfile %s dstfile %s",
#        paste(new_res, collapse = " "), wgs_crs, new_crs, infiles, outfiles)

# system(sprintf("gdalwarp -overwrite -ot Byte -te %s -r bilinear -tr %s -s_srs 'EPSG:4326' -t_srs '%s' srcfile %s dstfile %s", 
#                new_extent, new_res, new_crs, infile, outfile))

  ## Checks 
  n=62
  raster(infiles[n])
  raster(outfiles[n])
  
  file.remove(infiles)


## >> step three_Clip to make #cells equal between mask and layer ####
infiles <- outfiles
outfiles <- sub("_ee", "_clip2", outfiles)
mask_file <- file.path(output_dir, sprintf("globalmask_%sk_ee.tif", proj.res.km))
new_extent <- extent(raster(mask_file))

parallel::mclapply(seq_along(infiles),
                   function(x) system(paste0("gdalwarp -overwrite -ot Byte -te ",
                                             paste(new_extent[1], new_extent[3],
                                                   new_extent[2], new_extent[4]), " ",
                                             infiles[x], " ", outfiles[x])),
                   mc.cores = mc.cores, mc.preschedule = TRUE)

## Checks
raster(infiles[n])
raster(outfiles[n])
raster(mask_file)

file.remove(infiles)


## >> step four_Mask layers ####
## https://gitlab.unimelb.edu.au/garberj/gdalutilsaddons/-/blob/master/gdal_calc.R
infiles <- outfiles
outfiles <- sub("_clip2", "_eqar", outfiles)
mask_file <- file.path(output_dir, sprintf("globalmask_%sk_ee.tif", proj.res.km))

gdalUtils::gdalinfo(mask_file)
parallel::mclapply(seq_along(infiles),
                   function(x) system(paste0("gdal_calc.py -A ", infiles[x], " -B ", mask_file,
                                             " --calc='((B==1)*A)+(-9999*(B!=1))' --NoDataValue=-9999",
                                             " --outfile=", outfiles[x])),
                   mc.cores = mc.cores, mc.preschedule = TRUE)
gdalUtils::gdalinfo(outfiles[n])

## Checks
raster(infiles[n])
raster(outfiles[n])
unique(values(raster(outfiles[n])))
plot(raster(outfiles[n]))

file.remove(infiles)




## ------------------------------------------------------------------------- ##
## Find min non-NA set values across mask and covariates and sync NAs ####
## ------------------------------------------------------------------------- ##
rm(list=setdiff(ls(), c("output_dir", "data_dir", "proj.res.km")))
source("/home/payalb/gsdms_r_vol/tempdata/workdir/landuse_projects/land-use/gdal_calc.R") 

## Get input covariate filenames
infiles <- list.files(output_dir, full.names = TRUE, pattern = ".tif")
basename(infiles)
file.remove(grep(".xml", infiles, value = TRUE))
infiles <- list.files(output_dir, full.names = TRUE, pattern = ".tif")

length(infiles) - length(grep("mask|ESA", infiles, value = TRUE)) == 
  length(infiles[!infiles %in% grep("mask|ESA", infiles, value = TRUE)])
  ## Check

infiles <- infiles[!infiles %in% grep("mask|ESA", infiles, value = TRUE)]
basename(infiles)

mask_file <- file.path(output_dir, 
                       sprintf("globalmask_%sk_ee.tif", proj.res.km))
mask_file2 <- file.path(output_dir, 
                        sprintf("globalmask_%sk_ee_minNA.tif", proj.res.km))
file.remove(mask_file2)
file.copy(mask_file, mask_file2)

## Update mask based on covariate layers for syncing NAs
for(j in 1:length(infiles)){
  
  print(paste0("Updating mask using file = ", j, " of ", length(infiles) 
               ,":", infiles[j]))
  
  # gdalcalc(calc="((A>-9999) & (B==1))", infile = list(A=infiles[j], B=mask_file2),outfile = mask_file2,
  #          NoDataValue=-9999, overwrite=TRUE)
  
  system(paste0("gdal_calc.py -A ", infiles[j], " -B ", mask_file2,
                " --calc='((A>-9999) & (B==1))' --NoDataValue=-9999",
                " --outfile=", mask_file2))
}

## Update covariate layers based on updated mask for syncing NAs
for(j in 1:length(infiles)){
  
  print(paste0("Masking file = ", j, " of ", length(infiles) ,":", 
               infiles[j], " using updated mask"))
  
  gdalmask(infile = infiles[j], mask = mask_file2, outfile = infiles[j], 
           output_Raster = FALSE, overwrite=TRUE, verbose=TRUE)
  
}

## Checks
freq(raster(mask_file))
freq(raster(mask_file2))
## this might remove (if applicable) some ones from the original mask and made them zeros. From looking at the map, it seems like it is mostly picking up on large lakes? Looks like the great lakes, Lake victoria, and Lake Baikal are some of the removed cells.
gdalinfo(mask_file, stats = TRUE)
gdalinfo(mask_file2, stats = TRUE)

lapply(infiles, gdalinfo, stats = TRUE)
## check size for all
lapply(infiles, gdalsrsinfo, as.CRS = TRUE)

n = 103
gdalinfo(infiles[n])
raster(infiles[n])
unique(values(raster(infiles[n])))
plot(raster(infiles[n]), main = tools::file_path_sans_ext(basename(infiles[n])))




## ------------------------------------------------------------------------- ##
## Create seperate tifs for SRTM variables  ####
## ------------------------------------------------------------------------- ##

infile <- list.files(output_dir, pattern = "srtm_eqar.tif$", full.names = TRUE)
gdalinfo(infile)

outfile <- gsub("srtm", "slope", infile)
system(paste0("gdaldem slope ", infile, " ", outfile))

outfile <- gsub("srtm", "aspect", infile)
system(paste0("gdaldem aspect ", infile, " ", outfile))

outfile <- gsub("srtm", "terrain", infile)
system(paste0("gdaldem TRI ", infile, " ", outfile))

outfile <- gsub("srtm", "roughness", infile)
system(paste0("gdaldem roughness ", infile, " ", outfile))

file.rename(infile, gsub("srtm", "elevation", infile))



## ------------------------------------------------------------------------- ##
## Prepare CHELSA future layers - Interpolate between y_0 and y_end ####
## ------------------------------------------------------------------------- ##

biofuture <- list.files(file.path(data_dir, "CHELSA/future/mean"),
                        full.names = TRUE, recursive = TRUE)
biofuture_ssp1 <- grep("ssp1", biofuture, value = TRUE) 
## 19 vars * 3 time steps = 57 files
biofuture_ssp3 <- grep("ssp3", biofuture, value = TRUE) 
biofuture_ssp5 <- grep("ssp5", biofuture, value = TRUE) 




## --------------------------------------------------------- ##
## Create covariate data.tables for model and prediction ####
## -------------------------------------------------------- ##

infile <- list.files(output_dir, pattern = "_eqar(.*).tif$", full.names = TRUE)
## matching everything with _eqar and ending in .tif but also anythign in between the two using (.*)
infile <- infile[-grep('ESA_landuse_reclass_eqar.tif$', infile)]
## remove ESA_landuse_reclass_eqar.tif
lapply(infile, gdalinfo, stats = TRUE)
## check size for all


## >> Generate quadrature (background) points - TO MODIFY ####
## ***** Need to swap out rasterToPoints() and stack() when working with large rasters *****
## NOTE: Trial alternative approach using biogeoregions

mask_file <- file.path(input_dir, sprintf("globalmask_%sk_ee_minNA.tif", proj.res.km))
global_mask0 <- raster(mask_file); freq(global_mask0)
global_mask0[which(is.na(global_mask0[]))] <- 0; freq(global_mask0)
rpts <- rasterToPoints(global_mask0, spatial=TRUE)
mask_pts <- data.frame(rpts@data, X=coordinates(rpts)[,1], Y=coordinates(rpts)[,2])
mask_pts <- mask_pts[,-1]
fwrite(mask_pts, file = file.path(output_dir, sprintf("globalmask_%sk_eeXY.tif", proj.res.km)))
# mask_pts <- fread(file.path(output_dir, sprintf("globalmask_%sk_eeXY.tif", proj.res.km)))
dat <- stack(infile[-grep('bc45|bc85', infile)])
dat <- as.data.table(na.omit(cbind(mask_pts, as.matrix(dat))))


## >> Test for correlations in covariates ####
preds <- colnames(correlations(dat, thresh = 0.7, N = 200000))
preds <- sub("_current_", "", preds)
write.csv(preds, file.path(output_dir, "uncorrelated_covs.csv"), 
          row.names = FALSE)
# preds <- fread(file.path(output_dir, "uncorrelated_covs.csv"))$x

## >> Prepare background point data for model + remove correlated covariates ####
lucovs <- grep("landuse", names(dat), value = TRUE)
names(dat) <- sub("_current_", "", names(dat))
all(names(dat[,names(dat)[names(dat) %in% c(preds, lucovs, "X", "Y")], with = FALSE]) == preds)

dat <- dat[,names(dat)[names(dat) %in% c(preds, lucovs, "X", "Y")], with = FALSE] ## to ensure landuse, X, Y are retained
dim(dat)
names(dat) <- sub("_eqar", "", names(dat))
names(dat) <- sub("ESA_landuse_reclass_", "", names(dat))

## Rearrange covariate data accroding to static, dynamic and land use vars
## NOTE: Can be removed once we've figured out variable selection for models
static <- names(dat)[-grep("lu|bio", names(dat))]
dynamic <- names(dat)[grep("bio", names(dat))]
landuse <- names(dat)[grep("lu", names(dat))]
arrange_cols <- c(static, dynamic, landuse)
rm(static, dynamic, landuse)

dat <- dat[,..arrange_cols]
names(dat)
fwrite(dat, file = file.path(output_dir, "covariates_model.csv"))

dat_mod <- data.table::copy(dat)
rm(dat)

## >> Prepare prediction points data + remove correlated covariates ####
## RCP 45
dat <- stack(infile[-grep('bio_current|bc85', infile)])
dat <- as.data.table(na.omit(cbind(mask_pts, as.matrix(dat))))
names(dat) <- sub("bc45bi70", "bio", names(dat))
dat <- dat[,names(dat)[names(dat) %in% c(preds, lucovs, "X", "Y")], with = FALSE] ## retain landuse, X, Y 
dim(dat)
names(dat) <- sub("_eqar", "", names(dat))
names(dat) <- sub("ESA_landuse_reclass_", "", names(dat))
dat <- dat[,..arrange_cols]
all(names(dat) == names(dat_mod))

fwrite(dat, file = file.path(output_dir, "covariates_predict_rcp45.csv"))

dat_rcp45 <- data.table::copy(dat)
rm(dat)

## RCP 85
dat <- stack(infile[-grep('bio_current|bc45', infile)])
dat <- as.data.table(na.omit(cbind(mask_pts, as.matrix(dat))))
names(dat) <- sub("bc85bi70", "bio", names(dat))
dat <- dat[,names(dat)[names(dat) %in% c(preds, lucovs, "X", "Y")], with = FALSE] ## retain landuse, X, Y 
dim(dat)
names(dat) <- sub("_eqar", "", names(dat))
names(dat) <- sub("ESA_landuse_reclass_", "", names(dat))
dat <- dat[,..arrange_cols]
all(names(dat) == names(dat_mod))

fwrite(dat, file = file.path(output_dir, "covariates_predict_rcp85.csv"))

dat_rcp85 <- data.table::copy(dat)
rm(dat)

## Checks
dim(dat_mod); dim(dat_rcp45); dim(dat_rcp85)







## ------------------------------------------------------------------------- ##
## EXTRAS ---------------- ## ####

# ## >> >> step two_change NoData values to -9999 ####
change_values = "/home/payalb/gsdms_r_vol/tempdata/workdir/gsdms/scripts/change_values.sh"

luvars <- c(luvars[c(1:5, 11:12)], "crop")

hurttfiles <- list.files(hurtt_out, 
                         pattern = paste0(luvars, collapse = "|"), 
                         full.names = TRUE)
length(hurttfiles)
length(luvars)*12*3 # #variables * #time steps * #scenarios

gdalinfo(hurttfiles[1])
raster((grep("crop", hurttfiles, value = TRUE)[1]))

parallel::mclapply(seq_along(hurttfiles),
                   function(x) system(paste0("bash ", 
                                             change_values, 
                                             " ", x, " -9999")),
                   mc.cores = mc.cores, mc.preschedule = TRUE)

gdalinfo(hurttfiles[1])
system(paste0("bash ", change_values, " ", hurttfiles[1], " -9999"))
gdalUtils::gdalinfo(gsub(".tif", "_edt.tif", hurttfiles[1]))

file.copy(gsub(".tif", "_edt.tif", hurttfiles[1]), file.path(output_dir, "ESA_landuse_reclass.tif"))
file.remove(c(gsub(".tif", "_edt.tif", outfile), gsub(".tif", "_temp.tif", outfile)))


## ------------------------------------------------------------------------- ##
## Process CATEGORICAL covariate layers (land use) - Equal Area projection ####
## ------------------------------------------------------------------------- ## 

# ## OLD LAND USE DATA 
# ## >> Source: ESA (fractional land use) ####
# ## http://maps.elie.ucl.ac.be/CCI/viewer/download.php
# landuse <- list.files("/home/payalb/gsdms_r_vol/tempdata/research-cifs/uom_data/gsdms_data/landuse/ESA_CCI_Landcover",
#                       pattern = "landuse.tif$", full.names = TRUE)
# outfile <- gsub("ESA_landuse", "ESA_landuse_reclass", landuse)
# 
# ## >> >> step one_reclassify ####
# ## Recalssify according to flutes classes - See landuse_classifications.xlsx
# gdalUtils::gdalinfo(landuse)
# system(paste0("gdal_calc.py -A ", landuse,
#               " --calc='((A==10) + (A==11) + (A==12) + (A==20) + (A==30))*1 + ((A==40))*2 + ((A==50) + (A==60) + (A==61) + (A==62) + (A==70) + (A==71) + (A==72) + (A==80) + (A==81) + (A==82) + (A==90) + (A==100) + (A==160) + (A==170))*3 + ((A==110) + (A==130))*4 + ((A==120) + (A==121) + (A==122))*5 + ((A==180))*6 + ((A==190))*7 + ((A==140) + (A==150) + (A==151) + (A==152) + (A==153) + (A==200) + (A==201) + (A==202) + (A==220))*8 + ((A==210))*9' --NoDataValue=0",
#               " --outfile=", outfile))
# gdalUtils::gdalinfo(outfile)
# 
# 
# ## >> >> step two_change NoData values to -9999 ####
# change_values = "/home/payalb/gsdms_r_vol/tempdata/workdir/gsdms/scripts/change_values.sh"
# system(paste0("bash ", change_values, " ", outfile, " -9999"))
# gdalUtils::gdalinfo(gsub(".tif", "_edt.tif", outfile))
# file.copy(gsub(".tif", "_edt.tif", outfile), file.path(output_dir, "ESA_landuse_reclass.tif"))
# file.remove(c(gsub(".tif", "_edt.tif", outfile), gsub(".tif", "_temp.tif", outfile)))
# 
# 
# ## >> >> step three_create separate tif for each landuse class ####
# ## NOTES: get unique values: https://gis.stackexchange.com/questions/33388/python-gdal-get-unique-values-in-discrete-valued-raster
# ## NOTES: run python code: https://rstudio.github.io/reticulate/
# infile <- file.path(output_dir, "ESA_landuse_reclass.tif")
# vals <- 1:9
# class <- c("crop", "crop_mosaic", "forest", "grass",
#            "shrub", "wetland", "urban", "other", "water")
# outfile <- file.path(output_dir, paste0("ESA_lu", vals, class, ".tif"))
# 
# parallel::mclapply(seq_along(vals),
#                    function(x) system(paste0("gdal_calc.py -A ", infile,
#                                              " --calc='((A==", x, "))*1 + (",
#                                              paste0("(A==", vals[-x], ")",
#                                                     collapse = " + "),
#                                              ")*0' --NoDataValue=-9999",
#                                              " --outfile=", outfile[x])),
#                    mc.cores = length(vals), mc.preschedule = TRUE)
# gdalUtils::gdalinfo(outfile[1])
# 
# 
# ## >> >> step four_clip by e ####
# infile <- outfile
# outfile <- gsub(".tif", "_clip.tif", infile)
# outfile = gsub("\\..*", ".tif", outfile)
# new_extent <- "-180 -60 180 90"
# 
# parallel::mclapply(seq_along(vals),
#                    function(x) system(paste0("gdalwarp -overwrite -ot Byte",
#                                              " -te ", new_extent, " ",
#                                              infile[x], " ", outfile[x])),
#                    mc.cores = length(vals), mc.preschedule = TRUE)
# 
# # system(paste0("gdalwarp -overwrite -ot Byte",
# #               " -te ", new_extent, " ",
# #               infile, " ", outfile))
# gdalUtils::gdalinfo(outfile[1])
# 
# 
# 
#   # ## >> >> step five_create fractional layers at 0.1 degrees (= ~10k) -- NOT WORKING -- ####
#   # infile <- outfile
#   # outfile <- gsub("_clip.tif", "_frac.tif", infile)
#   # new_res <- c(0.1, 0.1)
#   # 
#   # system(paste0("gdalwarp -overwrite -r average -tr ",
#   #               paste(new_res, collapse = " "), " ",
#   #               infile[1], " ", outfile[1]))
# 
# 
# ## >> >> step six_reproject to Equal Earth ####
# ## resamplig method: near (https://support.esri.com/en/technical-article/000005606)
# ## see: https://gis.stackexchange.com/questions/352476/bilinear-resampling-with-gdal-leaves-holes
# infile <- outfile #list.files(output_dir, pattern = "_clip", full.names = TRUE)
# outfile <- gsub("_clip.tif", "_frac.tif", infile)
# new_res <- c(10000,10000)
# new_crs = '+proj=eqearth +ellips=WGS84 +wktext'
# wgs_crs = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
# 
#   # parallel::mclapply(seq_along(vals),
#   #                    function(x) system(paste0("gdalwarp -overwrite -ot Byte -r average -tr ",
#   #                                              paste(new_res, collapse = " "),
#   #                                              " -s_srs '", wgs_crs, "'",
#   #                                              " -t_srs '", new_crs, "' ",
#   #                                              infile[x], " ", outfile[x])),
#   #                    mc.cores = length(vals), mc.preschedule = TRUE)
# 
# system(paste0("gdalwarp -overwrite -r average -tr ",
#               paste(new_res, collapse = " "),
#               " -t_srs '", new_crs, "' ",
#               infile[1], " ", outfile[1]))
# 
# gdalUtils::gdalinfo(infile[1])
# sort(unique(values(raster(outfile[1]))))
# range(unique(values(raster(outfile[1]))), na.rm = TRUE)
# file.remove(infile)
# 
#   # ## >> >> step_clip: to make #cells equal between mask and layer ####
#   # infile <- outfile
#   # outfile <- sub("_ee", "_clip2", outfile)
#   # mask_file <- file.path(output_dir, sprintf("globalmask_%sk_ee.tif", proj.res.km))
#   # new_extent <- extent(raster(mask_file))
#   #
#   # system(paste0("gdalwarp -overwrite -ot Byte -te ",
#   #               paste(new_extent[1], new_extent[3],
#   #                     new_extent[2], new_extent[4]), " ",
#   #               infile, " ", outfile))
#   # file.remove(infile)
# 
# 
# ## >> >> step seven_mask ####
# ## https://gitlab.unimelb.edu.au/garberj/gdalutilsaddons/-/blob/master/gdal_calc.R
# infile <- outfile
# outfile <- sub("_ee", "_eqar", outfile)
# mask_file <- file.path(output_dir, sprintf("globalmask_%sk_ee.tif", proj.res.km))
# 
# gdalUtils::gdalinfo(mask_file)
# parallel::mclapply(seq_along(vals),
#                    function(x) system(paste0("gdal_calc.py -A ", infile[x],
#                                              " -B ", mask_file,
#                                              " --calc='((B==1)*A) + (-9999*(B!=1))' --NoDataValue=-9999",
#                                              " --outfile=", outfile[x])),
#                    mc.cores = length(vals), mc.preschedule = TRUE)
# gdalUtils::gdalinfo(outfile[3])
# file.remove(infile)




# ## >> Source: Copernicus (fractional land use) ####
# ## Link @ GEE: https://developers.google.com/earth-engine/datasets/catalog/COPERNICUS_Landcover_100m_Proba-V_Global
# ## Direct download link: https://zenodo.org/communities/copernicus-land-cover/search?page=1&size=20
# ## Data processing for global layer @ GEE co-authored by MC Loor: https://code.earthengine.google.com/?scriptPath=users%2Fpayalbal%2Fglobal_layers%3Afraclu_sklu. See processing in GEE for reclassification scheme.
# 
# ## Land use classes used:
# ## 1    urban
# ## 2    crop
# ## 3    forest
# ## 4    grass
# ## 5    other
# ## 6    NA
# 
# # ... gee.py script
# # ...download by tile
# 
# ## Stiching tiles by Maria del Mar Quiroga
# ## ref: https://stackoverflow.com/a/50235578
# 
# ## step one_Get a list of all tif tiles downloaded from Google Earth Engine
# tile_files <- list.files(path = file.path(data_dir, "landuse/Copernicus/global_frac/sklu_classes/tifs/"), pattern = "*.tif", full.names = TRUE)
# 
# ## step two_Build a virtual raster file stitching all tiles
# ## WARNING: This will fail if the file it is trying to write to (output.vrt) already exists
# gdalbuildvrt(gdalfile = tile_files, output.vrt = file.path(data_dir, "copernicus","lu_world.vrt"))
# 
# ## step three_Copy the virtual raster to an actual physical file
# ## WARNING: This takes ~5 minutes to run
# gdal_translate(src_dataset = file.path(data_dir, "copernicus", "lu_world.vrt"),
#                dst_dataset = file.path(data_dir, "copernicus", "lu_world.tif"),
#                output_Raster = FALSE,
#                options = c("BIGTIFF=YES", "COMPRESSION=LZW"))
# 
# ## step four_Save each band (i.e. land use class) as a tif file
# landuse <- file.path(file.path(data_dir, "copernicus", "lu_world.tif"))
# landuse <- brick(landuse)
# for (i in 1:nlayers(landuse)){
#   temp <- landuse[[i]]
#   writeRaster(temp, filename = file.path(data_dir, "copernicus", paste0("landuse_class", i, ".tif")))
# }
# rm(temp, landuse)
# ## landuse classes: urban, crop, forest, grass, other
# lu1 <- file.path(file.path(data_dir, "copernicus", "landuse_class1.tif"))
# lu2 <- file.path(file.path(data_dir, "copernicus", "landuse_class2.tif"))
# lu3 <- file.path(file.path(data_dir, "copernicus", "landuse_class3.tif"))
# lu4 <- file.path(file.path(data_dir, "copernicus", "landuse_class4.tif"))
# lu5 <- file.path(file.path(data_dir, "copernicus", "landuse_class5.tif"))
# 
# ## load landuse layer
# landuse <- c(lu1, lu2, lu3, lu4, lu5)



# ## Copy files from land-use projects
# lufiles <- list.files("/home/payalb/gsdms_r_vol/tempdata/workdir/landuse_projects/land-use/processed_layers/", 
#                      full.names = TRUE, all.files = TRUE)
# basename(tools::file_path_sans_ext(lufiles))
# x <- grep("landuse|lu", basename(tools::file_path_sans_ext(lufiles)), value = FALSE)
# lufiles <- lufiles[-x]
# 
# file.copy(lufiles, output_dir,
#           overwrite = FALSE, recursive = TRUE,
#           copy.mode = TRUE, copy.date = TRUE)



# ## OLD LAKES & RIVERS DATA
# ## >> Lakes ####
# ## Source: http://www.soest.hawaii.edu/wessel/gshhg/
# ## See SHAPEFILES.TXT for notes on downloaded files
# ## >> >> Read in files for resolution l, levels 2-4
# l2 <- readOGR(file.path(gsdms_data, "lakesrivers/GSHHS/GSHHS_shp/l", "GSHHS_l_L2.shp"))
# # ## same as
# # l2 <- sf::st_read(l2)
# # l2sp <- as(l2, "Spatial")
# l3 <- readOGR(file.path(gsdms_data, "lakesrivers/GSHHS/GSHHS_shp/l", "GSHHS_l_L3.shp"))
# l4 <- readOGR(file.path(gsdms_data, "lakesrivers/GSHHS/GSHHS_shp/l", "GSHHS_l_L4.shp"))
# 
# ## >> >> Combine levels 2-4 & write as shapefile
# temp <- rgeos::gUnion(l2,l3)
# temp2 <- rgeos::gUnion(temp,l4)
# writeOGR(temp2, dsn = file.path(gsdms_data, "lakesrivers/GSHHS/"), 
#          layer = "lakes_l24", driver="ESRI Shapefile")
# 
# ## >> >> Rasterise
# system(paste0("gdal_rasterize -at -burn 1 -ot Byte -tr .0025 .0025 ",
#               file.path(gsdms_data, "lakesrivers/GSHHS/", "lakes_l24.shp "),
#               file.path(gsdms_data, "lakesrivers/GSHHS/", "lakes_l24.tif")))
# 
# ## >> >> Caluclate distance raster
# infile <- file.path(gsdms_data, "lakesrivers/GSHHS/", "lakes_l24.tif")
# gdalinfo(infile)
# 
# outfile <- file.path(dst_folder, "dst_lakes.tif")
# 
# system(sprintf("gdal_proximity.py %s %s -values %d -of GTiff -distunits GEO -nodata %d -co compress=LZW -co BIGTIFF=YES", infile, outfile, ..., -9999))
# gdalinfo(outfile)
# 
# 
# 
# ## >> Rivers ####
# ## Source: http://www.soest.hawaii.edu/wessel/gshhg/
# ## See documentation here: /Volumes/uom_data/gsdms_data/lakesrivers/SHAPEFILES.TXT
# ## >> >> Read in files for resolution l, levels 2-9
# river_files <- list.files(file.path(gsdms_data, "lakesrivers/WDBII_shp/l"), full.names = TRUE)
# river_files <- grep(".shp$", river_files, value = TRUE)
# river_files <- grep("L02|L03|L04|L05|L06|L07|L08|L09", river_files, value = TRUE)
# 
# ## >> >> Combine levels 2-9 & write shapefile
# temp <- readOGR(river_files[1])
# for(i in 2:length(river_files)) {
#   temp <- rgeos::gUnion(temp,readOGR(river_files[i]))
# }
# writeOGR(temp2, dsn = file.path(gsdms_data, "lakesrivers/GSHHS/"), 
#          layer = "rivers_l29", driver="ESRI Shapefile")
# 
# ## >> >> Rasterise
# system(paste0("gdal_rasterize -at -burn 1 -ot Byte -tr .0025 .0025 ",
#               file.path(gsdms_data, "lakesrivers/GSHHS/", "rivers_l29.shp "),
#               file.path(gsdms_data, "lakesrivers/GSHHS/", "rivers_l29.tif")))
# 
# ## >> >> Caluclate distance raster
# infile <- file.path(gsdms_data, "lakesrivers/GSHHS/", "rivers_l29.tif")
# gdalinfo(infile)
# 
# outfile <- file.path(dst_folder, "dst_rivers.tif")
# 
# system(sprintf("gdal_proximity.py %s %s -values %d -of GTiff -distunits GEO -nodata %d -co compress=LZW -co BIGTIFF=YES", infile, outfile, ..., -9999))
# gdalinfo(outfile)