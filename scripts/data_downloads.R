## Download covariate data for gsdms

## Notes
## add data_url and target folder for all download
## target_file flag -O; target_folder flag -P


## Set working environment ---- ####
rm(list = ls())
gc()
# system("ps")
# system("pkill -f R")


## Load libraries
x <- c('bitops', 'RCurl', 'rstudioapi', 'job')
lapply(x, require, character.only = TRUE)
rm(x)

data_dir <- "/tempdata/research-cifs/6300-payalb/uom_data/gsdms_data"
setwd(data_dir)


## GBIF backbone taxonomy ---- ####
##  GBIF Secretariat (2017). GBIF Backbone Taxonomy. Checklist dataset https://doi.org/10.15468/39omei accessed via GBIF.org on 2019-08-26.
system(paste0("curl http://rs.gbif.org/datasets/backbone/backbone-current.zip -o ", data_dir, "/gbif_taxonomy.zip"))
unzip(file.path(data_dir, "gbif_taxonomy.zip"), list = TRUE)
unzip(file.path(data_dir, "gbif_taxonomy.zip"), files = 'Taxon.tsv', exdir = data_path)


## WorlClim ---- ####
## Source: https://www.worldclim.org/data/index.html
## >> Current Worlclim ####
data_url <- "https://biogeo.ucdavis.edu/data/worldclim/v2.1/base/wc2.1_30s_bio.zip"
target_folder <- file.path(data_dir, "worldclim", "current")
system(sprintf("wget -r -np -nH --reject 'index.html*' -e robots=off %s -P %s", data_url, target_folder))

x <- list.files(target_folder)
unzip(x, list = TRUE)
unzip(x, exdir = file.path(data_dir, "srtm"))
file.remove(x)

## >> Future Worldclim ####
models = c("ACCESS-ESM1-5", "BCC-CSM2-MR", "CanESM5", 
           "CanESM5-CanOE", "CMCC-ESM2", "CNRM-CM6-1", 
           "CNRM-CM6-1-HR", "CNRM-ESM2-1", "EC-Earth3-Veg",
           "EC-Earth3-Veg-LR", "FIO-ESM-2-0", "GFDL-ESM4",
           "GISS-E2-1-G", "GISS-E2-1-H", "HadGEM3-GC31-LL",
           "INM-CM4-8", "INM-CM5-0", "IPSL-CM6A-LR",
           "MIROC-ES2L", "MIROC6", "MPI-ESM1-2-HR",
           "MPI-ESM1-2-LR", "MRI-ESM2-0", "UKESM1-0-LL")
ssps = c("ssp126","ssp370","ssp585")
years=c("2021-2040", "2041-2060", "2061-2080","2081-2100")


url_file <- file.path(data_dir, "worldclim", "bio_future_urls.txt")
writeLines(c(""), url_file, sep = "")
url_failed <- file.path(data_dir, "worldclim", "bio_future_urls_failed.txt")
writeLines(c(""), url_failed, sep = "")
for (mod in models){
  for (ssp in ssps){
    for(year in years){
      name <- paste0(paste(c("wc2.1_30s_bioc", mod, ssp, year), collapse = "_"), ".tif")
      source_url <- file.path("https://geodata.ucdavis.edu/cmip6/30s", mod, ssp, name)
      if(url.exists(source_url)) {
        cat(paste(source_url, "\n"), file = url_file, append = T)
      } else {
        warning(c("url does not exist: ", source_url))
        cat(paste(source_url, "\n"), file = url_failed, append = T)
      }
    }
  }
}


## >> >> Download datasets for bioclim
job::job({
  input_files <- file.path(data_dir, "worldclim", "bio_future_urls.txt")
  target_folder <- file.path(data_dir, "worldclim", "future", "raw")
  system(sprintf("wget --no-host-directories --force-directories --input-file=%s -P %s", 
                 input_files, target_folder))
  list.files(target_folder)
})

  # ## Test with url
  # job::job({
  #   data_url <- "https://geodata.ucdavis.edu/cmip6/30s/BCC-CSM2-MR/ssp126/wc2.1_30s_bioc_BCC-CSM2-MR_ssp126_2021-2040.tif"
  #   target_folder <- file.path(data_dir, "worldclim", "future", "raw")
  #   system(sprintf("wget -r -np -nH --reject 'index.html*' -e robots=off %s -P %s", data_url, target_folder))
  # })


## CHELSA ---- ####
## >> Current CHELSA: 1981-2010 ####
## Source: https://chelsa-climate.org/downloads/ > Version 2.1/climatalogies/1981-2010/bio > select 'bio' and 'kg' variables > Download filepaths (RHS) to use with wget
input_files <- file.path( data_dir, "CHELSA/v2_cmip6_1981-2010_envidatS3paths.txt")
target_folder <- file.path( data_dir, "CHELSA")
system(sprintf("wget --no-host-directories --force-directories --input-file=%s -P %s", input_files, target_folder))


## >> Future CHELSA ####
## Source: https://chelsa-climate.org/downloads/ > Version 2.1/climatalogies/> create file paths by years + models + SSPs + 'bio' and 'kg' variables > Download filepaths (RHS) to use with wget
## ref for kg vars: https://www.nature.com/articles/sdata2018214
years=c("2011-2040","2041-2070","2071-2100")
models = c("GFDL-ESM4","IPSL-CM6A-LR","MPI-ESM1-2-HR","MRI-ESM2-0","UKESM1-0-LL")
ssps = c("ssp126","ssp370","ssp585")
bio.vars=1:19
kg.vars = 1:5

## >> >> Bioclim variables ####
## Create download paths and save as txt file
url_file <- file.path(data_dir, "CHELSA", "bio_future_urls.txt")
writeLines(c(""), url_file, sep = "")
url_failed <- file.path(data_dir, "CHELSA", "bio_future_urls_failed.txt")
writeLines(c(""), url_failed, sep = "")
job::job({
  for(year in years){
    for (mod in models){
      for (ssp in ssps){
        for(bio in bio.vars){
          var <- paste0("bio", bio)
          name <- paste(c("CHELSA", var, year, tolower(mod), ssp, "V.2.1.tif"), collapse = "_")
          source_url <- file.path("https://os.zhdk.cloud.switch.ch/envicloud/chelsa/chelsa_V2/GLOBAL/climatologies", year, mod, ssp, "bio",name)
          if(url.exists(source_url)) {
            cat(paste(source_url, "\n"), file = url_file, append = T)
          } else {
            warning(c("url does not exist: ", source_url))
            cat(paste(source_url, "\n"), file = url_failed, append = T)
          }
        }
      }
    }
  }
})

## Download datasets
job::job({
  input_files <- file.path(data_dir, "CHELSA", "bio_future_urls.txt")
  target_folder <- file.path(data_dir, "CHELSA", "future", "raw_bioclimvars")
  system(sprintf("wget --no-host-directories --force-directories --input-file=%s -P %s", 
                 input_files, target_folder))
  list.files(target_folder)
})


## >> >> Köppen-Geiger climate classification variables - NOT YET ####
## Create download paths and save as txt file
url_file <- file.path(data_dir, "CHELSA", "kg_future_urls.txt")
writeLines(c(""), url_file, sep = "")
url_failed <- file.path(data_dir, "CHELSA", "kg_future_urls_failed.txt")
writeLines(c(""), url_failed, sep = "")
job::job({
  for(year in years){
    for (mod in models){
      for (ssp in ssps){
        for(kg in kg.vars){
          var <- paste0("kg", kg)
          name <- paste(c("CHELSA", var, year, tolower(mod), ssp, "V.2.1.tif"), collapse = "_")
          source_url <- file.path("https://os.zhdk.cloud.switch.ch/envicloud/chelsa/chelsa_V2/GLOBAL/climatologies", year, mod, ssp, "bio", name)
          if(url.exists(source_url)) {
            cat(paste(source_url, "\n"), file = url_file, append = T)
          } else {
            warning(c("url does not exist: ", source_url))
            cat(paste(source_url, "\n"), file = url_failed, append = T)
          }
        }
      }
    }
  }
})

## Download datasets
job::job({
  input_files <- file.path(data_dir, "CHELSA", "kg_future_urls.txt")
  target_folder <- file.path(data_dir, "CHELSA", "future", "raw_kgvars")
  system(sprintf("wget --no-host-directories --force-directories --input-file=%s -P %s", input_files, target_folder))
  list.files(target_folder)
})



## SRTM ---- ####
## Source: https://topotools.cr.usgs.gov/gmted_viewer/gmted2010_global_grids.php
data_url <- "https://edcintl.cr.usgs.gov/downloads/sciweb1/shared/topo/downloads/GMTED/Grid_ZipFiles/mn30_grd.zip"
target_file <- file.path(data_dir, "srtm", "srtm.zip")
system(sprintf("wget %s -O %s", data_url, target_file))

unzip(target_file, list = TRUE)
unzip(target_file, exdir = file.path(data_dir, "srtm"))
file.remove(target_file)
## all files within the mn30_grd folder seem to be the same when in as raster
srtm <- file.path(data_dir, "srtm/mn30_grd/w001000.adf")
file.rename(srtm, sub(file_path_sans_ext(basename(srtm)), "srtm", srtm))
srtm <- file.path(data_dir, "srtm/mn30_grd/srtm.adf")


    ## Alternate SRTM source: http://srtm.csi.cgiar.org/srtmdata/ (extent= -180,180,-60,60)
    ##  on weblink select all 30 x 30 tiles and Geo TIFF format > Search
    ## Can only download by tile, so they'll need to be stiched after...
    ## To download a tile..
    # setwd(data_dir)
    # system("wget http://srtm.csi.cgiar.org/wp-content/uploads/files/srtm_30x30/TIFF/N00E090.zip -O srtm_csi.zip")
    # temp <- unzip("srtm_csi.zip", list = TRUE)
    # unzip("srtm_csi.zip", exdir = "./")
    # file.rename(paste0("./", temp$Name), sub(tools::file_path_sans_ext(basename(temp$Name)),"srtm_csi", paste0("./", temp$Name)))
    # file.remove("srtm_csi.zip")
    ## To stitch... [from SK's code, needs to be fixed]
    # tile_names <- as.vector(unlist(as.data.frame(outer(paste0(1:72, "_"), 1:24, FUN = "paste0"))))
    # name <- "srtm"
    # path <- paste0(data_raw, "srtm/")
    # if(is.null(path)) path <- getwd()
    # rs <- raster(nrows = 24, ncols = 72, xmn = -180, xmx = 180,
    #              ymn = -60, ymx = 60)
    # rs[] <- 1:length(rs)
    # points <- xyFromCell(rs, rs[])
    # tileout <- list()
    # for (i in 1:nrow(points)){
    #   # tileout[[i]] <- getData(name = toupper(name), path = path, lon = points[i,1], lat = points[i,2])
    #   tileout[[i]] <- SRTM(lon = points[i,1], lat = points[i,2])
    # }


  ## Alterate SRTM source: https://webmap.ornl.gov/ogc/dataset.jsp?dg_id=10008_1
  ## download by lower resolution and resample or download by smaller extents and stitch tiles together


## Soil ---- ####
## Source (needs login): https://daac.ornl.gov/SOILS/guides/igbp-surfaces.html
## On weblink check the files you want to download and place an order, then copy link from email recieved
## Click 'Get Data' on above link
## Select data, e.g. bulkdens,pawc,soilcarb,totaln and 'Add checked items'
## bulkdens - bulk density (weight of soil in a given volume)
## pawc - profile available water capacity (amount of water that can be stored in a soil profile)
## soilcarb - soil carbon density 
## totaln - total nitrogen density
## Click 'order checked items'
## Get the link to download data from email, which will look like
## this:  https://daac.ornl.gov/orders/<some_number>/
## Run the system(wget...) command with the retrieved link
system("wget -r -np -nH --reject 'index.html*' -e robots=off https://daac.ornl.gov/orders/71a04c738fde75ee64f57118ed466530/")

##  OR...[check if this works, if so use this instead of the above]
system("wget -r -np -nH --reject 'index.html*' -e robots=off https://daac.ornl.gov/orders/ca2c9c7d67f425992259251a578b8669/ -O data/soil.zip")
temp <- unzip("data/soil.zip", list = TRUE)
unzip("data/soil.zip", exdir = "/tempdata/workdir/data")
file.rename(paste0("/tempdata/workdir/data/", temp$Name), sub(file_path_sans_ext(basename(temp$Name)),"soil", paste0("/tempdata/workdir/data/", temp$Name)))
file.remove("data/soil.zip")



## Landuse ----- ####
## >> Source: Hurtt dataset ####
## https://luh.umd.edu/data.shtml - get download links for 'states.nc' from webpage
## >> >> Current: LUH2 v2h Release ####
job::job({
  target_folder <- file.path(data_dir, "landuse", "hurtt")
  system(sprintf("wget -r -np -nH --reject 'index.html*' -e robots=off https://luh.umd.edu/LUH2/LUH2_v2h/states.nc -P %s", target_folder))
})

## >> >> Future: LUH2 v2f Release ####
## >> >> NOTE: - Fix file path in data processing if downloaded this way ####
## RCP2.6 SSP1 , RCP7.0 SSP3, RCP8.5 SSP5
job::job({
  target_folder <- file.path(data_dir, "landuse", "hurtt")
  
  system(sprintf("wget -r -np -nH --reject 'index.html*' -e robots=off https://luh.umd.edu/LUH2/LUH2_v2f/IMAGE/multiple-states_input4MIPs_landState_ScenarioMIP_UofMD-IMAGE-ssp126-2-1-f_gn_2015-2100.nc -P %s", target_folder))
  
  system(sprintf("wget -r -np -nH --reject 'index.html*' -e robots=off https://luh.umd.edu/LUH2/LUH2_v2f/AIM/multiple-states_input4MIPs_landState_ScenarioMIP_UofMD-AIM-ssp370-2-1-f_gn_2015-2100.nc -P %s", target_folder))
  
  system(sprintf("wget -r -np -nH --reject 'index.html*' -e robots=off https://luh.umd.edu/LUH2/LUH2_v2f/MAGPIE/multiple-states_input4MIPs_landState_ScenarioMIP_UofMD-MAGPIE-ssp585-2-1-f_gn_2015-2100.nc -P %s", target_folder))
})



## >> Source: ESA - Climate Research Data Package (CRDP) ####
## http://maps.elie.ucl.ac.be/CCI/viewer/download.php
## Click on 'Copernicus Climate Change Service (C3S) Climate Data Store (CDS)' link
## > Download Data (login if required)
## Select options (here, 2018, v2.1.1, .tar.gz)
## Submit request > Right click Download button and Copy link location
system("wget -np -nH --reject 'index.html*' -e robots=off https://download-0009.copernicus-climate.eu/cache-compute-0009/cache/data3/dataset-satellite-land-cover-3c12d11c-1b01-4054-95ee-339179fbfa76.tar.gz -O /home/payalb/gsdms_r_vol/tempdata/research-cifs/uom_data/gsdms_data/landuse/CCI/dataset-satellite-land-cover-3c12d11c-1b01-4054-95ee-339179fbfa76.tar.gz")

ludir <- file.path( data_dir, "landuse/CCI")
landuse <- list.files(ludir, pattern = ".tar",
                      full.names = TRUE)
system("tar -tvf ", landuse)
system(paste0("tar -xzvf ", landuse, " --directory ", getwd())) ## does not work out of home directory

landuse <- list.files(getwd(), pattern = ".nc",
                      full.names = TRUE)
system(paste0("gdalwarp -of Gtiff -co COMPRESS=LZW -co TILED=YES -ot Byte -te -180.0000000 -90.0000000 180.0000000 90.0000000 -tr 0.002777777777778 0.002777777777778 -t_srs EPSG:4326 NETCDF:", landuse, ":lccs_class ", file.path(ludir, basename(tools::file_path_sans_ext(landuse))), ".tif"))

file.rename(paste0(file.path(ludir, basename(tools::file_path_sans_ext(landuse))), ".tif"), 
            file.path(ludir,"ESA_landuse.tif"))
file.remove(landuse)
## OR https://developers.google.com/earth-engine/datasets/catalog/COPERNICUS_Landcover_100m_Proba-V-C3_Global
## OR https://data.pnnl.gov/group/nodes/dataset/13192


  # ## >> Alternate source (discreet): Copernicus data ####
  # ## ------- reclassify, resample, download by tile and stich tiles
  # ## Link @ GEE: https://developers.google.com/earth-engine/datasets/catalog/COPERNICUS_Landcover_100m_Proba-V_Global (updated)
  # ## Direct download link: https://zenodo.org/communities/copernicus-land-cover/search?page=1&size=20
  # ## Data viewer: https://lcviewer.vito.be/download
  # ## Data processing for global layer @ GEE by MC Loor: https://code.earthengine.google.com/?scriptPath=users%2Fpayalbal%2Fglobal_layers%3Afraclu_mcloor_sklu. See processing in GEE for reclassification scheme.
  # 
  # ## Land use classes used:
  # ## 1	urban
  # ## 2	crop
  # ## 3	forest
  # ## 4	grass
  # ## 5	other 
  # ## 6	NA
  # 
  # ... gee.py script
  # ...download by title
  # 
  # ## Stiching tiles
  # ## Author: Maria del Mar Quiroga, based on: https://stackoverflow.com/a/50235578
  # library(gdalUtils)
  # 
  # # Get a list of all tif tiles downloaded from Google Earth Engine
  # tiffiles <- list.files(path="/tempdata/workdir/data/copernicus/tifs/", pattern="*.tif", full.names = TRUE)
  # 
  # # Build a virtual raster file stitching all tiles
  # # WARNING: This will fail if the file it is trying to write to (output.vrt) already exists
  # gdalbuildvrt(gdalfile = tiffiles, output.vrt = "/tempdata/workdir/troubleshooting/outputs/lu_world.vrt")
  # 
  # # Copy the virtual raster to an actual physical file
  # # WARNING: This takes ~5 minutes to run
  # gdal_translate(src_dataset = "/tempdata/workdir/troubleshooting/outputs/lu_world.vrt", 
  #                dst_dataset = "/tempdata/workdir/troubleshooting/outputs/lu_world.tif", 
  #                output_Raster = TRUE,
  #                options = c("BIGTIFF=YES", "COMPRESSION=LZW"))


  # ## >> Alternate source (fractional data) - Hoskins data ####
  # ... [follow up with CW for single fractional LU map with 12 classes]
  

  # ## >> Alternate source (discrete data 2015 only): Copernicus land use data ####
  # ## source: https://land.copernicus.eu/global/content/release-global-100m-land-cover-maps-2015
  # ## https://zenodo.org/communities/copernicus-land-cover/search?page=1&size=20
  # ## https://zenodo.org/record/3243509#.XiU6d1MzbOQ
  # ## https://lcviewer.vito.be/about
  # ## Reference: Marcel Buchhorn, Bruno Smets, Luc Bertels, Myroslava Lesiv, Nandin-Erdene Tsendbazar, Martin Herold, & Steffen Fritz. (2019). Copernicus Global Land Service: Land Cover 100m: epoch 2015: Globe (Version V2.0.2) [Data set]. Zenodo. http://doi.org/10.5281/zenodo.3243509
  # 
  # ## perhaps also download (for estimates of bias, uncertainty?):
  # ## probaility:
  # # https://zenodo.org/record/3243509/files/ProbaV_LC100_epoch2015_global_v2.0.2_discrete-classification-proba_EPSG-4326.tif
  # ## bias??
  # # https://zenodo.org/record/3243509/files/ProbaV_LC100_epoch2015_global_v2.0.2_DataDensityIndicator_EPSG-4326.tif
  # 
  # system("wget https://zenodo.org/record/3243509/files/ProbaV_LC100_epoch2015_global_v2.0.2_discrete-classification_EPSG-4326.tif -O /tempdata/workdir/data/raw/discrete-classification_EPSG4326.tif")


  # ## >> Alternate Source (discrete data): https://earthexplorer.usgs.gov/ ####
  # ## Click 'Data Sets' tab > Land Cover > GLCC > Click 'Results >>'
  # ## Select one pf the global datasets and download manually through the browser
  # ## Save zip folder @data_dir
  # ## Readme: https://www.usgs.gov/media/images/global-land-cover-characteristics-data-base-version-20
  # unzip(file.path(data_dir, "glccgbe20_tif.zip"), files = 'gbigbpgeo20.tif', exdir = data_dir)
  # landuse <- file.path(data_dir,"gbigbpgeo20.tif")
  # file.rename(landuse, sub(file_path_sans_ext(basename(landuse)), "landuse", landuse))
  # landuse <- file.path(data_dir,"landuse.tif")


  # ## >> Alternate source (discrete data) - Hurtt from Hoskins dataset: http://luh.umd.edu/ ####
  # crop <- raster(file.path(data_dir, "hoskins_landuse" , "CRP_2005/CRP_1km_2005_0ice.bil"))
  # pasture <- raster(file.path(data_dir, "hoskins_landuse" , "PAS_2005/PAS_1km_2005_0ice.bil"))
  # primaryforest <- raster(file.path(data_dir, "hoskins_landuse" , "PRI_2005/PRI_1km_2005_0ice.bil"))
  # secondforest <- raster(file.path(data_dir, "hoskins_landuse" , "SEC_2005/SEC_1km_2005_0ice.bil"))
  # urban <- raster(file.path(data_dir, "hoskins_landuse" , "URB_2005/URB_1km_2005_0ice.bil"))
  


## Roads ---- ####
## Source (requires login): https://sedac.ciesin.columbia.edu/data/set/groads-global-roads-open-access-v1
## Global geodatabase: https://sedac.ciesin.columbia.edu/downloads/data/groads/groads-global-roads-open-access-v1/groads-v1-global-gdb.zip
## Unable to find direct download link
...system("wget -r -np -nH --reject 'index.html*' -e robots=off https://sedac.ciesin.columbia.edu/downloads/data/groads/groads-global-roads-open-access-v1/groads-v1-global-gdb.zip -O data/groads.zip")
dir.create("/tempdata/workdir/data/groads")
unzip("data/groads.zip", exdir = "/tempdata/workdir/data/groads")

## Oceania East: https://sedac.ciesin.columbia.edu/downloads/data/groads/groads-global-roads-open-access-v1/groads-v1-oceania-east-shp.zip
## Notes: https://stackoverflow.com/questions/1324421/how-to-get-past-the-login-page-with-wget
## Unable to find direct download link
...system("wget -r -np -nH --reject 'index.html*' -e robots=off https://sedac.ciesin.columbia.edu/downloads/data/groads/groads-global-roads-open-access-v1/groads-v1-oceania-east-shp.zip -O data/roads_oceaniaeast.zip")
dir.create("/tempdata/workdir/data/roads_oceaniaeast")
unzip("data/roads_oceaniaeast.zip", exdir = "/tempdata/workdir/data/roads_oceaniaeast")

  ## Alternate source: https://www.globio.info/download-grip-dataset


## Lakes & rivers ---- ####
## Source: https://www.worldwildlife.org/pages/global-lakes-and-wetlands-database (Level 3)
target_folder <- "/home/payalb/gsdms_r_vol/tempdata/research-cifs/uom_data/gsdms_data/lakesrivers/GLWD/"
system(sprintf("wget -r -np -nH --reject 'index.html*' -e robots=off https://files.worldwildlife.org/wwfcmsprod/files/Publication/file/9slil0ww7t_GLWD_level3.zip -P %s", target_folder))
unzip(file.path(target_folder, "wwfcmsprod/files/Publication/file/9slil0ww7t_GLWD_level3.zip"), exdir = target_folder)


# ## Source: http://www.soest.hawaii.edu/wessel/gshhg/
# ## Dated from 1996 publication
# system("wget -r -np -nH --reject 'index.html*' -e robots=off http://www.soest.hawaii.edu/pwessel/gshhg/gshhg-shp-2.3.7.zip -O data/lakesrivers.zip")
# dir.create("/tempdata/workdir/data/lakesrivers")
# unzip("data/lakesrivers.zip", exdir = "/tempdata/workdir/data/lakesrivers")

## Other sources
## https://land.copernicus.eu/global/products/wb
## https://data.europa.eu/data/datasets/jrc-floods-floodmapgl_permwb-tif?locale=en based on GLWD
## https://global-surface-water.appspot.com/download

## Built-up areas ---- ####
## Source: https://ghsl.jrc.ec.europa.eu/download.php?ds=bu
## Selct options on LHS: resolution, projection 
## Right click on global download link & paste below
system("wget -r -np -nH --reject 'index.html*' -e robots=off https://cidportal.jrc.ec.europa.eu/ftp/jrc-opendata/GHSL/GHS_BUILT_LDSMT_GLOBE_R2018A/GHS_BUILT_LDS2014_GLOBE_R2018A_54009_1K/V2-0/GHS_BUILT_LDS2014_GLOBE_R2018A_54009_1K_V2_0.zip -O data/builtup.zip")
dir.create("/tempdata/workdir/data/builtup")
unzip("data/builtup.zip", exdir = "/tempdata/workdir/data/builtup")

  ## Alternatre source: https://sedac.ciesin.columbia.edu/data/set/ulandsat-hbase-v1
  ## Unable to find direct download link


## Protected areas ---- #### 
## Source: https://www.protectedplanet.net/en/thematic-areas/wdpa?tab=WDPA
target_folder <- "/home/payalb/gsdms_r_vol/tempdata/research-cifs/uom_data/gsdms_data/protectedareas"
system(sprintf("wget -r -np -nH --reject 'index.html*' -e robots=off https://d1gam3xoknrgr2.cloudfront.net/current/WDPA_Feb2022_Public_shp.zip -P %s", target_folder))
unzip(file.path(target_folder, "wwfcmsprod/files/Publication/file/9slil0ww7t_GLWD_level3.zip"), exdir = target_folder)
unzip(file.path(target_folder, "current/WDPA_Feb2022_Public_shp.zip"), exdir = target_folder)
unlink(file.path(target_folder, "current"), recursive = TRUE)
unlink(file.path(target_folder, "current"))

## >> Unzip tile files
infiles <- list.files(target_folder, pattern = ".zip", full.names = TRUE)
unzip(infiles[1], list = TRUE)

  # lapply(infiles, unzip, exdir = infiles)
  #   ## unnzipped files have same filenames, so files get overwritten
for (i in infiles) {
  unzip(i, exdir = file.path(target_folder, basename(tools::file_path_sans_ext(i)), ""))
}


## >> Stitch tiles
job_script <- "/home/payalb/gsdms_r_vol/tempdata/workdir/gsdms/scripts/tilestitch_wdpa_job.R"
rstudioapi::jobRunScript(job_script, encoding = "unknown", 
                         workingDir = "/home/payalb/gsdms_r_vol/tempdata/workdir/gsdms/",
                         importEnv = FALSE, exportEnv = "")

data_dir <- "/home/payalb/gsdms_r_vol/tempdata/research-cifs/uom_data/gsdms_data"
if(file.exists(file.path(data_dir, "protectedareas/wdpa.tif"))){
  file.remove(file.path(data_dir, "protectedareas/wdpa_temp.vrt"))
}


## Population ---- ####
## Source (requires login): https://sedac.ciesin.columbia.edu/data/set/gpw-v4-population-density-rev11
## Select year (2000), fileformat (GeoTiff), resolution (30 sec) 
## Click 'Create Download' and copy download links
## Unable to find direct download link
target_folder <- "/home/payalb/gsdms_r_vol/tempdata/research-cifs/uom_data/gsdms_data/population/SEDAC"
system(sprintf("wget -r -np -nH --reject 'index.html*' -e robots=off https://sedac.ciesin.columbia.edu/downloads/data/gpw-v4/gpw-v4-population-density-rev11/gpw-v4-population-density-rev11_2000_30_sec_tif.zip -P %s", target_folder))

unzip("data/population.zip", exdir = "/tempdata/workdir/data/population")


## Tree cover ---- ####
## Source: Hansen https://storage.googleapis.com/earthenginepartners-hansen/GFC-2020-v1.8/download.html
## Download urls for 'treecover2000'
input_files <- "/home/payalb/gsdms_r_vol/tempdata/research-cifs/uom_data/gsdms_data/Hansen_treecover/treecover2000_urls.txt"
target_folder <- "/home/payalb/gsdms_r_vol/tempdata/research-cifs/uom_data/gsdms_data/Hansen_treecover/"
system(sprintf("wget --no-host-directories --force-directories --input-file=%s -P %s", input_files, target_folder))

## >> Stitch tiles
job_script <- "/home/payalb/gsdms_r_vol/tempdata/workdir/gsdms/scripts/tilestitch_job.R"
rstudioapi::jobRunScript(job_script, encoding = "unknown", 
                         workingDir = "/home/payalb/gsdms_r_vol/tempdata/workdir/gsdms/",
                         importEnv = FALSE, exportEnv = "")

data_dir <- "/home/payalb/gsdms_r_vol/tempdata/research-cifs/uom_data/gsdms_data"
if(file.exists(file.path(data_dir, "Hansen_treecover/hansen_tree.tif"))){
  file.remove(file.path(data_dir, "Hansen_treecover/hansen_treetemp.vrt"))
}




## EXTRAS ---- ####
# ## WorldClim data - download tiles and stitch to make global layer ---- ####
# ## Source: https://www.worldclim.org/version1
# ## author: Chris Ware
# 
# ## set path to download each zip archive to (will be deleted once used)
# ## must have .zip ext
# zipdst = file.path(data_dir, 'bio_30s.zip')
# 
# ## set dir to extract raster files to
# rasterdst = file.path(data_dir, 'bio_30s/')
# if(!dir.exists(rasterdst)) {
#   dir.create(rasterdst)
# } # end if !dir.exists
# 
# ## set path to write the global raster for each bioclim variable to
# bioclim_dst = paste0(rasterdst, 'bio_current_') # or whatever
# 
# ## Looking at the getData function, it finds which tile a given lat lon  falls witin. 
# ## The tiles are then indexed by row/col which maps to a url. 
# ## Can cut to the chase and just create urls with all combos of row/cols (i.e. tile ids)
# urls_bioclim = NULL
# for (r in 1:5){
#   for (c in 1:12){
#     
#     rc = paste0(r-1, c-1)
#     fn = paste0('bio_', rc, '.zip')
#     thisurl = paste0('https://biogeo.ucdavis.edu/data/climate/worldclim/1_4/tiles/cur/', fn)
#     urls_bioclim = c(urls_bioclim, thisurl)
#     
#   }
# }
# 
# ## check all the urls exist (which hopefully means all the tiles can be downloaded)
#   # require(bitops)
#   # require(RCurl)
# test = lapply(urls_bioclim, url.exists)
# all(unlist(test))
# 
# ## Then, it's possible to loop over the urls, downloading, and unzipping them. 
# download_from_url <- function (urls, zipdst, rasterdst) {
#   
#   for (url in urls){
#     response = tryCatch({
#       download.file(url, destfile = zipdst)
#     }, 
#     error = function(e){e}) # shouldn't get called
#     
#     print(response) # should be 0
#     unzip(zipdst, exdir = rasterdst)
#     file.remove(zipdst)
#   }
# } 
# 
# download_from_url(urls = urls_bioclim, zipdst = zipdst, rasterdst = rasterdst)
# 
# ## rasterdst should now be full of bioclim 1-19 tiles (60 tiles for each var)
# ## loop over each bioclim variable, collect all tiles, and mosaic. 
# for (i in 1:19){
#   f = list.files(rasterdst, pattern = paste0('bio', i, '_'), full.names = TRUE)
#   f = f[grep('.bil$', f)]
#   
#   f = lapply(f, raster)
#   
#   ## specify a function for any overlapping cells (of which there won't be any, 
#   ## but the function requires a function...)
#   f$fun = mean
#   mos = do.call(mosaic, f)
#   
#   this_dst = paste0(bioclim_dst, i, '.tif') # or whatever
#   writeRaster(mos, this_dst)
# }
# 
# 
# ## WorldClim data - Future GCMs
# ## Year: 2070, GCM: BCC-CSM1-1
# url_bc85bi70 = "http://biogeo.ucdavis.edu/data/climate/cmip5/30s/bc85bi70.zip"
# url_bc26bi70 = "http://biogeo.ucdavis.edu/data/climate/cmip5/30s/bc26bi70.zip"
# download_from_url(urls = c(url_bc26bi70, url_bc85bi70), zipdst, rasterdst)
# 
#   # ## From regSSP scripts
#   # ## Global variables
#   # regions <- "global"
#   # rcps <- c("45", "60", "85")
#   # models <- c("BC", "CC", "GS", "HD", "HE", "IP", "MI", "MR", "MC", "MG", "NO")
#   # 
#   # global_mask <- raster(bio_current[1])
#   # global_mask[which(!is.na(global_mask[]))] <- 1
#   # crs(global_mask) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
#   # writeRaster(global_mask, filename = paste0(data_dir, "/mask_global.tif"))
#   # 
#   # 
#   # ## Downlaoad GCM model predictions
#   # urls_bioclim <- c()
#   # for (model_i in models){
#   #   for (rcp_j in rcps){
#   #     fn = paste0(tolower(model_i), rcp_j, 'bi70', '.zip')
#   #     thisurl = paste0('https://biogeo.ucdavis.edu/data/climate/cmip5/30s/', fn)
#   #     urls_bioclim = c(urls_bioclim, thisurl)
#   #   }
#   # }
#   # library(bitops)
#   # library(RCurl)
#   # test = lapply(urls_bioclim, url.exists)
#   # all(unlist(test))
#   # 
#   # zipdst = file.path(gsdms_data, 'gcm_30s.zip')
#   # rasterdst = file.path(gsdms_data, 'gcm_30s/')
#   # if(!dir.exists(rasterdst)) {
#   #   dir.create(rasterdst)
#   # } # end if !dir.exists
#   # for (i in 1:length(urls_bioclim)){
#   #   download_from_url(urls = urls_bioclim[i], zipdst, rasterdst)
#   # }
#   # 
#   # ## Mask data for study regions & create stacks by region and gcms
#   # files_gcm <- list.files(file.path(gsdms_data, 'gcm_30s'), full.names = T, recursive = T)
#   # gcm_reg_path <- file.path(data_dir, 'gcm_reg')
#   # dir.create(gcm_reg_path)
#   # 
#   # for(region in regions){
#   #   for(model_i in models){
#   #     reg_stack <- list()
#   #     file_mod <- files_gcm[grepl(model_i, files_gcm)]
#   #     for(j in 1:length(rcps)){
#   #       file_mod_rcp <- file_mod[grepl(rcps[j], file_mod)]
#   #       temp_stack <- list()
#   #       for(f in 1:length(file_mod_rcp)){
#   #         reg_mask <- readRDS(file.path(data_dir, paste0("mask_", region, ".rds")))
#   #         temp_stack[[f]] <- mask(crop(raster(file_mod_rcp[f]), reg_mask), reg_mask)
#   #       }
#   #       reg_stack[[j]] <- readAll(brick(temp_stack))
#   #     }
#   #     saveRDS(reg_stack, file = paste0(gcm_reg_path, "/", region, "_", model_i, ".rds"))
#   #   }
#   # }
#   # rm(temp_stack, reg_stack)
#   # 
#   # ## Extract cell-wise quartiles across GCM
#   # quartiles <- c("q1", "q2", "q3")
#   # 
#   # for(region in regions){
#   #   gcm <- list.files(gcm_reg_path, full.names = T, pattern = region)
#   #   reg_mask <- readRDS(file.path(data_dir, paste0("mask_", region, ".rds")))
#   #   r <- reg_mask
#   #   inds <- which(!is.na(r[]))
#   #   
#   #   j<-1
#   #   for(j in 1:length(rcps)){
#   #     saveRDS(stack(), file = paste0(data_dir, "/bio", "q1_", rcps[j], "_", region,  ".rds"))
#   #     saveRDS(stack(), file = paste0(data_dir, "/bio", "q2_", rcps[j], "_", region,  ".rds"))
#   #     saveRDS(stack(), file = paste0(data_dir, "/bio", "q3_", rcps[j], "_", region,  ".rds"))
#   #     print(paste0("processing rcp", rcps[j]))
#   #     for(k in 1:19){
#   #       print(paste0("processing bioclim var: ", k))
#   #       bio <- stack()
#   #       for(i in 1:length(models)){
#   #         print(paste0("processing model: ", i))
#   #         dat <- readRDS(gcm[[i]])[[j]]
#   #         bio <- stack(bio, dat[[k]])
#   #       }
#   #       
#   #       print(paste0("getting quartiles..."))
#   #       df1 <- na.omit(as.matrix(getValues(bio)))
#   #       c <-rowQuartiles(df1, probs = c(0.25, 0.5, 0.75))
#   #       for(m in 1:3){
#   #         bioclim <- readRDS(file = paste0(data_dir, "/bio", quartiles[m], "_", rcps[j], "_", region,  ".rds"))
#   #         r[inds] <- c[,m]
#   #         names(r) <- paste0("bio", k)
#   #         saveRDS(readAll(stack(bioclim, r)), file = paste0(data_dir, "/bio", quartiles[m], "_", rcps[j], "_", region,  ".rds"))
#   #       }
#   #     }
#   #   }
#   # }
#   # unlink(gcm_reg_path, recursive=T)
# 
# ## by quartile data in :
# ## files  <- list.files(file.path(rdata_path), pattern = region, full.names = TRUE)
# ## bioq <- files[grepl("bioq", files)]