## GBIF DATA PROCESSING - PART 3

## ------------------------------------------------------
## Run third filter on GBIF data
## ------------------------------------------------------

library(pacman)
p_load(DBI, RPostgreSQL, data.table, sp, raster)

data_path <- "/Volumes/discovery_data/data/"
backbone <-fread(paste0(data_path, "raw/", "biodiversity/GBIF/backbone/Taxon.tsv"))

global_mask <- raster::raster(paste0(data_path, "raw/climate/wc10/bio1.bil"))
global_mask[which(!is.na(global_mask[]))] <- 1
global_mask[which(is.na(global_mask[]))] <- NA

## Read in data from server
source("R/connect_to_server.R")
dat <- as.data.table(dbGetQuery(con,"
                  SELECT *
                  FROM gbif.iucn_species;
                  "))

## Run second filter
source("./R/filter.gbif.R")
gbifdat <- filter.gbif(dat, backbone, domain.mask = global_mask, 
                            remove_duplicates = FALSE, 
                            output_name = "gbif_iucnsp",
                            output_folder = "./output")
paste0("# species in data = ", length(unique(gbifdat$species)))
as.matrix(table(gbifdat$species))
