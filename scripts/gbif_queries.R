## gbif_queries
## Function to extract and prepare species data from GBIF db for ppm fitting
## $$ New package: Australian Mammals Species List as an excel sheet? ####

## Set working environment ####
rm(list = ls())
gc()

x = c("DBI", "RPostgreSQL", "data.table", "rpostgis", "sp")
lapply(x, require, character.only = TRUE)
rm(x)

## Connect to server
source("~/gsdms_r_vol/tempdata/workdir/gbifprocessing/scripts/connect_to_server.R")

## Check db has PostGIS installed
pgPostGIS(con)
pgListGeom(con, geog = TRUE)

## Get db column names
db_cols <- dbGetQuery(con, sprintf("
                      SELECT column_name
                      FROM information_schema.columns
                      WHERE table_schema = 'public'
                      AND table_name = '%s';", dbname))$column_name


## Get species list
splist <- fread("/tempdata/research-cifs/uom_data/gsdms_data/gbif/gbif_splist.csv")


## Record original number of records for species in clean db
sp <- splist[1]

dbname = paste0("spcounts_", gsub("gbif_", "", tolower(sp$taxclass)))
sp_N <- dbGetQuery(con, sprintf("
                    SELECT spcounts
                    FROM %s
                    WHERE species = '%s';",
                    dbname, sp$species))


## Get data for species fom GBIF db
dbname = paste0("gbif_", tolower(sp$taxclass))
spdat <- as.data.table(dbGetQuery(con, sprintf("
                    SELECT *
                    FROM %s
                    WHERE species = '%s';",
                    dbname, sp$species)))


dat2 <- dbReadDataFrame(dbGetQuery(con, sprintf("
                    SELECT *
                    FROM %s
                    WHERE species = '%s';",
                    dbname, sp$species)))



# filter by extent (acc to global/aus mask)
# filter by regions used for fitting/prediction
# remove duplicates
# create hold out versus test data (for fit vs prediction)
# if species with <=20 records (add to ppm code)
