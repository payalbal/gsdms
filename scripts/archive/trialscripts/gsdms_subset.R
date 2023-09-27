## Create subset GBIF table on postgresql server

## Specify parameters for subsetting GBIF data
n_species <- 100
n_records <- 1000


## Setup
x <- c("DBI", "RPostgreSQL", "foreach", "iterators", "parallel", "doParallel", "doMC")
lapply(x, require, character.only = TRUE)
setwd("./gbif_analysis")
source("./connect_to_server.R")

## Get species names to model
species_list <- dbGetQuery(con, paste0("SELECT * FROM gbif.species_counts WHERE spcounts >= ",
                                       n_records, ";", sep =""))$species
# species100 <- sample(species_list, n_species)
species_sub <- species_list[1:n_species]

## Get data for modelling
dbSendQuery(con, paste0("CREATE TABLE gbif.gsdms1000 AS SELECT * FROM gbif.gsdms_spdata WHERE species = '", 
                        species_sub[1], "';", sep =""))
registerDoMC(30)
foreach (i = species_sub[-1]) %dopar% {
  drv <- dbDriver("PostgreSQL")
  con <- dbConnect(drv, dbname="qaeco_spatial", user="qaeco", password="Qpostgres15",
                   host="boab.qaeco.com", port="5432")
  dbSendQuery(con, paste0("INSERT INTO gbif.gsdms1000 SELECT * FROM gbif.gsdms_spdata WHERE species = '",
                          i, "';", sep =""))
}



# ## Set primary key and index
# dbSendQuery(con,"
#             ALTER TABLE gbif.gsdms1000 ADD PRIMARY KEY (gbifid);
#             
#             CREATE INDEX IF NOT EXISTS index_gsdms1000_species
#             ON gbif.gsdms1000 USING btree (species ASC NULLS LAST);
#             ")