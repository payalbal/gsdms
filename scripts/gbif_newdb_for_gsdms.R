
## Load libraries
# install.packages("pacman")
pacman::p_load(DBI, RPostgreSQL, data.table, rgbif, foreach,
               iterators, parallel, doParallel, doMC) 

## Connect to server (named 'con' later in the code)
source("./scripts/connect_to_server.R")


## Get species counts
spcounts <- dbGetQuery(con,"
                       SELECT species, COUNT(*)
                       FROM gbif.clean_gbif
                       GROUP BY species;
                       ")
paste0("Total # species in db = ", dim(spcounts)[1])
paste0("# species to be retained (>= 20 records) = ", dim(spcounts[which(spcounts$count >= 20),])[1])
paste0("# species to be deleted (< 20 records) = ", dim(spcounts[which(spcounts$count < 20),])[1])


## Create species counts table
dbSendQuery(con,"
            DROP TABLE IF EXISTS gbif.species_counts;
            
            CREATE TABLE gbif.species_counts AS
            SELECT species, COUNT(*) AS spcounts
            FROM gbif.clean_gbif
            GROUP BY species;
            
            CREATE INDEX IF NOT EXISTS index_counts_species ON gbif.species_counts (species ASC NULLS LAST);
            ")

## Create table of species with >= 20 occurrences
dbSendQuery(con,"
            DROP TABLE IF EXISTS gbif.species_gte20;
            
            CREATE TABLE gbif.species_gte20 AS
            SELECT *
            FROM gbif.species_counts
            WHERE spcounts >= 20;
            
            CREATE INDEX IF NOT EXISTS index_gte20_species ON gbif.species_gte20 (species ASC NULLS LAST);
            ")

dbSendQuery(con,"
            ALTER TABLE gbif.species_more20
            ADD COLUMNS species_id SERIAL UNIQUE KEY;
            ")

## Create empty table to hold records for species with > = 20 records
dbSendQuery(con, "
            DROP TABLE IF EXISTS gbif.gsdms_gbif;

            CREATE TABLE gbif.gsdms_gbif AS SELECT * FROM gbif.clean_gbif WITH NO DATA;
          ")

## Insert records by species into new table
species_id <- dbGetQuery(con, "
                          SELECT species_id
                          FRPM gbif.species_more20;
                          ")

species_list <- dbGetQuery(con, "
                            SELECT species
                            FRPM gbif.species_more20;
                            ")

year_count <- dbGetQuery(con,"
            SELECT recyear, COUNT(*)
            FROM gbif.clean_gbif
            GROUP BY recyear;
            ")





temp_table <- dbGetQuery(con, 
            paste0("SELECT * FROM gbif.clean_gbif WHERE clean_gbif.species = (SELECT species FROM gbif.species_more20 WHERE species_id = ", i, ");", sep = ""))

dbSendQuery(con, 
            paste0("INSERT INTO gbif.gsdms_gbif SELECT * FROM", temp_table, " ;" , sep = ""))

## Get new species counts
spcounts <- dbGetQuery(con,"
                       SELECT species, COUNT(*)
                       FROM gbif.gsdms_gte20
                       GROUP BY species;
                       ")

## Delete unwanted tables
dbSendQuery(con,"
            DROP TABLE IF EXISTS gbif.species_counts;
            DROP TABLE IF EXISTS gbif.species_gte20;
            ")



## EXTRA
## dopar will not work for inserting records into same db...
  # cl <- makeCluster(30)
  # registerDoParallel(cl)
  # foreach (i=species_list, 
  #          .packages = c('DBI', 'RPostgreSQL', 'foreach',
  #                        'iterators', 'parallel', 'doParallel', 'doMC')) %dopar% {
  #                          drv <- dbDriver("PostgreSQL")
  #                          con <- dbConnect(drv, dbname="qaeco_spatial", user="qaeco", password="Qpostgres15",
  #                                           host="boab.qaeco.com", port="5432")
  #                          dbSendQuery(con, paste0("INSERT INTO gbif.gsdms_gte20 SELECT * FROM gbif.clean_gbif WHERE species = '",
  #                                                  i, "';", sep =""))
  #                        }
  # stopCluster(cl)

## Alternative - parallelisation problem..
spcounts <- dbGetQuery(con,"
                       SELECT species, COUNT(*)
                       FROM gbif.clean_gbif
                       GROUP BY species;
                       ")
species_list <- spcounts[which(spcounts$count >= 20),]$species
for (i in 0:(length(species_list)-1)) {
  dbSendQuery(con, 
              paste0("INSERT INTO gbif.gsdms_gbif SELECT * FROM gbif.clean_gbif WHERE clean_gbif.species = (SELECT species FROM gbif.species_gte20 ORDER BY species OFFSET ", i, " LIMIT 1);", sep = ""))
  
}







