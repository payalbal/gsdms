
## Load libraries
library(pacman)
p_load(DBI, RPostgreSQL) 

## Connect to server (named 'con' later in the code)
source("./R/connect_to_server.R")

## ------------------------------------------------------
## Catch original number of records
## ------------------------------------------------------

## Records (rows) in the download gbif data (GBIF data citation above)
paste0("# records (estd.) in downloaded GBIF csv file = ",
       as.numeric(dbGetQuery(con,"SELECT reltuples::bigint AS estimate
                                   FROM pg_class
                                   WHERE  oid = 'gbif.filtered'::regclass;
                                    ")))

## ------------------------------------------------------
## Run second filter
## ------------------------------------------------------
## Run query for second filter and create new table
dbSendQuery(con,"
                DROP TABLE IF EXISTS gbif.gsdms_spdata;
                
                CREATE TABLE gbif.gsdms_spdata AS
                SELECT 
                  gbifid,
                  species,
                  scientificname,
                  countrycode,
                  decimallatitude,
                  decimallongitude,
                  elevation,
                  elevationaccuracy,
                  recdepth,
                  depthaccuracy,
                  eventdate,
                  recyear,
                  taxonkey,
                  phylum,
                  taxclass,
                  taxorder,
                  taxfamily,
                  genus,
                  specieskey,
                  basisofrecord,
                  issue
                FROM gbif.filtered
                WHERE issue NOT IN ('ZERO_COORDINATE', 'COORDINATE_INVALID', 'COORDINATE_OUT_OF_RANGE', 
                                'COUNTRY_COORDINATE_MISMATCH', 'COORDINATE_REPROJECTION_FAILED',
                                'COORDINATE_REPROJECTION_SUSPICIOUS', 'GEODETIC_DATUM_INVALID',
                                'TAXON_MATCH_FUZZY', 'TAXON_MATCH_HIGHERRANK', 'TAXON_MATCH_NONE',
                                'BASIS_OF_RECORD_INVALID')
                AND species LIKE '% %'
                AND species NOT LIKE '%(sp|sp.|spp|spp.|spec|spec.)%';
                ") 

  ## About '% %' https://stackoverflow.com/questions/20701273/select-only-when-the-field-has-more-than-one-word


## ------------------------------------------------------
## Set primary key and index on species
## ------------------------------------------------------
dbSendQuery(con,"
            ALTER TABLE gbif.gsdms_spdata ADD PRIMARY KEY (gbifid);

            ALTER TABLE gbif.gsdms_spdata
            ALTER COLUMN gbifid TYPE INTEGER USING gbifid::integer;
            ")

dbSendQuery(con,"
            CREATE INDEX species_index
            ON gbif.gsdms_spdata USING btree (species ASC NULLS LAST);
            ")

#             CREATE INDEX gsdms_spdata_index
#             ON gbif.gsdms_spdata USING btree
#             (taxonkey COLLATE pg_catalog.default)
#             TABLESPACE pg_default

## Find & remove species  with single apostrophe in names: R does not recognise \'; run query directly in PGAdmin
SELECT *
FROM gbif.gsdms_spdata
WHERE species LIKE E'%\'%';
  ## Notes: https://stackoverflow.com/questions/34823158/whats-the-e-before-a-postgres-string

DELETE FROM gbif.gsdms_spdata
WHERE species LIKE E'%\'%';

## Display database information:
paste0("# COULUMNS IN DB: ",
       dbGetQuery(con,"
                  SELECT COUNT(*) 
                  FROM information_schema.columns
                  WHERE table_name='gsdms_spdata';
                  "))
paste0("# ROWS IN DB: ",
       dbGetQuery(con,"
                  SELECT reltuples::bigint AS estimate
                  FROM pg_class
                  WHERE  oid = 'gbif.gsdms_spdata'::regclass;
                  "))
paste0("SIZE OF DB : ",
       dbGetQuery(con,"
                  SELECT pg_size_pretty(pg_relation_size('gbif.gsdms_spdata'));
                  "))
  ## Notes: Columns = 21; Rows = 532516192; Size  = 128 GB


## ------------------------------------------------------
## Delete species with < 20 records
## ------------------------------------------------------
## Get species counts
spcounts <- dbGetQuery(con,"
                       SELECT species, COUNT(*)
                       FROM gbif.gsdms_spdata
                       GROUP BY species;
                       ")
paste0("Total # species in db = ", dim(spcounts)[1]) # 314876
paste0("# species retained (>= 20 records) = ", dim(spcounts[which(spcounts$count >= 20),])[1]) ## 108308
paste0("# species deleted (< 20 records) = ", dim(spcounts[which(spcounts$count < 20),])[1]) ## 206568

## List of species names with < 20 records
delete_species <- spcounts[which(spcounts$count < 20),]$species
  # delete_species <- dbGetQuery(con,"
  #                           SELECT species, COUNT(*)
  #                           FROM gbif.gsdms_spdata
  #                           GROUP BY species HAVING COUNT(*) < 20;
  #                           ")$species
saveRDS(delete_species, "./delete_species.rds")

## Remove remaining species with < 20 records
  # dbSendQuery(con, paste0("DELETE FROM gbif.gsdms_spdata WHERE species IN (", 
  #                        paste0(paste("'" , delete_species, "'", sep =""), 
  #                               collapse = ", "), ");", sep = ""))

  dbSendQuery(con, paste0("DELETE FROM gbif.gsdms_spdata WHERE species NOT LIKE '(",
                         paste0(paste("'" , delete_species, "'", sep =""),
                                collapse = ", "), ")';", sep = ""))....
  WHERE  species !~ '(Khairpur|Islamabad|Karachi)';
  where  name like any(array['%Khairpur%','%Islamabad%','%Karachi%']);
  
  
  # dbSendQuery(con, paste0("ALTER TABLE gbif.gsdms_spdata DISABLE TRIGGER ALL; DELETE FROM gbif.gsdms_spdata WHERE species IN (",
  #                         paste0(paste("'" , delete_species, "'", sep =""),
  #                                collapse = ", "), "); ",
  #                         "ALTER TABLE gbif.gsdms_spdata ENABLE TRIGGER ALL;", sep = ""))

  # ## Alternative to delete: create new table
  # dbSendQuery(con, paste0("CREATE TABLE gbif.gsdms AS SELECT * FROM gbif.gsdms_spdata WHERE species NOT IN (", 
  #                         paste0(paste("'" , delete_species, "'", sep =""), 
  #                                collapse = ", "), ");", sep = ""))

## Check for # of reminig species
newcounts <- dbGetQuery(con,"
                       SELECT species, COUNT(*)
                       FROM gbif.gsdms
                       GROUP BY species;
                       ")
paste0("Total # species in db = ", dim(newcounts)[1])
paste0("# species retained (>= 20 records) = ", dim(newcounts[which(newcounts$count >= 20),])[1])
paste0("# species deleted (< 20 records) = ", dim(newcounts[which(newcounts$count < 20),])[1])

## Display database information
paste0("# ROWS IN NEW DB: ",
       dbGetQuery(con,"
                  SELECT reltuples::bigint AS estimate
                  FROM pg_class
                  WHERE  oid = 'gbif.gsdms'::regclass;
                  "))
paste0("SIZE OF NEW DB: ",
       dbGetQuery(con,"
                  SELECT pg_size_pretty(pg_relation_size('gbif.gsdms'));
                  "))

## Get species list for new database
gbif_species <- newcounts$species
  # gbif_species <- dbGetQuery(con,"
  #                            SELECT species, COUNT(*)
  #                            FROM gbif.gsdms
  #                            GROUP BY species;
  #                            ")$species
saveRDS(gbif_species, "./gbif_species.rds")





## ------------------------------------------------------
## EXTRAS
## ------------------------------------------------------
## Test time to delete versus number of rows in a table
dbSendQuery(con,"
            DROP TABLE IF EXISTS gbif.trial;
            CREATE TABLE gbif.trial AS
            SELECT *
            FROM gbif.gsdms_spdata
            LIMIT 10000000;
            ")

## Get species counts
spcounts <- dbGetQuery(con,"
                       SELECT species, COUNT(*)
                       FROM gbif.trial
                       GROUP BY species;
                       ")

## List of species names with < 20 records
delete_species <- spcounts[which(spcounts$count < 20),]$species

## Delete species
starttime <- Sys.time()
dbSendQuery(con, paste0("DELETE FROM gbif.trial WHERE species IN (", 
                        paste0(paste("'" , delete_species, "'", sep =""), 
                               collapse = ", "), ");", sep = ""))
end.time <- Sys.time() - starttime
end.time

## Notes
  ## 1000000 ~ 50 secs
  ## 5000000 ~ 13 mins
  ## 2500000 ~ 6 mins
  ## 10000000 ~ 30 mins


## 
library(pacman)
p_load(DBI, RPostgreSQL)

setwd("./gbif_analysis")
source("./connect_to_server.R")


## Test time to delete based on query used
## Create trial tables
dbSendQuery(con,"
            DROP TABLE IF EXISTS gbif.trial;
            CREATE TABLE gbif.trial AS
            SELECT *
            FROM gbif.gsdms_spdata
            LIMIT 1000000;
            ")
dbSendQuery(con,"
            DROP TABLE IF EXISTS gbif.trial2;
            CREATE TABLE gbif.trial AS
            SELECT *
            FROM gbif.trial
            ")
dbSendQuery(con,"
            DROP TABLE IF EXISTS gbif.trial3;
            CREATE TABLE gbif.trial AS
            SELECT *
            FROM gbif.trial
            ")
delete_species <- dbGetQuery(con,"
                             SELECT species, COUNT(*)
                             FROM gbif.trial
                             GROUP BY species HAVING COUNT(*) <= 20;
                             ")$species

## Run delete queries
start.time <- Sys.time()
dbSendQuery(con, paste0("DELETE FROM gbif.trial WHERE species IN (",
                        paste0(paste("'" , delete_species[1:2], "'", sep =""),
                               collapse = ", "), ");", sep = ""))
end.time <- Sys.time() - start.time
paste0("species IN query:", end.time)

start.time <- Sys.time()
dbSendQuery(con, paste0("DELETE FROM gbif.trial2 WHERE species ~ '(", 
                        paste0(delete_species[1:2], collapse = "|"), ")';", sep = ""))
end.time <- Sys.time() - start.time
paste0("species ~ query:", end.time)

start.time <- Sys.time()
dbSendQuery(con, paste0("DELETE FROM gbif.trial3 WHERE species LIKE ANY(ARRAY[",
                        paste0(paste("'" , delete_species[1:2], "'", sep =""),
                               collapse = ", "), "]);", sep = ""))
end.time <- Sys.time() - start.time
paste0("species LIKE ANY ARRAY query:", end.time)

## Notes 
  ## IN query and LIKE ANY query take similar amoiunts of time
  ## ~ query takes much longer than the other two





## Check 
dbGetQuery(con,"
          SELECT species, COUNT(*)
          FROM gbif.trial
          GROUP BY species;
          ")

## Deletes everything...
# DELETE FROM gbif.trial
# WHERE EXISTS
# (SELECT
# FROM gbif.trial
# GROUP BY species HAVING COUNT(*) < 20);

## Doesn't delete anything...    
# DELETE FROM gbif.trial
# WHERE EXISTS
# (SELECT
# FROM gbif.trial
# HAVING COUNT(DISTINCT species) < 20);


## Try query to create table of species with >=20 - doesn't work, still get species with < 20 records
# CREATE TABLE gbif.trial AS
# SELECT *
# FROM gbif.gsdms_spdata
# WHERE EXISTS
# (SELECT
# FROM gbif.gsdms_spdata
# HAVING COUNT(DISTINCT species) >= 20)
# LIMIT 100;
# 
#SELECT species, COUNT(*)
# FROM gbif.trial
# GROUP BY species HAVING COUNT(*) < 20;



 


## Remove species  with single apostrophe in names: R does not recognise \'; run query directly in PGAdmin
# problem_species <- dbGetQuery(con,"
#                                SELECT *
#                                FROM gbif.gsdms_spdata 
#                                WHERE species LIKE E'%\'%';
#                                ")
SELECT *
  FROM gbif.gsdms_spdata
WHERE species LIKE E'%\'%';
## Notes: https://stackoverflow.com/questions/34823158/whats-the-e-before-a-postgres-string

problem_species <- c("Scalmicauda o'neili", "Perisama d'orbignyi")

## See number of records for species with single apostrophe in name: 
for (i in problem_species) {
  print(spcounts[grep(i, spcounts$species),])
}

## Delete records for species with with single apostrophe in name (if < 20 records): R does not recognise \'; run query directly in PGAdmin
# dbGetQuery(con, paste0("DELETE FROM gbif.gsdms_spdata WHERE species IN (", 
#                        paste0(paste("'" , problem_species, "'", sep =""), 
#                               collapse = ", "), ");", sep = ""))

DELETE FROM gbif.gsdms_spdata
WHERE species LIKE E'%\'%';

## Remove single apostrophe names from species name vector 
delete_species <- setdiff(delete_species, problem_species)
# delete_species <- setdiff(delete_species, delete_species[grep("Scalmicauda o'neili|Perisama d'orbignyi", delete_species)])