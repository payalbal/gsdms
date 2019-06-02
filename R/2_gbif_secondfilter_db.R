## GBIF DATA PROCESSING - PART 2
## Data citation: GBIF.org (12th December 2017) GBIF Occurrence Download https://doi.org/10.15468/dl.kaciwi

## ------------------------------------------------------
## Set up working environment
## ------------------------------------------------------

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

## Set primary key and index
dbSendQuery(con,"
            ALTER TABLE gbif.gsdms_spdata ADD PRIMARY KEY (gbifid);
            
            ALTER TABLE gbif.gsdms_spdata 
            ALTER COLUMN gbifid TYPE INTEGER USING gbifid::integer;
            
            CREATE INDEX gsdms_spdata_index
            ON gbif.gsdms_spdata USING btree
            (taxonkey COLLATE pg_catalog.default)
            TABLESPACE pg_default
            ")

# ## Remove species with < 20 occurreces
# dbSendQuery(con,"
#             DELETE FROM gbif.gsdms_spdata
#             WHERE EXISTS
#             (SELECT
#             FROM gbif.gsdms_spdata
#             GROUP BY species HAVING COUNT(*) < 20);
#             ")


## ------------------------------------------------------
## Checks
## ------------------------------------------------------
## Check number of columns
dbGetQuery(con,"
           SELECT COUNT(*) 
           FROM information_schema.columns
           WHERE table_name='gsdms_spdata';
           ")

## Check rows
dbGetQuery(con,"
           SELECT reltuples::bigint AS estimate
           FROM pg_class
           WHERE  oid = 'gbif.gsdms_spdata'::regclass;
           ")

## Check size of new table
dbGetQuery(con,"
           SELECT pg_size_pretty(pg_relation_size('gbif.gsdms_spdata'));
           ")

## Get species counts
spcounts <- dbGetQuery(con,"
                       SELECT species, COUNT(*)
                       FROM gbif.gsdms_spdata
                       GROUP BY species;
                       ")
dim(spcounts)[1]
dim(spcounts[which(spcounts$count >= 20),])[1]
dim(spcounts[which(spcounts$count < 20),])[1]

dbGetQuery(con,"
           SELECT species, COUNT(*)
           FROM gbif.gsdms_spdata
           GROUP BY species HAVING COUNT(*) < 20;
           ")



## ------------------------------------------------------
## EXTRAS - before removing species with <20 occurrences
## ------------------------------------------------------

## Get species counts
priorcounts <- dbGetQuery(con,"
                       SELECT species, COUNT(*)
                       FROM gbif.gsdms_spdata
                       GROUP BY species;
                       ")
  dim(priorcounts)[1] #314878
  dim(priorcounts[which(priorcounts$count >= 20),])[1] #108308
  dim(priorcounts[which(priorcounts$count < 20),])[1] #206570
  splist <- priorcounts[which(priorcounts$count >= 20),]$species
  length(splist) #108308

## Counts for species with >= 20 occurreces
newcounts <- dbGetQuery(con,"
                        SELECT species, COUNT(*)
                        FROM gbif.gsdms_spdata
                        GROUP BY species HAVING COUNT(*) >= 20;
                        ")
  dim(newcounts) #108308

## Counts for species with < 20 occurrences
lesscounts <- dbGetQuery(con,"
                         SELECT species, COUNT(*)
                         FROM gbif.gsdms_spdata
                         GROUP BY species HAVING COUNT(*) < 20;
                         ")
  dim(lesscounts) #206570

  ## Checks
  sum(sort(splist) == sort(newcounts$species)) #108308
  dim(priorcounts)[1] == sum(dim(newcounts)[1], dim(lesscounts)[1]) #TRUE
  

  
  
    dbSendQuery(con,"
                CREATE VIEW gbif.temp AS
                SELECT *
                FROM gbif.gsdms_spdata
                WHERE EXISTS
                (SELECT
                FROM gbif.gsdms_spdata
                GROUP BY species HAVING COUNT(*) >= 20)
                LIMIT 100;
                ")
 


