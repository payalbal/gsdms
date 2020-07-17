## gsdms table on postgresql server

## Setup
x <- c("DBI", "RPostgreSQL")
lapply(x, require, character.only = TRUE)
setwd("./gbif_analysis")
source("./connect_to_server.R")

## Check connection & display system info
dbGetQuery(con,"
            SHOW max_worker_processes;
            ")
dbGetQuery(con,"
           SHOW temp_buffers;
           ")
dbGetQuery(con,"
           SHOW shared_buffers;
           ")

## Run query for second filter and create new table
dbSendQuery(con,"
            DROP TABLE IF EXISTS gbif.gsdms;
            
            CREATE TABLE gbif.gsdms AS
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
            'BASIS_OF_RECORD_INVALID');
            ") 

## Set primary key
dbSendQuery(con,"
            ALTER TABLE gbif.gsdms
            ALTER COLUMN gbifid TYPE INTEGER USING gbifid::integer;

            ALTER TABLE gbif.gsdms ADD PRIMARY KEY (gbifid);
            ")



## Delete species with problematic names
## Create species counts table
dbSendQuery(con,"
            DROP TABLE IF EXISTS gbif.species_counts;
            
            CREATE TABLE gbif.species_counts AS
            SELECT species, COUNT(*) AS spcounts
            FROM gbif.gsdms
            GROUP BY species;
            
            CREATE INDEX IF NOT EXISTS index_counts_species ON gbif.species_counts (species ASC NULLS LAST);
            ")
## 2. Find & remove species  with single apostrophe in names: R does not recognise \'; run query directly in PGAdmin
SELECT *
FROM gbif.species_counts
WHERE species LIKE E'%\'%';
## Notes: https://stackoverflow.com/questions/34823158/whats-the-e-before-a-postgres-string

DELETE FROM gbif.gsdms
WHERE species LIKE E'%\'%';

## Find & remove species with incomplete names (e.g. sp|sp.|spp|spp.|spec|spec.)
oddnames <- dbGetQuery(con,"
                        SELECT species
                        FROM gbif.gsdms
                        WHERE species NOT LIKE '% %'
                        AND species LIKE '%(sp|sp.|spp|spp.|spec|spec.)%'
                        ")
if (!is.null(oddnames)){
  dbSendQuery(con,"
              DELETE FROM gbif.gsdms
              WHERE species LIKE '% %'
              AND species NOT LIKE '%(sp|sp.|spp|spp.|spec|spec.)%'
              ")  
}
## Note: ## About '% %' https://stackoverflow.com/questions/20701273/select-only-when-the-field-has-more-than-one-word  

## Create index on species
dbSendQuery(con,"
            CREATE INDEX IF NOT EXISTS index_gsdms_species
            ON gbif.gsdms USING btree (species ASC NULLS LAST);
            ")



## Create species counts table
dbSendQuery(con,"
            DROP TABLE IF EXISTS gbif.species_counts;
            
            CREATE TABLE gbif.species_counts AS
            SELECT species, COUNT(*) AS spcounts
            FROM gbif.gsdms
            GROUP BY species;

            CREATE INDEX IF NOT EXISTS index_counts_species ON gbif.species_counts (species ASC NULLS LAST);
            ")

## Create list of species with < 20 occurrences
dbSendQuery(con,"
            DROP TABLE IF EXISTS gbif.species_less20;
            
            CREATE TABLE gbif.species_less20 AS
            SELECT *
            FROM gbif.species_counts
            WHERE spcounts < 20;

            CREATE INDEX IF NOT EXISTS index_less20_species ON gbif.species_less20 (species ASC NULLS LAST);
            ")

## Create list of species with >= 20 occurrences
dbSendQuery(con,"
            DROP TABLE IF EXISTS gbif.species_more20;

            CREATE TABLE gbif.species_more20 AS
            SELECT *
            FROM gbif.species_counts
            WHERE spcounts >= 20;

            CREATE INDEX IF NOT EXISTS index_more20_species ON gbif.species_more20 (species ASC NULLS LAST);
            ")

## Create new table with species >= 20 occurrences
# dbSendQuery(con,"
#             CREATE TABLE gbif.gsdmsv2 AS
#             SELECT * 
#             FROM gbif.gsdms
#             WHERE gsdms.species IN
#             (SELECT species
#             FROM gbif.species_more20);
#             ")

species_list <- dbGetQuery(con, "
                              SELECT *
                              FROM gbif.species_counts
                              WHERE spcounts >= 20;
                              ")$species

species100 <- sample(species_list, 100)

dbSendQuery(con, paste0("CREATE TABLE gbif.gsdmsv2 AS SELECT * FROM gbif.gsdms WHERE species = '", 
                        species_list[1], "';", sep =""))
registerDoMC(30)
foreach (i = species_list[-1]) %dopar% {
  drv <- dbDriver("PostgreSQL")
  con <- dbConnect(drv, dbname="qaeco_spatial", user="qaeco", password="Qpostgres15",
                   host="boab.qaeco.com", port="5432")
  dbSendQuery(con, paste0("INSERT INTO gbif.gsdmsv2 SELECT * FROM gbif.gsdms WHERE species = '",
                          i, "';", sep =""))
}

## OR
-- Create an empty table
CREATE TABLE species_gte_20 AS SELECT * FROM gsdms_spdata WITH NO DATA;

-- insert in batches for each species, change OFFSET 0 incrementally to OFFSET [no of species > 20]
INSERT INTO spdata_gte_20
SELECT *
  FROM gbif.gsdms_spdata
WHERE gsdms_spdata.species = (SELECT species FROM species_more20 ORDER BY species OFFSET 0 LIMIT 1); // get one species at a time


## Set Primary key
dbSendQuery(con,"
            ALTER TABLE gbif.gsdmsv2 ADD PRIMARY KEY (gbifid);
            ")

## Create index on species
dbSendQuery(con,"
            CREATE INDEX IF NOT EXISTS index_gsdmsv2_species
            ON gbif.gsdmsv2 USING btree (species ASC NULLS LAST);
            ")

## Consider partitioning - manual or automatic (Ivo's suggestion)     



## Get info from nnew table...
## EXTRAS
## Counting by class
CREATE TABLE gbif.birdcounts AS
SELECT species, COUNT(*) AS spcounts
FROM gbif.gsdms_spdata
WHERE taxclass = 'Aves'
GROUP BY species HAVING COUNT(*) >= 20;

aves_N <- dbGetQuery(con,"
                      SELECT COUNT(DISTINCT species) 
                      FROM gbif.birdcounts;
                      ")

## Get list of taxa in db
taxlist <- dbGetQuery(con,"
                        SELECT DISTINCT taxclass AS taxnames
                        FROM gbif.gsdms_spdata
                        ORDER BY taxnames;
                        ")
