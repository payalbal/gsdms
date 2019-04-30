## GBIF DATA PROCESSING - PART 1

## GBIF data download reference: GBIF.org (12th December 2017) GBIF Occurrence Download https://doi.org/10.15468/dl.kaciwi

## ------------------------------------------------------
## Set up working environment
## ------------------------------------------------------

## Load libraries
library(pacman)
p_load(DBI, RPostgreSQL, data.table, rgbif, foreach, tidyverse, 
       iterators, parallel, doParallel, doMC) 

## Connect to server (named 'con' later in the code)
source("./R/connect_to_server.R")

## ------------------------------------------------------
## Catch original number of records
## ------------------------------------------------------

## Records in the online GBIF database
paste("# records in GBIF database = ", occ_count())

## Records (rows) in the download gbif data (GBIF data citation above)
paste0("# records (estd.) in downloaded GBIF csv file = ",
        as.numeric(dbGetQuery(con,"SELECT reltuples::bigint AS estimate
                                   FROM pg_class
                                   WHERE  oid = 'gbif.raw_data'::regclass;
                                    ")))

## ------------------------------------------------------
## Run first filter on raw GBIF data
## ------------------------------------------------------
  ## Run query for first filter and create new table
dbSendQuery(con,"
            CREATE TABLE gbif.filtered AS
            SELECT
              gbifid,
              occurrenceid, 
              species,
              scientificname,
              countrycode,
              decimallatitude,
              decimallongitude,
              coordinateuncertaintyinmeters,
              coordinateprecision,
              elevation,
              elevationaccuracy,
              \"depth\",
              depthaccuracy,
              eventdate,
              \"year\",
              taxonkey,
              phylum,
              \"class\",
              \"order\",
              \"family\",
              genus,
              specieskey,
              basisofrecord,
              issue
            FROM gbif.raw_data
            WHERE kingdom = 'Animalia'
            AND year BETWEEN '1950' AND '2018'
            AND basisofrecord IN ('HUMAN_OBSERVATION','PRESERVED_SPECIMEN', 'OBSERVATION', 'MATERIAL_SAMPLE', 'MACHINE_OBSERVATION');

            ALTER TABLE gbif.filtered
            RENAME COLUMN \"depth\" TO recdepth;
            ALTER TABLE gbif.filtered
            RENAME COLUMN \"year\" TO recyear;
            ALTER TABLE gbif.filtered
            RENAME COLUMN \"class\" TO taxclass;
            ALTER TABLE gbif.filtered
            RENAME COLUMN \"order\" TO taxorder;
            ALTER TABLE gbif.filtered
            RENAME COLUMN \"family\" TO taxfamily;
            ") 


## Convert empty characters in gbif.filtered to NULL
gbif.fields <- dbGetQuery(con,"
                           SELECT column_name 
                           FROM information_schema.columns 
                           WHERE table_schema = 'gbif'
                           AND table_name = 'filtered';
                           ")$column_name

## Set empty characters in columns within gbif.filtered as NULL
# registerDoMC(detectCores() - 1) # (72-1) -  uses up too much
registerDoMC(30)
foreach (i=gbif.fields) %dopar% { 
  drv <- dbDriver("PostgreSQL")
  con <- dbConnect(drv, dbname="", user="", 
                   password="", host="", port="")
  dbSendQuery(con, paste("UPDATE gbif.filtered SET ",i," = NULL WHERE ",i,
                         "::char = '';", sep=""))
}

  ## Error 1: deadlock detected
  ## see: https://medium.com/@clairesimmonds/postgresql-decoding-deadlocks-183e6a792fd3
  ## Error 2: Error in mclapply(argsList, FUN, mc.preschedule = preschedule, 
  ##                            mc.set.seed = set.seed,: (list) object cannot 
  ##                            be coerced to type 'integer'
  ## see:  https://stackoverflow.com/questions/17355288/parallel-foreach-loops-produce-mclapply-error

## Discard records with NULL values for 'species' or 'scientificname'
dbSendQuery(con,"
            DELETE FROM gbif.filtered
            WHERE species IS NULL
            OR scientificname IS NULL;
            ")

## Discard records with NULL lat long values
dbSendQuery(con,"
            DELETE FROM gbif.filtered
            WHERE decimallatitude IS NULL
            OR decimallongitude IS NULL;
            ")

## Convert lat and long to numeric
dbSendQuery(con,"
            ALTER TABLE gbif.filtered 
            ALTER COLUMN decimallatitude TYPE NUMERIC USING decimallatitude::numeric;
            ")

dbSendQuery(con,"
            ALTER TABLE gbif.filtered
            ALTER COLUMN decimallongitude TYPE NUMERIC USING decimallongitude::numeric;
            ")

## Set primary key and index
dbSendQuery(con,"
            ALTER TABLE gbif.filtered ADD PRIMARY KEY (gbifid);
             
            ALTER TABLE gbif.filtered 
            ALTER COLUMN gbifid TYPE INTEGER USING gbifid::integer;
            
            CREATE INDEX filtered_index
            ON gbif.filtered USING btree
            (taxonkey COLLATE pg_catalog.default)
            TABLESPACE pg_default;
            ")



  # ------------------------------------------------------
  ## Checks
  # ------------------------------------------------------
  ## Check number of columns
  length(gbif.fields) == dbGetQuery(con,"
              SELECT COUNT(*) 
              FROM information_schema.columns
              WHERE table_name='filtered';
              ")
  
  ## Check rows
  as.numeric(dbGetQuery(con,"
                        SELECT reltuples::bigint AS estimate
                        FROM pg_class
                        WHERE  oid = 'gbif.filtered'::regclass;
                        "))
  
  ## Check that NOT NULL lat long values == number of rows
  as.numeric(dbGetQuery(con,"
              SELECT *
              FROM gbif.temp
              WHERE decimallatitude IS NOT NULL
              OR decimallongitude IS NOT NULL;
              "))
  
  ## Check for NULL values in lat-long
  dbGetQuery(con,"
              SELECT *
              FROM gbif.filtered
              WHERE decimallatitude IS NULL
              OR decimallongitude IS NULL;
              ")
  
  ## Check for specific values, e.g. 
  dbGetQuery(con,"
              SELECT * FROM gbif.filtered
              WHERE basisofrecord = 'HUMAN_OBSERVATION';
              ")



  

## -------------------------------------------------------------------------- ##
## -------------------------------------------------------------------------- ##
## SUBSETS
  # ## Example by species: AND species = 'Diceros bicornis'
  # 
  # ## Example by CLASS: AVES
  # dbSendQuery(con,"
  #             CREATE TABLE gbif.aves AS
  #             SELECT *
  #             FROM gbif.filtered
  #             WHERE class = 'Aves';
  # 
  #             ALTER TABLE gbif.aves ADD PRIMARY KEY (gbifid);
  #            
  #             ALTER TABLE gbif.aves 
  #             ALTER COLUMN gbifid TYPE INTEGER USING gbifid::integer;
  #             
  #             CREATE INDEX aves_index
  #             ON gbif.aves USING btree
  #             (taxonkey COLLATE pg_catalog.default)
  #             TABLESPACE pg_default;
  #             ")
  # 
  # ## Example by ORDER: CARNIVORA
  # dbSendQuery(con,"
  #             CREATE TABLE gbif.carnivora AS
  #             SELECT *
  #             FROM gbif.filtered
  #             WHERE taxorder = 'Carnivora';
  #             
  #             ALTER TABLE gbif.carnivora ADD PRIMARY KEY (gbifid);
  #             
  #             ALTER TABLE gbif.carnivora 
  #             ALTER COLUMN gbifid TYPE INTEGER USING gbifid::integer;
  #             
  #             CREATE INDEX carnivora_index
  #             ON gbif.carnivora USING btree
  #             (taxonkey COLLATE pg_catalog.default)
  #             TABLESPACE pg_default;
  #             ")
  # ## Issue: Remove records where family = "" (2642 records)
