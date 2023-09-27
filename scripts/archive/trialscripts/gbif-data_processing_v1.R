
## ------------------ R script to accompany gbif_datacleaning.Rmd ------------------ ##

## GBIF download reference: GBIF.org (12th December 2017) GBIF Occurrence Download https://doi.org/10.15468/dl.kaciwi
## This is a data processing script for GBIF data downloaded as a .csv and saved as a database on a local server. 

## ------------------------------------------------------
## Set up working environment
## ------------------------------------------------------

# Load libraries
library(pacman)
p_load(DBI, RPostgreSQL, data.table, rgbif, foreach, tidyverse, foreach, iterators, parallel, doParallel, doMC) 

## Specify a driver for postgreSQL type database
drv <- dbDriver("PostgreSQL")

## ## Connect to server (named 'con' later in the code)
source("R/connect_to_server.R")

## ------------------------------------------------------
## Catch original number of records
## ------------------------------------------------------

## Records in the online GBIF database
occ_count()

## Records (rows) in the download gbif data (GBIF data citation above)
as.numeric(dbSendQuery(con,"SELECT reltuples::bigint AS estimate
                            FROM pg_class
                            WHERE  oid = 'gbif.raw_data'::regclass;
                            "))


## ------------------------------------------------------
## Run first filter on GBIF data
## ------------------------------------------------------
## Note: Downlaoded GBIF data is stored as a database on a local server.

## Run query for first filter and create new table: gbif.firstfilter
dbSendQuery(con,"
           CREATE TABLE gbif.firstfilter AS
           SELECT *
           FROM gbif.raw_data
           WHERE kingdom = 'Animalia'
           AND decimallatitude IS NOT NULL
           AND decimallongitude IS NOT NULL
           AND year BETWEEN '1950' AND '2017';
           
           ALTER TABLE gbif.firstfilter 
           DROP COLUMN catalognumber,
           DROP COLUMN collectioncode,
           DROP COLUMN identifiedby,
           DROP COLUMN institutioncode,
           DROP COLUMN mediatype,
           DROP COLUMN publishingorgkey,
           DROP COLUMN recordedby,
           DROP COLUMN rightsholder,
           DROP COLUMN typestatus;
           
           ALTER TABLE gbif.firstfilter ADD PRIMARY KEY (gbifid);
           
           ALTER TABLE gbif.firstfilter 
           ALTER COLUMN gbifid TYPE INTEGER USING gbifid::integer;

           ALTER TABLE gbif.firstfilter
           RENAME COLUMN \"order\" TO taxorder;
           
           CREATE INDEX firstfilter_index
           ON gbif.firstfilter USING btree
           (taxonkey COLLATE pg_catalog.default)
           TABLESPACE pg_default;
           ") 


## ------------------------------------------------------
## Convert empty characters in gbif.firstfilter to NULL
##  ------------------------------------------------------

## Obtain field names for gbif.firstfilter
gbif.fields <- dbGetQuery(con,"
                           SELECT column_name 
                           FROM information_schema.columns 
                           WHERE table_schema = 'gbif'
                           AND table_name = 'firstfilter';
                           ")$column_name

## Set empty characters in columns within gbif.firstfilter as NULL
p_load(foreach, iterators, parallel, doParallel, doMC)

registerDoMC(detectCores() - 1) #set up cores for parallel processing (72-1)

foreach (i=gbif.fields) %dopar% { 
  drv <- dbDriver("PostgreSQL")  #Specify a driver for postgreSQL type database
  con <- dbConnect(drv, dbname="qaeco_spatial", user="qaeco", password="Qpostgres15", host="boab.qaeco.com", port="5432")
  dbSendQuery(con, paste("UPDATE gbif.firstfilter SET ",i," = NULL WHERE ",i,"::char = '';", sep=""))
}


# ------------------------------------------------------
## Check new table
# ------------------------------------------------------
# Check number of columns
length(gbif.fields)
dbGetQuery(con,"
            SELECT COUNT(*) from information_schema.columns
            where table_name='firstfilter';
            ")

# Check for NULL values in lat-long
dbSendQuery(con,"
            SELECT *
            FROM gbif.firstfilter
            WHERE decimallatitude IS NOT NULL
            OR decimallongitude IS NOT NULL;
            ")

# Check for NULL values in table...??

## Check min max values for year - too long...
#SELECT MIN(year), MAX(year)
#FROM gbif.firstfilter;


# ------------------------------------------------------
## Change lat-long to numeric
# ------------------------------------------------------
dbSendQuery(con,"
                ALTER TABLE gbif.firstfilter 
                ALTER COLUMN decimallatitude TYPE NUMERIC USING decimallatitude::numeric
                ")

dbSendQuery(con,"
                ALTER TABLE gbif.firstfilter 
                ALTER COLUMN decimallongitude TYPE NUMERIC USING decimallongitude::numeric
                ")


## ------------------------------------------------------
## Create subset data tables
## ------------------------------------------------------

## CLASS: AVES
dbSendQuery(con,"
            CREATE TABLE gbif.aves AS
            SELECT *
            FROM gbif.firstfilter
            WHERE class = 'Aves';

            ALTER TABLE gbif.aves ADD PRIMARY KEY (gbifid);
           
            ALTER TABLE gbif.aves 
            ALTER COLUMN gbifid TYPE INTEGER USING gbifid::integer;
            
            CREATE INDEX aves_index
            ON gbif.aves USING btree
            (taxonkey COLLATE pg_catalog.default)
            TABLESPACE pg_default;
            ")

dbGetQuery(con,"
              SELECT pg_size_pretty(pg_relation_size('gbif.aves'));
              ")

dbGetQuery(con,"
              SELECT reltuples::bigint AS estimate
              FROM pg_class
              WHERE oid = 'gbif.aves'::regclass;
              ")

## ORDER: CARNIVORA
dbSendQuery(con,"
            CREATE TABLE gbif.carnivora AS
            SELECT *
            FROM gbif.firstfilter
            WHERE taxorder = 'Carnivora';
            
            ALTER TABLE gbif.carnivora ADD PRIMARY KEY (gbifid);
            
            ALTER TABLE gbif.carnivora 
            ALTER COLUMN gbifid TYPE INTEGER USING gbifid::integer;
            
            CREATE INDEX carnivora_index
            ON gbif.carnivora USING btree
            (taxonkey COLLATE pg_catalog.default)
            TABLESPACE pg_default;
            ")

## Issue: Remove records where family = "" (2642 records)


## ------------------------------------------------------
## Run second filter on GBIF data
## ------------------------------------------------------

## how to subset > by class ...see example below. 
## download > filter > save; dowload > filter > append...how? ask casey



backbone <-fread("data/gbif_backbone_taxonomy.tsv")
backbone <- backbone[, c("taxonID", "scientificName","taxonRank", "taxonomicStatus", "kingdom", "phylum", "class","order", "family", "genus")]
backbone <- backbone[kingdom == "Animalia",]
backbone <- backbone[phylum == "Chordata",]
chordata_classes <- unique(backbone$class)
temp <- lapply(chordata_classes,function(x)name_suggest(q = x, rank = "class"))
temp2 <- lapply(temp,function(x)occ_search(taxonKey = x$key, return = "meta")$count) ## slow



# Count for kingdom = Animalia
occ_data(taxonKey = 1, hasCoordinate = TRUE, hasGeospatialIssue = FALSE)$meta$count
# Compare with https://www.gbif.org/occurrence/search?has_coordinate=true&has_geospatial_issue=false&taxon_key=1
# notes: using year = c(1950:2017) in the functions does not work...

# Records for phylum = Chordata
occ_data(taxonKey = 44, hasCoordinate = TRUE, hasGeospatialIssue = FALSE, year = c(1950:2017))$meta$count

# Count for family = Felidae using names_suggest
# https://github.com/azizka/Using_biodiversity_data_for_biogeography/wiki/01_Downloading-geographic-occurrence-data-from-GBIF
felidae <- name_suggest(q = "Felidae", rank = "family")
occ_search(taxonKey = felidae$key, return = "meta")$count
felidae_data <- occ_data(taxonKey = felidae$key, hasCoordinate = TRUE, hasGeospatialIssue = FALSE, limit = 5)$meta$count

# Count for multiple families
# file:///Users/payalb/Google%20Drive/UoM/Discovery-trade/doc/ref/setup_occurrence_database.html 
birdfams <- c("Phasianidae","Hydrobatidae","Spheniscidae")
birdfams_gbif <- lapply(birdfams,function(x)name_suggest(q = x, rank = "family"))
birdfams_counts <- lapply(birdfams_gbif,function(x)occ_search(taxonKey = x$key, return = "meta")$count)






## ------------------------------------------------------
## NOTES: 
## ------------------------------------------------------

# GBIF originally has a field named 'order' which needs to be renamed. 
# This is because "order" is an operator in PostgreSQL (e.g. ORDER BY x). So pgAdmin does not recognise "order"/'order'/$$order$$ as a column name. This might work in a query: \"order\"
# e.g. dbSendQuery(con, "UPDATE gbif.firstfilter SET \"order\" = NULL WHERE \"order\"::char = '';")


## Example of subsetting data by family
# backbone <-fread("data/gbif_backbone_taxonomy.tsv")
# backbone <- backbone[, c("taxonID", "scientificName","taxonRank", "taxonomicStatus", "kingdom", "phylum", "class","order", "family", "genus")]
# carnivora.N <- dim(backbone[which(backbone$order == "Carnivora")])
# carnivora.fam <- unique(backbone[which(backbone$order == "Carnivora")]$family)
# foreach (i=carvivora.fam) %do% {
#   dbSendQuery(con, paste("CREATE TABLE gbif.carnivora AS SELECT * FROM gbif.firstfilter WHERE family = '", carnivora.fam[i], "';", sep="") 
