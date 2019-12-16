## ---- [UNDER DEV WITH MDAP TEAM] ---- ##
## Contibutors:
## 
## Add script: gbif_newdb_for_gsdms.R


## ------------------------------------------------------
## Set work environment
## ------------------------------------------------------
# install.packages("pacman")
pacman::p_load(DBI, RPostgreSQL, data.table, rgbif, foreach,
       iterators, parallel, doParallel, doMC) 

## Connect to server (named 'con' later in the code)
source("./R/connect_to_server.R")




## ------------------------------------------------------
## Catch original number of records
## ------------------------------------------------------
## All records in the online GBIF database
paste("# Total records in GBIF database as of",  Sys.Date(), "=", occ_count())

  ## Records for mammals, repltiles, amphibians and birds from GBIF db
  ## To search for taxon key, see: https://github.com/gbif/portal-feedback/issues/1820
    # paste("# records in GBIF database for Aves = ", 
    #       occ_count(taxonKey = 212, georeferenced = TRUE))
    # paste("# records in GBIF database for Aves with coordinates = ", 
    #       occ_data(taxonKey = 212, hasCoordinate = TRUE, 
    #                hasGeospatialIssue = FALSE)$meta$count)
    # 
    # paste("# records in GBIF database for Mammalia = ", 
    #       occ_count(taxonKey = 359, georeferenced = TRUE))
    # paste("# records in GBIF database for Mammalia with coordinates = ",
    #       occ_data(taxonKey = 359, hasCoordinate = TRUE, 
    #                hasGeospatialIssue = FALSE)$meta$count)
    # 
    # paste("# records in GBIF database for Reptilia = ", 
    #       occ_count(taxonKey = 358, georeferenced = TRUE))
    # paste("# records in GBIF database for Reptilia with coordinates = ", 
    #       occ_data(taxonKey = 358, hasCoordinate = TRUE, 
    #                hasGeospatialIssue = FALSE)$meta$count)
    # ... add amphibians


## Records (rows) in the download gbif data
paste0("# records (estd.) in downloaded GBIF csv file = ",
       as.numeric(dbGetQuery(con,"SELECT reltuples::bigint AS estimate
                             FROM pg_class
                             WHERE  oid = 'gbif.raw_data'::regclass;
                             ")))




## ------------------------------------------------------
## Processing - Drop unwanted columns and rows
## ------------------------------------------------------
## Run query for first filter and create new table
dbSendQuery(con,"
            DROP TABLE IF EXISTS gbif.clean_gbif;

            CREATE TABLE gbif.clean_gbif AS
            SELECT
            gbifid,
            species,
            scientificname,
            countrycode,
            decimallatitude,
            decimallongitude,
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
            AND basisofrecord IN ('HUMAN_OBSERVATION','PRESERVED_SPECIMEN', 
                                  'OBSERVATION', 'MATERIAL_SAMPLE', 
                                  'MACHINE_OBSERVATION')
            AND issue NOT IN ('ZERO_COORDINATE', 'COORDINATE_INVALID', 
                              'COORDINATE_OUT_OF_RANGE', 
                              'COUNTRY_COORDINATE_MISMATCH', 
                              'COORDINATE_REPROJECTION_FAILED',
                              'COORDINATE_REPROJECTION_SUSPICIOUS',
                              'GEODETIC_DATUM_INVALID',
                              'TAXON_MATCH_FUZZY', 'TAXON_MATCH_HIGHERRANK', 
                              'TAXON_MATCH_NONE',
                              'BASIS_OF_RECORD_INVALID')
            AND species LIKE '% %'
            AND species NOT LIKE '%(sp|sp.|spp|spp.|spec|spec.)%';

            ALTER TABLE gbif.clean_gbif
            RENAME COLUMN \"depth\" TO recdepth;
            ALTER TABLE gbif.clean_gbif
            RENAME COLUMN \"year\" TO recyear;
            ALTER TABLE gbif.clean_gbif
            RENAME COLUMN \"class\" TO taxclass;
            ALTER TABLE gbif.clean_gbif
            RENAME COLUMN \"order\" TO taxorder;
            ALTER TABLE gbif.clean_gbif
            RENAME COLUMN \"family\" TO taxfamily;
            ") 

## Check columns
dbGetQuery(con,"
           SELECT column_name 
           FROM information_schema.columns 
           WHERE table_schema = 'gbif'
           AND table_name = 'clean_gbif';
           ")$column_name


## Set primary key and index on taxonkey
dbSendQuery(con,"
            ALTER TABLE gbif.clean_gbif ADD PRIMARY KEY (gbifid);
            
            ALTER TABLE gbif.clean_gbif 
            ALTER COLUMN gbifid TYPE INTEGER USING gbifid::integer;
            ")

## Create index on species
dbSendQuery(con,"
            CREATE INDEX index_gbif_species
            ON gbif.clean_gbif USING btree (species ASC NULLS LAST);
            ")

## Create index on taxonkey
dbSendQuery(con,"
            CREATE INDEX index_gbif_taxonkey
            ON gbif.clean_gbif USING btree (taxonkey ASC NULLS LAST);
            ")

## Create index on recyear
dbSendQuery(con,"
            CREATE INDEX index_gbif_year
            ON gbif.clean_gbif USING btree (recyear ASC NULLS LAST);
            ")




## ------------------------------------------------------
## Processing - Split GBIF db by year
## ------------------------------------------------------
## Create table for counts by year
dbSendQuery(con,"
            DROP TABLE IF EXISTS gbif.year_counts;
            
            CREATE TABLE gbif.year_counts AS
            SELECT recyear, COUNT(*) AS yearcounts
            FROM gbif.clean_gbif
            GROUP BY recyear;
            ")

year_count <- dbGetQuery(con,"
                        SELECT *
                         FROM gbif.year_counts;
                         ")

year_idx <- c(40, 51, 57, 60, 62, 63, 64, 65, 66, 68)
year_table <- matrix(NA,length(year_idx),1)
rownames(year_table) <- c("1950_89", "1990_00", "2001_06", "2007_09", "2010_11", "2012", "2013", "2014", "2015", "2016_17")
for (i in 1:length(year_idx)){
  if (i == 1){
    year_table[i,] <- sum(year_count[1:year_idx[i],]$yearcounts)
  } else{
    year_table[i,] <- sum(year_count[(year_idx[i-1]+1):year_idx[i],]$yearcounts)
    
  }
}


## Split database by year and create new tables
create.db.year <- function(db_con, schema, parent_db, year1, year2) {
  db_name <- paste0("year", year1, "_", (year2 %% 100))
  RPostgreSQL::dbSendQuery(con, 
                           paste0("DROP TABLE IF EXISTS ", schema, ".", db_name,
                                  "; ", "CREATE TABLE ", schema, ".", db_name, 
                                  " AS SELECT * FROM ", schema, ".", parent_db, 
                                  " WHERE recyear BETWEEN '" , year1,
                                  "' AND '", year2, "'; ", 
                                  "CREATE INDEX index_", db_name, " ON ", schema, 
                                  ".", db_name, " USING btree (recyear ASC NULLS LAST); ",
                                  "CREATE INDEX index_species_", db_name, " ON ", schema, 
                                  ".", db_name, " USING btree (species ASC NULLS LAST); ")
  )
  msg1= "-------------------------------------------------"
  msg2 = paste0("Table : ", db_name, " created at ", Sys.time())
  msg3= "-------------------------------------------------"
  cat(paste(msg1, msg2, msg3, sep = '\n'))
}

create.db.year(db_con = con, schema = "gbif", parent_db = "clean_gbif", year1 = 1950, year2 = 1989)
create.db.year(db_con = con, schema = "gbif", parent_db = "clean_gbif", year1 = 1990, year2 = 2000)
create.db.year(db_con = con, schema = "gbif", parent_db = "clean_gbif", year1 = 2001, year2 = 2006)
create.db.year(db_con = con, schema = "gbif", parent_db = "clean_gbif", year1 = 2007, year2 = 2009)
create.db.year(db_con = con, schema = "gbif", parent_db = "clean_gbif", year1 = 2010, year2 = 2011)
create.db.year(db_con = con, schema = "gbif", parent_db = "clean_gbif", year1 = 2012, year2 = 2012)
create.db.year(db_con = con, schema = "gbif", parent_db = "clean_gbif", year1 = 2013, year2 = 2013)
create.db.year(db_con = con, schema = "gbif", parent_db = "clean_gbif", year1 = 2014, year2 = 2014)
create.db.year(db_con = con, schema = "gbif", parent_db = "clean_gbif", year1 = 2015, year2 = 2015)
create.db.year(db_con = con, schema = "gbif", parent_db = "clean_gbif", year1 = 2016, year2 = 2017)

  


## ------------------------------------------------------
## Processing - Specify empty char as NULL [NOT WORKING]
## ------------------------------------------------------

## Notes
  ## Parallelizing columns versus row (Ivo in email): Postgres does row-level lock for update op. So, if you parallelise operation based on fields (columns), all the processes try to lock the same row at the same time. You need to parallelise the operation based on the rows, by splitting on primary key or dates for example. Construct clean up SQL function for all the (char/text) columns in a single row and apply that on group of rows.
  ## So this won't work....
  # registerDoMC(detectCores() - 1) # (72-1) -  uses up too much
  # cl <- makeCluster(30)
  # registerDoParallel(cl)
  # foreach (i=index_gbif_taxonkey, 
  #          .packages = c('DBI', 'RPostgreSQL', 'foreach',
  #                        'iterators', 'parallel', 'doParallel', 'doMC')) %dopar% {
  #                          drv <- dbDriver("PostgreSQL")
  #                          con <- dbConnect(drv, dbname="qaeco_spatial",
  #                                           user="qaeco", password="", 
  #                                           host="boab.qaeco.com", port="5432")
  #                          dbSendQuery(con, paste("UPDATE gbif.clean_gbif SET ",
  #                                                 i," = NULL WHERE ",i,
  #                                                 "::char = '';", sep=""))
  #                        }
  # stopCluster(cl)
  
  ## Because we get this warning...
  ## Error 1: deadlock detected
  ## see: https://medium.com/@clairesimmonds/postgresql-decoding-deadlocks-183e6a792fd3
  ## Error 2: Error in mclapply(argsList, FUN, mc.preschedule = preschedule, 
  ##                            mc.set.seed = set.seed,: (list) object cannot 
  ##                            be coerced to type 'integer'
  ## see:  https://stackoverflow.com/questions/17355288/parallel-foreach-loops-produce-mclapply-error


## Alternative approach as per year_tables [UNDER DEV]
## get tables name for year_tables
table_names <- dbGetQuery(con, "
                          SELECT * FROM information_schema.tables
                          WHERE table_schema = 'gbif'
                          AND table_name ~ 'year_';
                          ")$table_name

gbif.fields <- dbGetQuery(con,"
           SELECT column_name 
                          FROM information_schema.columns 
                          WHERE table_schema = 'gbif'
                          AND table_name = 'clean_gbif';
                          ")$column_name



registerDoMC(detectCores() - 1)
cl <- makeCluster(30)
registerDoParallel(cl)

foreach (name in table_names,
         .packages = c('DBI', 'RPostgreSQL', 'foreach',
                       'iterators', 'parallel', 'doParallel', 'doMC')) %dopar% {
                         drv <- dbDriver("PostgreSQL")
                         con <- dbConnect(drv, dbname="qaeco_spatial",
                                          user="qaeco", password="",
                                          host="boab.qaeco.com", port="5432")
                         for gbifid... { # keep simple... no need to subset by species/year here
                         for (i in gbif.fields) {
                           drv <- dbDriver("PostgreSQL")
                           con <- dbConnect(drv, dbname="qaeco_spatial",
                                            user="qaeco", password="Qpostgres15", 
                                            host="boab.qaeco.com", port="5432")
                           dbSendQuery(con, paste("UPDATE gbif.clean_gbif SET ",
                                                  i," = NULL WHERE ",i,
                                                  "::char = '';", sep=""))
                         }
                         }
stopCluster(cl)

UPDATE gbif.clean_gbif SET column_name = NULL WHERE column_name::char = '' AND index_gbif_taxonkey = i;




## ------------------------------------------------------
## Processing - Clean by species field
## ------------------------------------------------------
## Discard records with NULL values for 'species' or 'scientificname'
dbSendQuery(con,"
            DELETE FROM gbif.clean_gbif
            WHERE species IS NULL
            OR scientificname IS NULL;
            ")

## Find & remove species  with single apostrophe in names
dbSendQuery(con,"
            DELETE FROM gbif.clean_gbif
            WHERE species LIKE E'%\\'%';
            ")
    ## R may not recognise \'
    ## If so, run query directly in PGAdmin
      # SELECT * FROM gbif.clean_gbif
      # WHERE species LIKE E'%\'%';
    ## See: https://stackoverflow.com/questions/34823158/whats-the-e-before-a-postgres-string




## ------------------------------------------------------
## Processing - Clean by lat & long fields
## ------------------------------------------------------
## Discard records with NULL lat long values
dbSendQuery(con,"
            DELETE FROM gbif.clean_gbif
            WHERE decimallatitude IS NULL
            OR decimallongitude IS NULL;
            ")

## Convert lat and long to numeric
dbSendQuery(con,"
            ALTER TABLE gbif.clean_gbif 
            ALTER COLUMN decimallatitude TYPE NUMERIC USING decimallatitude::numeric;
            ")

dbSendQuery(con,"
            ALTER TABLE gbif.clean_gbif
            ALTER COLUMN decimallongitude TYPE NUMERIC USING decimallongitude::numeric;
            ")




## ------------------------------------------------------
## Processing - Vacuum & reindex
## ------------------------------------------------------
## Clean up space & recreate indices
## See: https://confluence.atlassian.com/kb/optimize-and-improve-postgresql-performance-with-vacuum-analyze-and-reindex-885239781.html
dbSendQuery(con,"
            VACUUM(FULL, ANALYZE, VERBOSE) gbif.clean_gbif
            ")

dbSendQuery(con,"
            REINDEX TABLE gbif.clean_gbif
            ")




## ------------------------------------------------------
## Db checks
## ------------------------------------------------------
## Run ANALYSE
dbGetQuery(con,"
            ANALYZE VERBOSE gbif.clean_gbif
            ")


## Check rows
dbGetQuery(con,"
           SELECT reltuples::bigint AS estimate
           FROM pg_class
           WHERE  oid = 'gbif.clean_gbif'::regclass;
           ")

## Check that NOT NULL lat long values == number of rows
dbGetQuery(con,"
           SELECT COUNT(*)
           FROM gbif.clean_gbif
           WHERE decimallatitude IS NOT NULL
           OR decimallongitude IS NOT NULL;
           ")

## Check for NULL values in lat-long
dbGetQuery(con,"
           SELECT COUNT(*)
           FROM gbif.clean_gbif
           WHERE decimallatitude IS NULL
           OR decimallongitude IS NULL;
           ")

## Check for specific values, e.g. 
dbGetQuery(con,"
           SELECT * FROM gbif.clean_gbif
           WHERE basisofrecord = 'HUMAN_OBSERVATION';
           ")

dbGetQuery(con,"
           SELECT * FROM gbif.clean_gbif
           WHERE basisofrecord = 'LITERATURE';
           ")










