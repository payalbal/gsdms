## Deleteing rows from gbif.data
dbSendQuery(con,"
            DROP TABLE IF EXISTS gbif.speciescounts;
            
            CREATE TABLE gbif.speciescounts AS
            SELECT species, COUNT(*) AS spcounts
            FROM gbif.gsdms_spdata
            GROUP BY species;
            
            CREATE TABLE gbif.species_less20 AS
            SELECT *
            FROM gbif.speciescounts
            WHERE spcounts < 20;
            
            CREATE TABLE gbif.species_more20 AS
            SELECT *
            FROM gbif.speciescounts
            WHERE spcounts >= 20;
            ")

dbSendQuery(con,"
            CREATE INDEX IF NOT EXISTS index_less20_species ON gbif.species_less20 (species);
            ")
dbSendQuery(con,"
            CREATE INDEX IF NOT EXISTS index_more20_species ON gbif.species_more20 (species);
            ")

dbGetQuery(con,"
            SHOW max_worker_processes;
            SHOW temp_buffers;
            SHOW shared_buffers;
            ")

# dbSendQuery(con,"
#             SET max_worker_processes
#             TO 20;
#             ")
# >> ERROR:  parameter "max_worker_processes" cannot be changed without restarting the server

dbSendQuery(con,"
            CREATE INDEX IF NOT EXISTS index_gsdms_species ON gbif.gsdms_spdata (species);
            ")

# ## Delete...
# dbSendQuery(con,"
#             DELETE FROM gbif.gsdms_spdata
#             USING gbif.species_less20
#             WHERE gsdms_spdata.species = species_less20.species;
#             ")


## Create new table of species with >=20 occurrence
dbSendQuery(con,"
            CREATE TABLE gbif.gsdms AS
            SELECT * 
            FROM gbif.gsdms_spdata
            WHERE gsdms_spdata.species IN
            (SELECT species
            FROM gbif.species_more20);
            ")


## https://stackoverflow.com/questions/8290900/best-way-to-delete-millions-of-rows-by-id
## https://stackoverflow.com/questions/47402098/postgresql-slow-delete-from-where-exists

## Trial table
DROP TABLE IF EXISTS gbif.trial; 
CREATE TABLE gbif.trial AS
SELECT *
  FROM gbif.gsdms_spdata
LIMIT 1000000;

DROP TABLE IF EXISTS gbif.trialcounts; 
CREATE TABLE gbif.trialcounts AS
SELECT species, COUNT(*) AS spcounts
FROM gbif.trial
GROUP BY species;

DROP TABLE IF EXISTS gbif.trial_less20; 
CREATE TABLE gbif.trial_less20 AS
SELECT *
  FROM gbif.trialcounts
WHERE spcounts < 20;

DROP TABLE IF EXISTS gbif.trial_more20; 
CREATE TABLE gbif.trial_more20 AS
SELECT *
  FROM gbif.trialcounts
WHERE spcounts >= 20;

DELETE FROM gbif.trial
USING gbif.trial_less20
WHERE trial.species = trial_less20.species;











