## ---- [UNDER DEV WITH MDAP TEAM] ---- ##
## Contibutors: Casey Visitn, Usha Nattala




## Load libraries
pacman::p_load(DBI, RPostgreSQL)

## Provide details to connect to server (named 'con' later in the code)
source("./R/connect_to_server.R")

## Set up schema for GBIF data
dbSendQuery(con,"
            CREATE SCHEMA IF NOT EXISTS gbif;
            SET SEARCH_PATH TO gbif;
            ")

## Set up raw GBIF csv as a database...get complete code from Casey
dbSendQuery(con,"
            DROP TABLE IF EXISTS gbif.raw_gbif;
            
            CREATE TABLE gbif.raw_gbif
            (
            gbifid character varying COLLATE pg_catalog.'default' NOT NULL,
            datasetkey character varying COLLATE pg_catalog.'default'
            occurrenceid character varying COLLATE pg_catalog.'default',
            kingdom character varying COLLATE pg_catalog.'default',
            phylum character varying COLLATE pg_catalog.'default',
            class character varying COLLATE pg_catalog.'default',
            'order character varying COLLATE pg_catalog.'default',
            family character varying COLLATE pg_catalog.'default',
            genus character varying COLLATE pg_catalog.'default',
            species character varying COLLATE pg_catalog.'default',
            infraspecificepithet character varying COLLATE pg_catalog.'default',
            taxonrank character varying COLLATE pg_catalog.'default',
            scientificname character varying COLLATE pg_catalog.'default',
            countrycode character varying COLLATE pg_catalog.'default',
            locality character varying COLLATE pg_catalog.'default',
            publishingorgkey character varying COLLATE pg_catalog.'default',
            decimallatitude character varying COLLATE pg_catalog.'default',
            decimallongitude character varying COLLATE pg_catalog.'default',
            coordinateuncertaintyinmeters character varying COLLATE pg_catalog.'default',
            coordinateprecision character varying COLLATE pg_catalog.'default',
            elevation character varying COLLATE pg_catalog.'default',
            elevationaccuracy character varying COLLATE pg_catalog.'default',
            depth character varying COLLATE pg_catalog.'default',
            depthaccuracy character varying COLLATE pg_catalog.'default',
            eventdate character varying COLLATE pg_catalog.'default',
            day character varying COLLATE pg_catalog.'default',
            month character varying COLLATE pg_catalog.'default',
            year character varying COLLATE pg_catalog.'default',
            taxonkey character varying COLLATE pg_catalog.'default',
            specieskey character varying COLLATE pg_catalog.'default',
            basisofrecord character varying COLLATE pg_catalog.'default',
            institutioncode character varying COLLATE pg_catalog.'default',
            collectioncode character varying COLLATE pg_catalog.'default',
            catalognumber character varying COLLATE pg_catalog.'default',
            recordnumber character varying COLLATE pg_catalog.'default',
            identifiedby character varying COLLATE pg_catalog.'default',
            license character varying COLLATE pg_catalog.'default',
            rightsholder character varying COLLATE pg_catalog.'default',
            recordedby character varying COLLATE pg_catalog.'default',
            typestatus character varying COLLATE pg_catalog.'default',
            establishmentmeans character varying COLLATE pg_catalog.'default',
            lastinterpreted character varying COLLATE pg_catalog.'default',
            mediatype character varying COLLATE pg_catalog.'default',
            issue character varying COLLATE pg_catalog.'default'
            )
            WITH (
            OIDS = FALSE
            )
            TABLESPACE pg_default;
            ")


## Insert values into empty table...
dbSendQuery(con,"
            INSERT INTO gbif.raw_gbif (field1, field2, field3) VALUES...
            ...
            ")

## Create indices on database
dbSendQuery(con,"
            DROP INDEX gbif.gbif_index;
            
            CREATE INDEX gbif_index
            ON gbif.raw_gbif USING btree
            (taxonkey COLLATE pg_catalog.'default')
            TABLESPACE pg_default;
            
            DROP INDEX gbif.gbif_index_species;
            
            CREATE INDEX gbif_index_species
            ON gbif.raw_gbif USING btree
            (species COLLATE pg_catalog.'default')
            TABLESPACE pg_default;
            ")