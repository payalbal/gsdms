## GBIF DATA PROCESSING - PART 2
## Extract IUCN species from GBIF and create new database

library(pacman)
p_load(DBI, RPostgreSQL, foreach, iterators, parallel, doParallel, doMC)
source("R/connect_to_server.R")

## Get IUCN species list
# data_path <- "/Volumes/payal_umelb/data/"
# source("./R/filter_iucn_data.R")
# dat <- readOGR(dsn=paste0(data_path, "raw/biodiversity/IUCN/rangemaps/TERRESTRIAL_MAMMALS"), layer = "TERRESTRIAL_MAMMALS")
# iucn_species <- as.vector(filter.iucn(dat, output_folder = "./output", species_list = TRUE))


## Load IUCN species list created above. 
iucn_species <- as.vector(read.table("./output/2019-04-03_filtered_iucn.txt")[[1]])
  # iucn_species <- c("Panthera pardus" , "Lophuromys melanonyx”, "Macaca arctoides”, “Macaca nigra")

## Extract GBIF data based on IUCN species list
iucn_species <- iucn_species[1:4]
dbSendQuery(con, paste0("CREATE TABLE gbif.iucn_species AS SELECT * FROM gbif.filtered WHERE species = '", 
                        iucn_species[1], "';", sep =""))
registerDoMC(30)
foreach (i = iucn_species[-1]) %dopar% {
  drv <- dbDriver("PostgreSQL")
  con <- dbConnect(drv, dbname="qaeco_spatial", user="qaeco", password="Qpostgres15",
                   host="boab.qaeco.com", port="5432")
  dbSendQuery(con, paste0("INSERT INTO gbif.iucn_species SELECT * FROM gbif.filtered WHERE species = '",
                          i, "';", sep =""))
}

## Set primary key and index
dbSendQuery(con,"
            ALTER TABLE gbif.iucn_species ADD PRIMARY KEY (gbifid);
            
            CREATE INDEX iucn_species_index
            ON gbif.iucn_species USING btree
            (taxonkey COLLATE pg_catalog.default)
            TABLESPACE pg_default;
            ")


## Find # records in gbif.filtered per species in iucn species list
registerDoMC(30)
species_counts <- foreach (i = iucn_species, .combine=rbind) %dopar% {
  drv <- dbDriver("PostgreSQL")
  con <- dbConnect(drv, dbname="qaeco_spatial", user="qaeco", password="Qpostgres15",
                   host="boab.qaeco.com", port="5432")
  dbGetQuery(con, paste0("SELECT species, COUNT(*) FROM gbif.filtered WHERE species = '",
                         i, "' GROUP BY species;", sep =""))
}


## Load and add zero count species
nospecies <- which(!(as.factor(iucn_species) %in% species_counts[,1]))
temp <- as.data.frame(cbind(iucn_species[nospecies], rep(0,length(nospecies))))
colnames(temp) <- colnames(species_counts)
temp <- rbind(species_counts,temp)
iucn_species_counts <- temp[match(iucn_species,temp$species),]
write.csv(iucn_species_counts, "./output/iucn_species_counts.csv", row.names=FALSE)
rm(temp, nospecies, species_counts)

## Checks
## Number of species in list vs filtered iucn_species db
paste0("# IUCN species in list = ", length(iucn_species))
paste0("# species in gbif.iucn_species = ", 
       dbGetQuery(con, "SELECT COUNT(DISTINCT species) 
                  FROM gbif.iucn_species;
                  "))



## Name matching in GBIF...







