 

## -------------------------- EXPLORING GBIF DATA -------------------------- ##

library(rgbif)
# https://ropensci.org/tutorials/rgbif_tutorial/



## ----------------------------------
## GBIF fields:
## ----------------------------------

# https://www.gbif.org/developer/occurrence
# http://rs.gbif.org/core/dwc_occurrence.xml
# ?occ_data

temp <- occ_data(hasCoordinate = TRUE, start=0, limit=10)  
# temp <- occ_search(hasCoordinate = TRUE, fields = "all", return = "all", start = 0, limit =10)
# occ_data for better speed, returns same $data and $meta as occ_search but wothout $hier and $media
names(temp)

fields1 <- sort(names(temp$data)) # field names using occ_data in rgbif
fields2 <- sort(names(read.csv("R/gbifcsv_header.csv", sep="\t")))  # field names in csv, provided by Casey

fields1 <- tolower(fields1) # convert character elements to lower case 
fields2 <- tolower(fields2)
setdiff(fields1, fields2) # field names in occ_data BUT NOT in csv
length(setdiff(fields1, fields2)) # 1 additional field counted i.e. 'issues' as opposed ot 'issue' in csv
setdiff(fields2, fields1) # field names in csv BUT NOT in occ_data
length(setdiff(fields1, fields2)) # 1 additional field counted i.e. 'issue' as opposed ot 'issues' in occ_data
#write.table(fields1, "inprogress/explout/rgbif_header.txt", sep=",", col.names = F, row.names = F, quote = FALSE)
#write.table(fields2, "inprogress/explout/gbif_csv_header.txt", sep=",", col.names = F, row.names = F, quote = FALSE)



## ----------------------------------
## Counting records using occ_count
## ----------------------------------

occ_count(georeferenced=TRUE) # = 889554826 [best]
# same as: occ_data(hasCoordinate = TRUE)$meta$count [better than occ_search]
# same as: occ_search(hasCoordinate = TRUE, fields = "all", return = "meta")$count [worst]
occ_data(hasCoordinate = TRUE, hasGeospatialIssue = FALSE)$meta$count # = 889554805

system.time(occ_count(georeferenced=TRUE))
system.time(occ_data(hasCoordinate = TRUE)$meta$count)
# help: https://www.r-bloggers.com/5-ways-to-measure-running-time-of-r-code/

occ_count(georeferenced=TRUE) # does not work with year
occ_count(type = 'basisOfRecord')
occ_count(georeferenced=TRUE, basisOfRecord = 'OBSERVATION') # combination not supported in occ_count
occ_count(type='year', from=1900, to=2017) # specified range does not work. Default range: 1500-2019
by.year <- unlist(occ_count(type = "year"))
by.year.georef <- unlist(occ_count(georeferenced=TRUE, type = "year"))
which((by.year == by.year.georef) == FALSE)  # gives same number of records; combination with georef not supported in occ_count
sum((by.year > 200000) == TRUE, na.rm = TRUE) # years with >200 000 records
which((by.year > 200000) == TRUE)

by.country <- unlist(occ_count(type = 'countries'))
by.country.georef <- unlist(occ_count(type = 'countries', georeferenced=TRUE))
which((by.country == by.country.georef) == FALSE) # gives same number of records; combination with georef not supported in occ_count
sum((by.country > 200000) == TRUE, na.rm = TRUE)
which((by.country > 200000) == TRUE)

# NOTES:
# occ_count: occ_count is not very versitile to combine arguments, alternative is to use occ_search/occ_data



## ----------------------------------
## Counting records using occ_search/occ_data
## ----------------------------------

# ISSUE: Takes much longer to run compared to occ_count

library(tidyverse) # includes dplyr so no need to load that separately
library(foreach)

# counts using geographic filters
occ_data(hasCoordinate = TRUE, hasGeospatialIssue = FALSE)$meta$count # = 889554805
# counts by other filters - basis of record (note: can only provide one value for basisofrecord at a time)
occ_data(hasCoordinate = TRUE, hasGeospatialIssue = FALSE, basisOfRecord = c("OBSERVATION"))$meta$count

# counts by year
yr.start <- 2016 # as trial; 1950 for full data download 
yr.end <- 2017

by.year.filtered <- foreach(i = c(yr.start:yr.end)) %do% {
  occ_data(hasCoordinate = TRUE, hasGeospatialIssue = FALSE, year = i)$meta$count
}
names(by.year.filtered) <- c(yr.start:yr.end)

# counts by countries
countries <- isocodes$code[1:2]

by.country.filtered <- foreach(j=countries, combine = rbind) %do% {
    occ_data(hasCoordinate = TRUE, hasGeospatialIssue = FALSE, country = j)$meta$count
}
names(by.country.filtered) <- countries

# counts by year & country
by.year.country <- foreach(i = c(yr.start:yr.end), combine = rbind) %do% {
  foreach(j=countries, combine = rbind) %do% {
    occ_data(hasCoordinate = TRUE, hasGeospatialIssue = FALSE, year = i, country = j)$meta$count
  }
}
names(by.year.country) <- c(yr.start:yr.end)
for(k in 1:length(yr.start:yr.end)){names(by.year.country[[k]]) <- countries}
for(k in yr.start:yr.end){assign(paste(k,"_list",""), unlist(by.year.country[[1]], use.names = TRUE))}



## ------------------------------------------
## Exploring taxonomic info in gbif backbone
## ------------------------------------------

library(data.table)
# https://github.com/Rdatatable/data.table/wiki/Getting-started

gbifdat.dir <- '/Volumes/payal_umelb/data/biodiversity_data/GBIF/'
backbone.raw <-fread(paste(gbifdat.dir,"backbone-current/Taxon.tsv", sep =""))
# backbone.raw <-as.data.frame(fread(paste(gbifdat.dir,"backbone-current/Taxon.tsv", sep =""))) # read in data as a data frame instead
backbone <- backbone.raw[, c("taxonID", "scientificName", "canonicalName","taxonRank", "taxonomicStatus", "kingdom", "phylum", "class","order", "family", "genus")]
# unique(is.na(backbone_sub)) # to check if there are NAs, takes some time to run

sum(backbone$kingdom == "Animalia") # to count number of rows for Animalia, but this doesn't tell us much
sum(backbone$kingdom =="Animalia" & backbone$taxonRank == "species") # to count number of species in Animalia in the backbone
# note: this is more than the number of species listed on gbif webpage for Animalia (https://www.gbif.org/species/1/metrics?root=true) perhaps because not all species will have recorded occurrences in the database.
sum(backbone$kingdom =="Animalia" & backbone$phylum == "Chordata" & backbone$taxonRank == "species")


length(backbone[backbone$kingdom =="Animalia" & backbone$taxonRank == "species", ]$canonicalName)
length(unique(backbone[backbone$kingdom =="Animalia" & backbone$taxonRank == "species", ]$canonicalName))
backbone.splist <- unique(backbone[backbone$kingdom =="Animalia" & backbone$taxonRank == "species", ]$canonicalName) 
## to list all species for Animalia listed in gbif backbone
## 'canonicalName' is the 'scientificName' without authorship

# Taxon IDs and subsetting based on taxrank()
kingdom <- backbone[taxonRank == "kingdom"]
# "%w/o%" <- function(x, y) x[!x %in% y] #--  x without y
# kingdom[, c(taxrank() %w/o% c("kingdom", "species", "infraspecific")) := NULL]
kingdom[, .(taxonID, kingdom)]

phylum <- backbone[kingdom == "Animalia" & taxonRank == "phylum"]
phylum[, .(taxonID, kingdom, phylum)]

unique(backbone$kingdom)
aves.N <- dim(backbone[which(backbone$class == "Aves")])
aves.sp <- unique(backbone[which(backbone$class == "Aves")]$canonicalName)

carnivora.N <- dim(backbone[which(backbone$order == "Carnivora")])
carnivora.sp <- unique(backbone[which(backbone$order == "Carnivora")]$canonicalName)
carnivora.fam <- unique(backbone[which(backbone$order == "Carnivora")]$family)
  
## ----------------------------------
## Counting by taxonomic info
## ----------------------------------

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



## ----------------------------------
## Species list from GBIF backbone csv using dplyr
## ----------------------------------

# Example for family Felidae using download key from gbif website
# https://poldham.github.io/abs/gbif.html#reviewing_gbif_data
felidae_csv <- fread(paste(gbifdat.dir,"Temp/0002279-180412121330197.csv", sep =""), na.strings = c("", NA))
unique(felidae_csv$species)
species_names <- felidae_csv %>% filter(taxrank == "SPECIES") %>% distinct(species, .keep_all = FALSE)


# Other ways to access species lists
# https://recology.info/2012/10/rgbif-newfxns/



## ----------------------------------
## Connecting to localserver via R
## ----------------------------------

library(DBI)
library(RPostgreSQL)

# drv <- dbDriver('PostgreSQL')  
# db <- 'gbifrecords'  
# host_db <- 'localhost'  
# db_port <- '5432'  
# db_user <- 'payal'  
# db_password <- ''
# conn <- dbConnect(drv, dbname=db, host=host_db, port=db_port, user=db_user, password=db_password)

rv <- dbDriver("PostgreSQL")  #Specify a driver for postgreSQL type database
con <- dbConnect(drv, dbname="qaeco_spatial", user="qaeco", password="Qpostgres15", host="boab.qaeco.com", port="5432")  #Connection to database server on Boab


## ----------------------------------
## Counts relevant to gbif_firstfilter
## ----------------------------------
# Total gbif records >> 982293333
occ_count()

# Recordss for Kingdom Animalia >> 737988477
occ_data(taxonKey = 1)$meta$count

# Records for Animalia with lat-long coordinates >> 701405414
occ_data(taxonKey = 1, hasCoordinate = TRUE)$meta$count

# Records in Animalia with lat-long coordinates within the years 1950-2017 using >> 652772682
occ_data(taxonKey = 1, hasCoordinate = TRUE, year = '1950,2017')$meta$count 
# check against https://www.gbif.org/occurrence/search?has_coordinate=true&taxon_key=1&year=1950,2017

# Records in the downloaded gbif csv file, i.e. rows in gbif.raw_data on boab
raw.n <- as.numeric(dbGetQuery(con,"SELECT reltuples::bigint AS estimate
                               FROM   pg_class
                               WHERE  oid = 'gbif.raw_data'::regclass;
                               "))
  # count obtained from db metadata 

# Records returened from running the first filter query on gbif data
as.numeric(dbGetQuery(con,"
                      SELECT COUNT(*)
                      FROM gbif.raw_data
                      WHERE kingdom = 'Animalia'
                      AND decimallatitude IS NOT NULL
                      AND decimallongitude IS NOT NULL
                      AND year BETWEEN '1950' AND '2017';
                      "))

# Records in new table, i.e. rows in gbif.firstfilter
firstfilter.n <- as.numeric(dbGetQuery(con,"
                                       SELECT reltuples::bigint AS estimate
                                       FROM   pg_class
                                       WHERE  oid = 'gbif.firstfilter'::regclass;
                                       "))
  # count obtained from db metadata (compare to count above)



# Parallelizing 
# registerDoMC(detectCores() - 1)
# https://discuss.ropensci.org/t/rgbif-occ-search-and-occ-download-questions/406


