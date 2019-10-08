
## Load libraries
pacman::p_load(tools, DBI, RPostgreSQL)


## Set path to data folder
data_path <- "/Volumes/discovery_data/gsdms_data/gbif" 
# data_path <- "data" # for boab

data_raw <- file.path(data_path, "raw")
data_clean <- file.path(data_path, "processed")

if(!dir.exists(data_raw)){dir.create(data_raw)}
if(!dir.exists(data_clean)){dir.create(data_clean)}


## GBIF backbone taxonomy
##  Ref: GBIF Secretariat (2017). GBIF Backbone Taxonomy. Checklist dataset https://doi.org/10.15468/39omei accessed via GBIF.org on 2019-08-26.
system(paste0("curl http://rs.gbif.org/datasets/backbone/backbone-current.zip -o ", 
              file.path(data_path, "gbif_taxonomy.zip")))
unzip(file.path(data_path, "gbif_taxonomy.zip"), list = TRUE)
unzip(file.path(data_path, "gbif_taxonomy.zip"), files = 'Taxon.tsv', 
      exdir = "./data")
backbone <-fread(file.path(data_path, "Taxon.tsv"))


## GBIF occurrence data
## Source (needs login): https://www.gbif.org/occurrence/download
## - Sign in on gbif.org and click: Get Data > Occurrences > Download. 
## - You can also select filters before downloading, e.g. by year/region/tax
## - Click to download the csv
## - Youll get an email saying your download is ready, with a link that looks 
##   something like this: http://api.gbif.org/v1/occurrence/download/request/<some_number>.zip
## - From a terminal on the server you want it on, do something like:
##   curl http://api.gbif.org/v1/occurrence/download/request/<some_number>.zip  -o all_of_gbif.zip

system(paste0("curl http://api.gbif.org/v1/occurrence/download/request/<some_number>.zip -o all_of_gbif.zip"))
unzip(file.path(data_path, "all_of_gbif.zip"), list = TRUE)
unzip(file.path(data_path, "all_of_gbif.zip"), exdir = data_path)




##  Here, records for Aves, Mammalis and Reptilia with coordinates are downloaded
taxa <- c("gbif_aves", "gbif_mammals", "gnif_reptiles")
gbif_aves <- "http://api.gbif.org/v1/occurrence/download/request/0008696-190813142620410.zip"
gbif_mammals <- "http://api.gbif.org/v1/occurrence/download/request/0008700-190813142620410.zip"
gbif_reptiles <- "http://api.gbif.org/v1/occurrence/download/request/0008702-190813142620410.zip"
urls <- c(gbif_aves, gbif_mammals, gbif_reptiles)

for(i in 1:length(urls)) {
  
  system(paste0("curl ", urls[i], " -o ", paste0(data_path, "/", taxa[i], ".zip")))
  temp <- unzip(paste0(data_path, "/", taxa[i], ".zip"), list = TRUE)
  unzip(paste0(data_path, "/", taxa[i], ".zip"), exdir = data_path)
  file.rename(file.path(data_path, temp$Name), 
              sub(tools::file_path_sans_ext(basename(temp$Name)),
                  taxa[i], file.path(data_path, temp$Name)))
  file.remove(paste0(data_path, "/", taxa[i], ".zip"))
  
}







## ------------------------------------------------------------------------
## Notes:

# there are other ways of downloading files from the terminal, and curl might not be installed.
# 
# Currently, the entire csv is 471GB, but 69GB zipped. It says Estimated 
# download times: 100 Mbit connection : 01:32 hours
# Which is about how long it took for me.