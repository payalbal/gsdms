
rm(list = ls())
gc()
x <- c('bitops', 'RCurl', 'rstudioapi')
lapply(x, require, character.only = TRUE)
rm(x)

data_dir <- "/tempdata/research-cifs/6300-payalb/uom_data/gsdms_data"

input_files <- file.path(data_dir, "CHELSA", "bio_future_urls.txt")
target_folder <- file.path(data_dir, "CHELSA", "future", "raw")
system(sprintf("wget --no-host-directories --force-directories --input-file=%s -P %s", input_files, target_folder))
list.files(target_folder)