# LOADING AND EXPLORING PREDICTS DATABASE

# Ref: https://inbo.github.io/tutorials/data-handling-large-files-R.html

setwd("/Volumes/UoM_DATA/OneDrive - The University of Melbourne/biodiversity_data/PREDICTS")
wd <- getwd()  

library("data.table")

# read.table.timing <- system.time(read.table(paste(wd, "/resource.csv", sep=""), header=TRUE, sep=",", na.strings = "NA"))
    # Cannot use read.table perhaps because there are NA entries..?
# readr.timing <- system.time(read_delim(paste(wd, "/resource.csv", sep=""), ",", col_names = TRUE))
# data.table.timing <- system.time(allData <- fread(paste(wd, "/resource.csv", sep=""), showProgress = FALSE))
# dat <- data.frame(method = c('read.table', 'readr', 'fread'), 
#                    timing = c(NA, readr.timing[3], data.table.timing[3]))
# dat


# Loading data
allData <- fread(paste(wd, "/resource.csv", sep=""), showProgress = FALSE)

typeof(allData)
str(allData)

colnames(allData)

unique(allData$Diversity_metric)
unique(allData$Diversity_metric_unit)
unique(allData$Diversity_metric_type)

