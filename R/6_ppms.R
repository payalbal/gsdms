## Fitting PPMS

library(ppmlasso)

## Load design martix as per MBH desing/random sample
Xsamp <- readRDS("./output/designmatrix.rds")

## Load filtered biodiversity data
gbif <- read.table("/Volumes/payal_umelb/data/processed/biodiv/2019-04-12_gbif_iucnsp.csv", sep = ",", header = TRUE)
length(unique(gbif$species))