## ------
## TO ADD TO PPM CODE ####

## Select regions by species occurence data - in PPM code
dbSendQuery(con, "DROP TABLE IF EXISTS temp;
                  CREATE TABLE temp AS
                  SELECT t1.*, t2.*
                  FROM gbif_arachnida AS t1, ecoregions_dinerstein_2017 AS t2
                  WHERE ST_Contains(t2.wkb_geometry, t1.points_geom)
                  LIMIT 100;")


## Filter by spatial domain
if(!is.null(domain.mask)){

  ## Filter by extent
  dat <- dat[which(decimallongitude > domain.mask@extent@xmin)]
  dat <- dat[which(decimallongitude < domain.mask@extent@xmax)]
  dat <- dat[which(decimallatitude > domain.mask@extent@ymin)]
  dat <- dat[which(decimallatitude < domain.mask@extent@ymax)]

  ## Filter by location on spatial grid (remove points falling outside of mask)
  sp <- SpatialPoints(dat[,.(decimallongitude,decimallatitude)])
  grd.pts<-extract(domain.mask, sp)
  dat <- dat[!is.na(grd.pts),]

} else {
  warning("domain.mask not provided")
}

## Remove spatial duplicates if TRUE
## identified by species name and coordinates to give only one record for a location
if(remove_duplicates == TRUE){
  dat <- unique(dat, by =c("species", "decimallongitude", "decimallatitude"))
}


## Retain species with >= 20 occurrences
dat <- dat[dat$species %in% names(which(table(dat$species) >= 20)),]

## Create cleaned data file for modelling + log file of how it was created
if(is.null(select_fields)){
  select_fields = c("gbifid", "species", "decimallatitude", "decimallongitude", "taxonkey")
} else {
  check_fields = all(select_fields %in% names(dat))
  if(!check_fields){
    select_fields = c("gbifid", "species", "decimallatitude", "decimallongitude", "taxonkey", "issue")
    warning("Specified select_fields were not found in the dataset - returning default fields instead")
  } else {
    select_fields = c(select_fields)
  }
}
dat <- dat[, select_fields, with = FALSE]

## Write the data to file, if an output folder is specified
if(!is.null(output_folder))
{
  if(!dir.exists(output_folder))
  {
    dir.create(output_folder)
  } # end if !dir.exists

  output_path <- file.path(output_folder,paste0(Sys.Date(), "_", output_name,".csv"))
  write.csv(dat, output_path, row.names=FALSE)

  ## Write a log file describing how the data was created *************************************
  fileConn<-file(file.path(output_folder,paste0(output_name,"_",Sys.Date(),"_log_file.txt")),'w')
  writeLines("#######################################################################",con = fileConn)
  writeLines("###",con = fileConn)
  writeLines("### GBIF data filtration log file ",con = fileConn)
  writeLines("###",con = fileConn)
  writeLines(paste0("### Created ",Sys.time()," using the filter_gbif_data() function."),con = fileConn)
  writeLines("###",con = fileConn)
  writeLines("#######################################################################",con = fileConn)
  writeLines("",con = fileConn)
  writeLines(paste0("Output data file = ", output_path),con = fileConn)
  writeLines(paste0("Domain mask applied = ", domain.mask@file@name),con = fileConn)
  writeLines(paste0("Data restricted to after ", start.year),con = fileConn)
  writeLines(paste0("Data restricted to spatial uncertainty < ",spatial.uncertainty.m, "m"),con = fileConn)
  writeLines(paste0("Spatial duplicates ", if(remove_duplicates == TRUE){"removed"}else{"retained"},con = fileConn))
  writeLines(paste0("Number of records before filtering = ", n.rec.start),con = fileConn)
  writeLines(paste0("Number of records after filtering = ", nrow(dat)),con = fileConn)
  writeLines("#######################################################################",con = fileConn)
  close(fileConn)
  ## *****************************************************************************************
} # end !is.null(output_folder)

# write some feedback to the terminal
if(verbose)
{
  msg1 = 'Returned object is a data.table'
  msg2 = paste('These data have been also been written to ', output_path)
  msg3 = paste("# records in raw data = ", n.rec.start)
  msg4 = paste("# records in filtered data = ", dim(dat)[1])
  msg5 = paste("# records removed =", n.rec.start-dim(dat)[1])
  msg6 = paste0("Spatial duplicates ", if(remove_duplicates == TRUE){"removed"}else{"retained"})
  cat(paste(msg1, msg2, msg3, msg4, msg5, msg6, sep = '\n'))
} # end if(verbose)

return(dat)



## NOT VALID BECAUSE FIELD NOT RETAINED IN DATABASE
## Filter by coordinate uncertainty (if 'coordinateuncertaintyinmeters' field is provided)
if("coordinateuncertaintyinmeters" %in% filter_fields){
  dat <- dat[!which(coordinateuncertaintyinmeters > spatial.uncertainty.m)]}
## note: !which() retains NAs unline which()

## Filter records by coordinate precision i.e. records with less than 2 decimal places
dat <- dat[!which(sapply((strsplit(sub('0+$', '', as.character(dat$decimallongitude)), ".", fixed = TRUE)), function(x) nchar(x[2])) < 2)]
dat <- dat[!which(sapply((strsplit(sub('0+$', '', as.character(dat$decimallatitude)), ".", fixed = TRUE)), function(x) nchar(x[2])) < 2)]