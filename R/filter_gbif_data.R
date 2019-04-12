
##'@title filter_gbif_data.R
##'
##'@description Processing for GBIF data from processed gbif database on server.
##'
##'@param gbif.download.data (data.table) A data.table holding the downloaded GBIF data. GBIF.org 15th October 2018 GBIF Occurrence Download https://doi.org/10.15468/dl.g2zaxo
##'@param gbif.nub.taxonomy (data.table) A data.table holding the downloaded GBIF backbone taxonomy. 5 February 2018 https://www.gbif.org/dataset/d7dddbf4-2cf0-4f39-9b2a-bb099caae36c#description
##'@param subset.gbifnubtaxonomy.byclass (string) Optional - Character string to subset gbif.nub.taxonomy by class. Default: NULL
##'@param output_folder (string) A folder to save the outputs to. If none specified, no file is written.
##'@param output_name (string) A name to use in saving the outputs. Default: 'filtered_data'.
##'@param domain.mask (raster layer) A raster layer specifying the analysis domain
##'@param start.year (integer) The earliest year for which data will be retained  Default: 1970.
##'@param end.year (integer) The latest year for which data will be retained  Default: 2018
##'@param remove_duplicates (boolean) To include or exclude spatial duplicates. Deafult TRUE excludes. 
##'@param spatial.uncertainty.m (float) Distance (m) threshold applied to field coordinateuncertaintyinmeters if provided for inclusion of records. Default: 100000 mts.
##'@param filter_fields (string) GBIF data fields, i.e. columns, that will be reatined for the purpose of filtering
##'@param filter_basisofrecord (string)
##'@param issue_geospatial (string) Geospatial issues to exclude from the output dataset. Default: "ZERO_COORDINATE", "COORDINATE_INVALID",
##' "COORDINATE_OUT_OF_RANGE", "COUNTRY_COORDINATE_MISMATCH", "COORDINATE_REPROJECTION_FAILED", "COORDINATE_REPROJECTION_SUSPICIOUS", 
##' "GEODETIC_DATUM_INVALID" 
##'@param issue_taxonomic (string) Taxonomic issues to exclude from output dataset. Default: "TAXON_MATCH_FUZZY", "TAXON_MATCH_HIGHERRANK", 
##'"TAXON_MATCH_NONE"
##'@param issue_basisofrecord (string) Issues related to 'basisofrecord' to exclude from output dataset. Default: 
##'@param issue_date (string) Issues related to 'eventDate' to exclude from output dataset. Deafult: NULL. If required, use MODIFIED_DATE_UNLIKELY and ignore others
##'@param select_fields (string) List of data fields to be returned in the output dataset. Default: "gbifid", "species", "decimallatitude", 
##'"decimallongitude", "taxonkey"
##'@param verbose (boolean) Print messages to console. Default TRUE.
##'
##'@return data.table
##'
##'@examples output = filter_gbif_data(My.GBIF.data, output.folder = 'C:/Users/data_processed', output.name = My.filtered.GBIF.data, domain.mask = Aust.ras)
##'
##'@export
##'
##'@note description of gbif issues can be found here: https://gbif.github.io/gbif-api/apidocs/org/gbif/api/vocabulary/OccurrenceIssue.html
##'@note eventDate is not required for finding duplicates because MaxEnt only uses one point per lat-long for a species
##'@note code adapted from https://github.com/cwarecsiro/gdmEngine/blob/master/gdmEngine/R/filter_ALA_data.R.


filter_gbif_data = function (gbif.downloaded.data,
                             gbif.nub.taxonomy,
                             subset.gbifnubtaxonomy.byclass = NULL,
                             output_folder = NULL,
                             output_name = "filtereted_gbif",
                             domain.mask = NULL,
                             start.year = 1950,
                             end.year = 2018,
                             remove_duplicates = TRUE,
                             spatial.uncertainty.m = 100000,
                             filter_fields = c("gbifid", "species", "scientificname", 
                                               "countrycode", "decimallatitude", 
                                               "decimallongitude", "coordinateuncertaintyinmeters",
                                               "year","taxonkey","phylum","class","order","family",
                                               "genus","specieskey","basisofrecord","issue"),
                             filter_basisofrecord = c("HUMAN_OBSERVATION", "PRESERVED_SPECIMEN", "OBSERVATION", 
                                                      "MATERIAL_SAMPLE", "MACHINE_OBSERVATION"),
                             issue_geospatial = c("ZERO_COORDINATE", "COORDINATE_INVALID", "COORDINATE_OUT_OF_RANGE", 
                                                  "COUNTRY_COORDINATE_MISMATCH", "COORDINATE_REPROJECTION_FAILED", 
                                                  "COORDINATE_REPROJECTION_SUSPICIOUS", "GEODETIC_DATUM_INVALID"),
                             issue_taxonomic = c("TAXON_MATCH_FUZZY", "TAXON_MATCH_HIGHERRANK", "TAXON_MATCH_NONE"),
                             issue_basisofrecord = c("BASIS_OF_RECORD_INVALID"),
                             issue_date = NULL,
                             select_fields = NULL,
                             verbose = TRUE)
{
  
  ## Load required packages
  if("pacman"%in%installed.packages()){
    library(pacman)
    p_load("data.table","sp","raster")
  } else {
    install.packages("pacman")
    library(pacman)
    p_load("data.table","sp","raster")
  }
  
  ## Read in the data
  dat <- gbif.downloaded.data
  
  ## catch the original number of records 
  n.rec.start <- nrow(dat)
  
  ## Rename columns
  names(dat)[names(dat) == "recyear"] = "year"
  names(dat)[names(dat) == "taxclass"] = "class"
  names(dat)[names(dat) == "taxorder"] = "order"
  names(dat)[names(dat) == "taxfamily"] = "family"

  ## Drop unwanted columns
  names(dat) <- sapply(names(dat), tolower)
  dat <- dat[, filter_fields, with = FALSE]
  
  ## Checks - already removed in first filter ------------------------------- ##
  ## Delete records without spatial coordinates
  dat <- na.omit(dat, cols = c("decimallatitude", "decimallongitude"))
  ## Delete records without species name
  dat <- na.omit(dat, cols = "species")
  ## Filter by date range 
  dat <- dat[year >= start.year & year <= end.year]
  ## Filter by basis of record
  dat <- dat[basisofrecord %in% filter_basisofrecord]
  ## ------------------------------------------------------------------------ ##
  
  ## Filter by issues in data
  dat <- dat[!grep(paste0(c(issue_geospatial,issue_taxonomic, issue_basisofrecord, issue_date), collapse = "|"), dat$issue, perl = TRUE, value = FALSE)]
  
  ## Filter by coordinate uncertainty (if 'coordinateuncertaintyinmeters' field is provided)
  if("coordinateuncertaintyinmeters" %in% filter_fields){
    dat <- dat[!which(coordinateuncertaintyinmeters > spatial.uncertainty.m)]}
  ## note: !which() retains NAs unline which()
  
  ## Filter records by coordinate precision i.e. records with less than 2 decimal places
  dat <- dat[!which(sapply((strsplit(sub('0+$', '', as.character(dat$decimallongitude)), ".", fixed = TRUE)), function(x) nchar(x[2])) < 2)]
  dat <- dat[!which(sapply((strsplit(sub('0+$', '', as.character(dat$decimallatitude)), ".", fixed = TRUE)), function(x) nchar(x[2])) < 2)]
  
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
  
  ## Remove records with incomplete scientific names i.e. only genus or only species provided
  dat <- dat[!which(sapply(strsplit(dat$species, " "), length) < 2)]
  
  ## Remove records with incomplete species names i.e. when species is named as sp., spp. sp
  dat <- dat[!grep(paste0(c("sp.", "spp", "sp"), "$", collapse = "|"), dat$species, perl = TRUE, value = FALSE)]
  
  ## Remove records with scientific names not included in gbif backbone taxonomy
  ## 'canonicalName' is the 'scientificName' without authorship
  gbif.nub.taxonomy <- gbif.nub.taxonomy[, .(taxonID, canonicalName,taxonRank, taxonomicStatus, 
                                             kingdom, phylum, class, order, family, genus)]
  if(!is.null(subset.gbifnubtaxonomy.byclass)){
    gbif.nub.taxonomy <- gbif.nub.taxonomy[class==subset.gbifnubtaxonomy.byclass]
  } else {
    dat_species_list <- unique(dat$species)
    check_list <- dat_species_list[which(!(dat_species_list %in% unique(gbif.nub.taxonomy$canonicalName)), arr.ind = TRUE)]
    if(!identical(check_list, character(0))){
      dat <- dat[!species %in% check_list]
    } # end !identical(check_list, character(0))    ## cannot catch character(0) with is.null
  } # end !is.null(gbif.nub.taxonomy)
  
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
    msg5 = paste("# records removed (including spatial duplicates) =", n.rec.start-dim(dat)[1])
    cat(paste(msg1, msg2, msg3, msg4, msg5, sep = '\n'))
  } # end if(verbose)
  
  return(dat)
  
} # end filter_gbif_data function