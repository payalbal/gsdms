##'@title filter_iucn_data.R
##'
##'##'@description Filters IUCN spatial data
##'
##'@param iucn.downloaded.data IUCN spatial data for TERRETRIAL MAMMALS. https://www.iucnredlist.org/resources/spatial-data-download
##'@param output_folder (string) A folder to save the outputs to.
##'                              If none specified, no file is written.
##'@param output_name (string) A name to use in saving the outputs.
##'                            Default: 'filtered_iucn'
##'@param species_list (boolean) If TRUE, function returns a 1D array of 
##'                             species names from filtered dataset. If FALSE, 
##'                             functions returns entire filtered dataset, 
##'                             as per filter_fields. Default is TRUE
##'@param filter_fields (string) IUCN data fields that will be retained for 
##'                              the purpose of filtering.
##'@param filter_code (string) IUCN categories for threatened status. 
##'                            Possible values: VU, EN, CR, EW, EX, DD, NT, LC
##'                            Default: "VU", "EN", "CR".
##'@param subspecies (boolean) To include/exclude subspecies from filtered data.
##'                            Default: FALSE to exclude
##'@param marine (boolean) To include/exclude species classified as marine. 
##'                        Default: FALSE to exclude
##'@param freshwater (boolean) To include/exclude species classified as freshwater
##'                            Default: FALSE to exclude
##'@param filter_origin (string) Filter species based on 'origin' information 
##'                              in data i.e. 1-native,2-reintroduced,3-introduced,
##'                              4-vagrant,5-origin uncertain,6-assisted colonisation
##                               Default: native
##'@param extinct (boolean) Filter species based in 'presence' information in data
##'                         i.e. 1-extant,2-possibly extant discontinued,
##'                         3-possibly extant,4-possibly extinct,5-extinct post1500,
##'                         6-presence uncertain
##'                         Default: possibly extinct,extinct post1500,
##'                         presence uncertain
##'@param verbose (boolean) Print messages to console. Default TRUE.
##'
##'@return .txt
##'

filter_iucn_data = function(iucn.downloaded.data,
                           output_folder = NULL,
                           output_name = "filtered_iucn",
                           species_list = TRUE,
                           filter_fields = c("id_no", "binomial", "presence", "origin", 
                                             "seasonal", "year","kingdom","phylum",
                                             "class", "order_","family","genus", "code",
                                             "marine","terrestial","freshwater",
                                             "shape_Leng", "shape_Area"),
                           filter_code = c("VU", "EN", "CR"),
                           subspecies = FALSE,
                           marine = FALSE,
                           freshwater = FALSE,
                           filter_origin = c("native"),
                           extinct = FALSE,
                           verbose = TRUE)
{
  
  ## Read in the data
  dat <- iucn.downloaded.data
  
  ## Catch the original number of records 
  n.rec.start <- nrow(dat)
  
  ## Drop unwanted columns
  dat <- dat[, c(filter_fields)]
  
  ## Filter data threatened status: include VU, EN, CR; exclude EW, EX, DD, NT, LC
  dat <- dat[which(dat$code %in% filter_code),]
  
  ## Remove subspecies, i.e. records with a third name
  if (subspecies == FALSE) {
    dat <- dat[which(sapply(strsplit(as.character(dat$binomial), " "), length) == 2), ]
  }
  
  ## Remove species classified as marine and freshwater
  if (marine == FALSE) {
    dat<- dat[which(dat$marine == "f"), ]
  }
  
  if (freshwater == FALSE) {
    dat <- dat[which(dat$freshwater == "f"),]
  }
  
  # Filter by origin
  if(!is.null(filter_origin)){
    filter_origin[filter_origin == "native"] <- 1
    filter_origin[filter_origin == "reintroduced"] <- 2
    filter_origin[filter_origin == "introduced"] <- 3
    filter_origin[filter_origin == "vagrant"] <- 4
    filter_origin[filter_origin == "originuncertain"] <- 5
    filter_origin[filter_origin == "assistedcolonisation"] <- 6
    dat <- dat[which(dat$origin %in% filter_origin),]
  } else {
    dat <- dat
  }
  
  if (extinct == FALSE) {
    ## Remove species classified as 'extinct'
    dat <- dat[which(dat$presence == 4 | dat$presence == 5 | dat$presence == 6),]
  }
  
  ## Output
  if(!is.null(output_folder)) {
    if(!dir.exists(output_folder))
    {
      warning("output folder not provided")
    }
    output_path <- file.path(output_folder,paste0(Sys.Date(), "_", output_name,".txt"))
    
    if (species_list == TRUE) {
      outdata <- as.data.frame(unique(dat$binomial))
      write.table(outdata, output_path,  row.names = FALSE, col.names = FALSE)
      
    } else {
      outdata <- dat
      write.table(dat, output_path,  row.names = FALSE, col.names = FALSE)
    }
  }
  
  # write some feedback to the terminal
  if(verbose)
  {
    msg1 = paste('Returned object is a text file of dimensions =', dim(outdata)[1], 
                 "rows and ", dim(outdata)[2], "columns")
    msg2 = paste('These data have been also been written to ', output_path)
    msg3 = paste("# records in raw data = ", n.rec.start)
    msg4 = paste("# records in filtered data = ", dim(dat)[1])
    msg5 = paste("# records removed (including spatial duplicates) =", n.rec.start-dim(dat)[1])
    msg6 = paste("# unique species in filtered data =", length(unique(dat$binomial)))
    cat(paste(msg1, msg2, msg3, msg4, msg5, msg6, sep = '\n'))
  } # end if(verbose)
  
  return(outdata)
}

