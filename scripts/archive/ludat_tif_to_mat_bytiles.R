library(sp)
library(raster)

# data_covs <- "./data" # server
data_covs <- "/Volumes/discovery_data/gsdms_data" # local
files <- list.files(file.path(data_covs, "copernicus", "global_frac", "sklu_classes", "tifs"), full.names = TRUE)
sklu_classes <- c("urban", "crop", "forest", "grass", "other")


for(n in 1:length(files)){
  
  lu_stack <- stack(files[n])
  
  lu <- matrix(data = NA, ncell(lu_stack), nlayers(lu_stack)) ## replaced ncell(reg_mask)
  for(i in 1: nlayers(lu_stack)){
    lu[,i] <- getValues(lu_stack[[i]])
  }
  lu <- na.omit(lu)
  colnames(lu) <- sklu_classes

  ## Checks  
  print(paste0("Checking file: ", basename(files[n]), " ...."))
  print("Total cells - NAs in stack = dim of matrix: ")
  print(ncell(lu_stack[[1]]) - sum(is.na(getValues(lu_stack[[1]]))) == nrow(lu))
  print( " Values in each row add up to 1: ")
  rs <- round(rowSums(lu, na.rm = TRUE, dims = 1))
  print(all(rs==1))
  ## Note some minute residual decimal places, therefore rounding to 1
  
  ## Saving
  print(paste0("Saving file: ", basename(files[n]), " ...."))
  saveRDS(lu, file = sub("tifs" , "mats",  paste0(tools::file_path_sans_ext(files[n]), ".rds")))
  
}
