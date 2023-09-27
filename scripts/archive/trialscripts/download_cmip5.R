## Description: Download CMIP5 data from Worldclim and prepare mean/SD climate layers
## Author: Simon Kapitza
library("matrixStats")

startime <- Sys.time()

## Download
download_path <- file.path(".", "temp")
dir.create(download_path)
gcm_path <- file.path(".", "GCM")
dir.create(gcm_path)
layer_path <- file.path(".", "RData")
dir.create(layer_path)

global_mask <- raster("global_mask.tif")
res <- 0.5
scens <- c("26", "85")
# models <- c("GF", "MP", "BC", "CC", "CN", "GS", "HD", "HE", "IP", "MI", "MR", "MC", "MG", "NO")
models <- c("GF")


for(i in 1:length(models)){
  temp_stack <- list()
  for (j in 1:length(scens)){
    f <- tryCatch(getData("CMIP5", rcp = scens[j], year = 70, model = models[i], res = res, var = "bio", path = download_path), error = function (e) NA)
    if(class(f) != "RasterStack"){next}
    temp_stack[[j]] <- crop(f, global_mask)
  }
  save(temp_stack, file = paste0(download_path, models[i], ".rda"))
  unlink(paste0(download_path, "/cmip5"), recursive = TRUE)
  rm(temp_stack)
}

Sys.time() - startime ## ~21 mins

load("./GCM/tempGF.rda")

## Calculate mean and mean +-SD
startime <- Sys.time()

gcm <- list.files(gcm_path, full.names = T)

r <-global_mask
inds <- which(!is.na(r[]))
quartiles <- c("q1", "q2", "q3")

k <- 1
for(k in 1:length(scens)){
  saveRDS(stack(), file = paste0(layer_path, "/bio", "q1_", scens[k], ".rds"))
  saveRDS(stack(), file = paste0(layer_path, "/bio", "q2_", scens[k], ".rds"))
  saveRDS(stack(), file = paste0(layer_path, "/bio", "q3_", scens[k], ".rds"))
  print(paste0("processing rcp", scens[k]))
  for(j in 1:19){
    print(paste0("processing cov: ", j))
    bio <- stack()
    for(i in 1:length(models)){
      print(paste0("processing model: ", models[i]))
      dat <- get(load(gcm[[i]]))[[k]]
      bio <- stack(bio, dat[[j]])
    }
    
    writeRaster(bio, filename = paste0("download_path", "bio", ".tif"), format="GTiff", overwrite = TRUE, bylayer = TRUE)
    rm(bio)
    
    print(paste0("getting quantiles..."))
    df1 <- na.omit(as.matrix(getValues(bio)))
    c <-rowQuantiles(df1, probs = c(0.25, 0.5, 0.75))
    for(m in 1:length(quartiles)){
      bioclim <- readRDS(file = paste0(layer_path, "/bio", quartiles[m], "_", scens[k], ".rds"))
      r[inds] <- c[,m]
      names(r) <- paste0("bio", j)
      saveRDS(stack(bioclim, r), file = paste0(layer_path, "/bio", quartiles[m], "_", scens[k], ".rds"))
    }
  }
}

Sys.time() - startime


## doparallel for getValues()
require(foreach,doParallel,raster,rgdal)
stack_in <- bio
#Determine optimal block size for loading in MODIS stack data
block_width = 15
nrows = dim(stack_in)[1]
nblocks <- nrows%/%block_width
bs_rows <- seq(1,nblocks*block_width+1,block_width)
bs_nrows <- rbind(matrix(block_width,length(bs_rows)-1,1),nrows-bs_rows[length(bs_rows)]+1)
print('Working on the following rows')
print(paste(bs_rows))

#Register the parallel backend
library("doParallel", "iterators", "parallel", "foreach")
library("raster")
cl <- makeCluster(30) # or detectCores() to use the default number of cores available
registerDoParallel(cl)

result <- foreach(i = 1:length(bs_rows), .combine = rbind, .packages = c("raster")) %dopar% {
  print(paste("Working on row",i))
  stack_values = getValues(stack_in, bs_rows[i], bs_nrows[i])
  return(FUN(stack_values))
}

stopImplicitCluster()
