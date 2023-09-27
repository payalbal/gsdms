## Download bioclim tiles and stitch to make global layer
## author: Chris Ware


# output file paths to set
# -----------------------------------------------------------------------------
# just somewhere to download each zip archive to (will be deleted once used)
# must have .zip ext
zipdst = './data/bio_30s.zip'

# a dir to extract raster files to
rasterdst = './data/bioclim/'

# somewhere to write the global raster for each bioclim variable to
bioclim_dst = paste0(rasterdst, 'global_bio_') # or whatever

# Looking at the getData function, it finds which tile a given lat lon  falls witin. 
# The tiles are then indexed by row/col which maps to a url. 
# Can cut to the chase and just create urls with all combos of row/cols (i.e. tile ids)
# -----------------------------------------------------------------------------
urls = NULL
for (r in 1:5){
  for (c in 1:12){
    
    rc = paste0(r-1, c-1)
    fn = paste0('bio_', rc, '.zip')
    thisurl = paste0('https://biogeo.ucdavis.edu/data/climate/worldclim/1_4/tiles/cur/', fn)
    urls = c(urls, thisurl)
    
  }
}

# check all the urls exist (which hopefully means all the tiles can be downloaded)
# -----------------------------------------------------------------------------
require(bitops)
require(RCurl)
test = lapply(urls, url.exists)
all(unlist(test))

# Then, it's possible to loop over the urls, downloading, and unzipping them. 
# -----------------------------------------------------------------------------
for (url in urls){
  
  response = tryCatch({
    download.file(url, destfile = zipdst)
  }, 
  error = function(e){e} # shouldn't get called
  )
  
  print(response) # should be 0
  
  unzip(zipdst, exdir = rasterdst)
  
  file.remove(zipdst)
} 
  
# rasterdst should now be full of bioclim 1-19 tiles (60 tiles for each var)
# loop over each bioclim variable, collect all tiles, and mosaic. 
# -----------------------------------------------------------------------------
for (i in 1:19){
  f = list.files(rasterdst, pattern = paste0('bio', i, '_'), full.names = TRUE)
  f = f[grep('.bil$', f)]
  
  f = lapply(f, raster)
  
  # specify a function for any overlapping cells (of which there won't be any, 
  # but the function requires a function...)
  f$fun = mean
  mos = do.call(mosaic, f)
  
  this_dst = paste0(bioclim_dst, i, '.tif') # or whatever
  writeRaster(mos, this_dst)
}

