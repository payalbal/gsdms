
x <- c('rgdal', 'geostuff', 'tools', 'bitops', 
       'RCurl', 'gdalUtils', 'usethis',
       'parallel', 'doMC')
lapply(x, require, character.only = TRUE)
rm(x)


mc.cores = future::availableCores()-2
proj.res.km <- 10
data_dir <- "/home/payalb/gsdms_r_vol/tempdata/research-cifs/uom_data/gsdms_data"
output_dir <- file.path(data_dir, "outputs", sprintf("layers_%sk", proj.res.km))

biofuture <- list.files("/home/payalb/gsdms_r_vol/tempdata/research-cifs/proj-2200_nature_futures21-1128.4.411/global/cmip6/bio/mean", full.names = TRUE, recursive = TRUE)


infiles = biofuture
outfiles = file.path(output_dir, paste0(tools::file_path_sans_ext(basename(infiles)), "_clip.", tools::file_ext(infiles)))
new_extent <- "-180 -60 180 90"

parallel::mclapply(seq_along(infiles),
                   function(x) system(paste0("gdalwarp -overwrite -ot Byte",
                                             " -te ", new_extent, " ",
                                             infiles[x], " ", outfiles[x])),
                   mc.cores = mc.cores, mc.preschedule = TRUE)