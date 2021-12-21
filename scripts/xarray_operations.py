import xarray as xr
import itertools
import re
from pathlib import Path


## Specify extract function
def extract_nc(infile, outdir, lu_var, time_slice):
  
  ## >> Load file
  sspnc = xr.open_dataset(infile, decode_times=False)

  ## >> Specify extent
  extr = sspnc.sel( lon = slice( -180, 180 ), lat = slice( 90, -90 ) )

  ## >> Select variable: primf
  extr = sspnc[lu_var]

  ## >> Select time slices
  extrsub = extr.sel( time = time_slice, method='nearest', tolerance = 0.5, drop=True )

  ## >> Write output file
  outfile = [outdir, re.sub("_2015-2100", "", Path(infile).stem), "_", lu_var, "_t", time_slice, ".nc"]
  outfile = "".join(outfile)
  extrsub.to_netcdf(outfile)



## Call function
outdir = '/tempdata/research-cifs/uom_data/gsdms_data/outputs/hurtt/'

luvars = ['primf', 'primn', 'secdf', 'secdn', 'urban', 'c3ann', 'c4ann', 'c3per', 'c4per', 'c3nfx', 'pastr', 'range']
time_steps = ['1', '6', '11', '16', '21', '26', '31', '36', '41', '46', '51', '56']
# luvars = ['primf', 'urban']
# time_steps = ['1']

## >> Run loops for each SSP
infile = '/tempdata/research-cifs/uom_data/gsdms_data/landuse/Hurtt/ssp1-rcp26_2015-2100.nc'
for i in luvars:
  for j in time_steps:
    extract_nc(infile = infile, outdir = outdir, lu_var = i, time_slice = j)
    

infile = '/tempdata/research-cifs/uom_data/gsdms_data/landuse/Hurtt/ssp3-rcp70_2015-2100.nc'
for i in luvars:
  for j in time_steps:
    extract_nc(infile = infile, outdir = outdir, lu_var = i, time_slice = j)
    
    
infile = '/tempdata/research-cifs/uom_data/gsdms_data/landuse/Hurtt/ssp5-rcp85_2015-2100.nc'
for i in luvars:
  for j in time_steps:
    extract_nc(infile = infile, outdir = outdir, lu_var = i, time_slice = j)
  
    

  





# ## Example run - detailed
# import xarray as xr
# 
# ## >> Load file
# fp = '/tempdata/research-cifs/uom_data/gsdms_data/landuse/Hurtt/ssp5-rcp85_2015-2100.nc'
# sspnc = xr.open_dataset(fp, decode_times=False)
# print( sspnc )
# 
# ## >> Specify extent
# print('subset extent')
# extr = sspnc.sel( lon = slice( -180, 180 ), lat = slice( 90, -90 ) )
# 
# ## >> Select variable: primf
# print('select var primf')
# extr = sspnc['primf']
# print( extr )
# 
# ## >> Select time slices
# print('select time step 1')
# extrsub = extr.sel( time = 1.0, method='nearest', tolerance = 0.5, drop=True )
# print( extrsub )
# 
# extrsub.to_netcdf('/tempdata/research-cifs/uom_data/gsdms_data/outputs/hurtt/primf_t2015.nc')

