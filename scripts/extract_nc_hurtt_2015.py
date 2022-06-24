## Extract layers from Hurtt dataset - historical layers: LUH2 v2h Release (https://luh.umd.edu/)
## Layer extracted years (2015)


import xarray as xr
import itertools
import re
from pathlib import Path
import numpy as np
from numpy import dstack


## Specify extraction parameters
luvars = ['primf', 'primn', 'secdf', 'secdn', 'urban', 'c3ann', 'c4ann', 'c3per', 'c4per', 'c3nfx', 'pastr', 'range']
time_steps = '1165' ## last time step for 2015 i.e. 1166 in R, but python indices start from 0


## Specify extract function
def extract_nc(infile, outdir, lu_var, time_slice):

  ## >> Load file
  sspnc = xr.open_dataset(infile, decode_times=False)

  ## >> Specify extent
  extr = sspnc.sel( lon = slice( -180, 180 ), lat = slice( 90, -60 ) )

  ## >> Select variable: primf
  extr = extr[lu_var]

  ## >> Select time slices
  extrsub = extr.sel( time = time_slice, method='nearest', tolerance = 0.5, drop=True )


  ## >> Write output file
  outfile = [outdir, re.sub("states", "historic", Path(infile).stem), "_", lu_var, "_t", time_slice, ".nc"]
  outfile = "".join(outfile)
  extrsub.to_netcdf(outfile)



## Run extract function
outdir = '/tempdata/research-cifs/6300-payalb/uom_data/gsdms_data/outputs/hurtt_2015/'

infile = '/tempdata/research-cifs/6300-payalb/uom_data/gsdms_data/landuse/hurtt/LUH2/LUH2_v2h/states.nc'
for i in luvars:
  extract_nc(infile = infile, outdir = outdir, lu_var = i, time_slice = time_steps)

  
  
  
  
  
# ## Debug run - select statements
# import xarray as xr
# 
# ## Specify extraction parameters
# luvars = ['primf']
# time_steps = '1165'
#
# ## >> Load file
# fp = '/tempdata/research-cifs/6300-payalb/uom_data/gsdms_data/landuse/hurtt/LUH2/LUH2_v2h/states.nc'
# sspnc = xr.open_dataset(fp, decode_times=False)
# print( sspnc )
# 
# ## >> Specify extent
# print('subset extent')
# extr = sspnc.sel( lon = slice( -180, 180 ), lat = slice( 90, -60 ) )
# 
# ## >> Select variable: primf
# print('select var primf')
# extr = sspnc['primf']
# print( extr )
# 
# ## >> Select time slices
# print('select time step 1')
# extrsub = extr.sel( time = 1165, method='nearest', tolerance = 0.5, drop=True )
# print( extrsub )
#
# extrsub.to_netcdf('/tempdata/research-cifs/6300-payalb/uom_data/gsdms_data/outputs/hurtt_2015/primf_t2015.nc')



# ## Debug run - string join statements
# lu_var = 'primf'
# time_slice = '1165'
# infile = '/tempdata/research-cifs/6300-payalb/uom_data/gsdms_data/landuse/hurtt/LUH2/LUH2_v2h/states.nc'
# outdir = '/tempdata/research-cifs/6300-payalb/uom_data/gsdms_data/outputs/hurtt_2015/'
# outfile = [outdir, re.sub("states", "historic", Path(infile).stem), "_", lu_var, "_t", time_slice, ".nc"]
# outfile = "".join(outfile)
# print(outfile)
