## Extract layers from Hurtt dataset - future layers: LUH2 v2f Release (https://luh.umd.edu/)
## Layers extracted by SSP, land use class and years


import xarray as xr
import itertools
import re
from pathlib import Path
import numpy as np
from numpy import dstack

## Specify extraction parameters
luvars = ['primf', 'primn', 'secdf', 'secdn', 'urban', 'c3ann', 'c4ann', 'c3per', 'c4per', 'c3nfx', 'pastr', 'range']
time_steps = ['0', '5', '10', '15', '20', '25', '30', '35', '40', '45', '50', '55'] ## python indices start from 0


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
  outfile = [outdir, re.sub("_2015-2100", "", Path(infile).stem), "_", lu_var, "_t", time_slice, ".nc"]
  outfile = "".join(outfile)
  extrsub.to_netcdf(outfile)



## Run extract function (loop through extraction parameters for each input file)
outdir = '/tempdata/research-cifs/6300-payalb/uom_data/gsdms_data/outputs/hurtt_future/'

infile = '/tempdata/research-cifs/6300-payalb/uom_data/gsdms_data/landuse/hurtt/ssp1-rcp26_2015-2100.nc'
for i in luvars:
  for j in time_steps:
    extract_nc(infile = infile, outdir = outdir, lu_var = i, time_slice = j)


infile = '/tempdata/research-cifs/6300-payalb/uom_data/gsdms_data/landuse/hurtt/ssp3-rcp70_2015-2100.nc'
for i in luvars:
  for j in time_steps:
    extract_nc(infile = infile, outdir = outdir, lu_var = i, time_slice = j)


infile = '/tempdata/research-cifs/6300-payalb/uom_data/gsdms_data/landuse/hurtt/ssp5-rcp85_2015-2100.nc'
for i in luvars:
  for j in time_steps:
    extract_nc(infile = infile, outdir = outdir, lu_var = i, time_slice = j)



  
  

# ## Debug run - select statements
# import xarray as xr
# 
# ## Specify extraction parameters
# luvars = ['primf', 'urban']
# time_steps = ['1']
# 
# ## >> Load file
# fp = '/tempdata/research-cifs/6300-payalb/uom_data/gsdms_data/landuse/hurtt/ssp5-rcp85_2015-2100.nc'
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
# extrsub.to_netcdf('/tempdata/research-cifs/uom_data/gsdms_data/outputs/hurtt_future/primf_t2015.nc')

