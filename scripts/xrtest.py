import xarray as xr

fp = '/tempdata/research-cifs/uom_data/gsdms_data/landuse/Hurtt/ssp5-rcp85_2015-2100.nc'
sspnc = xr.open_dataset(fp, decode_times=False)

print( sspnc )

print('PRIMEF')
print( sspnc['primf'] )

print('select time 80')
extr = sspnc.sel( time = 80.0, method='nearest', tolerance = 0.5, drop=True )
print( extr['primf'] )

print('subset lats and longs')
extrsub = extr.sel( lon = slice( 130, 140 ), lat = slice( -34, -36 ) )
print( extrsub['primf'] )

extr.to_netcdf('/tempdata/research-cifs/uom_data/gsdms_data/landuse/Hurtt/newtest.nc')
