## Aggregate crop layers
# https://stackoverflow.com/questions/10129561/how-to-find-max-values-from-multiple-lists

import xarray as xr
import itertools
import re
from pathlib import Path
import numpy as np
from numpy import dstack
import os
from os import listdir
from os.path import isfile, join


# ## >> Specify max function
# def max_nc(scenario, t, infile1, infile2, infile3, infile4, infile5, outdir):
# 
#     ## Read in .nc files
#     A = xr.open_dataset(infile1)
#     B = xr.open_dataset(infile2)
#     C = xr.open_dataset(infile3)
#     D = xr.open_dataset(infile4)
#     E = xr.open_dataset(infile5)
# 
#     ## Convert to numpy arrays
#     A = A.to_array()
#     B = B.to_array()
#     C = C.to_array()
#     D = D.to_array()
#     E = E.to_array()
# 
#     ## Create max array
#     stacked_arrays = dstack((A, B, C, D, E))
#     max_of_stack = stacked_arrays.max(2)
# 
#     ## Convert back to xarray...
#     max_of_stack = max_of_stack.to_dataset()
# 
#     ## Write output .nc file
#     outfile = [outdir, scenario, "_t", t,  ".nc"]
#     outfile = "".join(outfile)
#     max_of_stack.to_netcdf(outfile)
# 
# 
# 
# ## Specify parameters
# ssp = ['ssp1-rcp26', 'ssp3-rcp70', 'ssp5-rcp85']
# years = ['2015', '2020', '2025', '2030', '2035', '2040', '2045', '2050', '2055', '2060', '2065', '2070']
# inpath = '/tempdata/research-cifs/uom_data/gsdms_data/outputs/hurtt'
# outdir = '/tempdata/research-cifs/uom_data/gsdms_data/outputs/hurtt/'
# 
# ## Get input files
# infiles = [os.path.join(inpath, f) for f in os.listdir(inpath) if
# os.path.isfile(os.path.join(inpath, f))]
# 
# ## Loop over parameters
# for i in ssp:
#   for j in years:
# 
#     pattern1 = i
#     reobj1 = re.compile(pattern1)
#     pattern2 = j
#     reobj2 = re.compile(pattern2)
#     pattern3 = ["c3", "c4"]
#     reobj3 = '|'.join(pattern3)
#     reobj3 = re.compile(reobj3)
# 
#     temp = []
# 
#     for file in infiles:
#       if re.findall(reobj1, file):
#         if re.findall(reobj2, file):
#           if re.findall(reobj3, file):
#             temp.append(file)
#             
#     
#     
#     ## >> Run function
#     max_nc(i, j, temp[0], temp[1], temp[2], temp[3], temp[4], outdir)





## >> Debugging
ssp = ['ssp5-rcp85']
years = ['2070']
inpath = '/tempdata/research-cifs/uom_data/gsdms_data/outputs/hurtt'
outdir = '/tempdata/research-cifs/uom_data/gsdms_data/outputs/hurtt/'

infiles = [os.path.join(inpath, f) for f in os.listdir(inpath) if
os.path.isfile(os.path.join(inpath, f))]

for i in ssp:
  for j in years:

    pattern1 = i
    reobj1 = re.compile(pattern1)
    pattern2 = j
    reobj2 = re.compile(pattern2)
    pattern3 = ["c3", "c4"]
    reobj3 = '|'.join(pattern3)
    reobj3 = re.compile(reobj3)

    temp = []

    for file in infiles:
      if re.findall(reobj1, file):
        if re.findall(reobj2, file):
          if re.findall(reobj3, file):
            temp.append(file)
            
    print("View an input file : \n", temp[0])

    A = xr.open_dataset(temp[0])
    B = xr.open_dataset(temp[1])
    C = xr.open_dataset(temp[2])
    D = xr.open_dataset(temp[3])
    E = xr.open_dataset(temp[4])
    print("\n\n data type : \n", type(A))
    print("\n data variables : \n", A.keys)
    
    A = xr.Dataset.to_array(A)
    B = xr.Dataset.to_array(B)
    C = xr.Dataset.to_array(C)
    D = xr.Dataset.to_array(D)
    E = xr.Dataset.to_array(E)
    print("\n\n data type : \n", type(A))
    print("\n data variabes : \n", A.dtype.names)

    # A = A.to_dataset(dim="c3nfx")
    # print("dataset variables : \n", A.keys)

    stacked_arrays = dstack((A, B, C, D, E))
    max_of_stack = stacked_arrays.max(2)

    print("\n\n max_of_stack data type : \n", type(max_of_stack))
    print("\n max_of_stack data variabes : \n", A.dtype.names)

    # ## Convert numpy.ndarray to xarray
    # ## https://xarray.pydata.org/en/v0.20.1/generated/xarray.DataArray.html
    # da = xr.DataArray(
    #   data=A,
    #   dims=["lon", "lat", "c3ann"],
    #   coords=dict(
    #     lon=(["x", "y"], lon),
    #     lat=(["x", "y"], lat),
    #     time=time,
    #     reference_time=reference_time),
    #     attrs=dict(
    #       description="Ambient temperature.",
    #       units="degC",
    #       )
    #       )

    ## Write output .nc file
    ## https://stackoverflow.com/questions/55956270/convert-a-numpy-dataset-to-netcdf/55973923#55973923
    outfile = [outdir, i, "_t", j,  ".nc"]
    outfile = "".join(outfile)
    # max_of_stack.to_netcdf(outfile)
    
    import numpy as np
    from netCDF4 import Dataset
    # -----------------------
    nyears = 16;
    unout = 'days since 2000-01-01 00:00:00'
    # -----------------------
    ny, nx = (250, 186)
    lon = np.linspace(9,30,nx);
    lat = np.linspace(50,60,ny);
    
    dataout = np.random.random((nyears,ny,nx)); # create some random data
    datesout = [datetime.datetime(2000+iyear,1,1) for iyear in range(nyears)]; # create datevalues
    # =========================
    ncout = Dataset('myfile.nc','w','NETCDF3'); # using netCDF3 for output format 
    ncout.createDimension('lon',nx);
    ncout.createDimension('lat',ny);
    ncout.createDimension('time',nyears);
    lonvar = ncout.createVariable('lon','float32',('lon'));lonvar[:] = lon;
    latvar = ncout.createVariable('lat','float32',('lat'));latvar[:] = lat;
    timevar = ncout.createVariable('time','float64',('time'));timevar.setncattr('units',unout);timevar[:]=date2num(datesout,unout);
    myvar = ncout.createVariable('myvar','float32',('time','lat','lon'));myvar.setncattr('units','mm');myvar[:] = dataout;
    ncout.close();
    
    
    
    
    ## https://ecco-v4-python-tutorial.readthedocs.io/ECCO_v4_Saving_Datasets_and_DataArrays_to_NetCDF.html
    
    ## numpy.ndarray to netcdf ... under dev
    ## http://pyhogs.github.io/intro_netcdf4.html
    # import netCDF4 as nc4
    # f = nc4.Dataset(outfile,'w', format='NETCDF4') #'w' stands for write
    
    
    # ## Write numpy array as geotiff
    # ## https://gis.stackexchange.com/questions/37238/writing-numpy-array-to-raster-file
    # ## https://borealperspectives.org/2014/01/16/data-type-mapping-when-using-pythongdal-to-write-numpy-arrays-to-geotiff/
    # import sys
    # from osgeo import gdal
    # from osgeo.gdalconst import *
    # 
    # # register all of the GDAL drivers
    # gdal.AllRegister()
    #     
    # # open the image
    # inDs = gdal.Open("/home/payalb/gsdms_r_vol/tempdata/research-cifs/uom_data/gsdms_data/outputs/hurtt/temp.tif")
    # if inDs is None:
    #   print ('Could not open image file')
    #   sys.exit(1)
    # 
    # # read in the crop data and get info about it
    # band1 = inDs.GetRasterBand(1)
    # rows = inDs.RasterYSize
    # cols = inDs.RasterXSize
    # 
    # # create the output image
    # driver = inDs.GetDriver()
    # #print driver
    # outDs = driver.Create(outfile, cols, rows, 1, GDT_Int32)
    # if outDs is None:
    #     print ('Could not create output tif')
    #     sys.exit(1)
    # 
    # outBand = outDs.GetRasterBand(1)
    # 
    # 
    # # write the data
    # outBand.WriteArray(max_of_stack, 0, 0)
    # 
    # # flush data to disk, set the NoData value and calculate stats
    # outBand.FlushCache()
    # outBand.SetNoDataValue(-99)
    # 
    # # georeference the image and set the projection
    # outDs.SetGeoTransform(inDs.GetGeoTransform())
    # outDs.SetProjection(inDs.GetProjection())
    # 
    # del outData
    
    
  
    
    ## >> Notes
    ## https://ecco-v4-python-tutorial.readthedocs.io/ECCO_v4_Saving_Datasets_and_DataArrays_to_NetCDF.html
    ## https://xarray.pydata.org/en/stable/user-guide/reshaping.html
    
    ## Install osgeo: https://gis.stackexchange.com/questions/356574/modulenotfounderror-no-module-named-osgeo-in-virtualenvironment
    

