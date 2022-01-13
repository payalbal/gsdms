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

import sys
from osgeo import gdal
from osgeo.gdalconst import *

# register all of the GDAL drivers
gdal.AllRegister()

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

    ## Write output .nc file
    outfile = [outdir, i, "_t", j,  ".nc"]
    outfile = "".join(outfile)
    # max_of_stack.to_netcdf(outfile)
    
    
    # open the image
    inDs = gdal.Open("/home/payalb/gsdms_r_vol/tempdata/research-cifs/uom_data/gsdms_data/outputs/hurtt/temp.tif")
    if inDs is None:
      print ('Could not open image file')
      sys.exit(1)
    
    # read in the crop data and get info about it
    band1 = inDs.GetRasterBand(1)
    rows = inDs.RasterYSize
    cols = inDs.RasterXSize

    # create the output image
    driver = inDs.GetDriver()
    #print driver
    outDs = driver.Create(outfile, cols, rows, 1, GDT_Int32)
    if outDs is None:
        print ('Could not create output tif')
        sys.exit(1)
    
    outBand = outDs.GetRasterBand(1)
    
    
    # write the data
    outBand.WriteArray(max_of_stack, 0, 0)
    
    # flush data to disk, set the NoData value and calculate stats
    outBand.FlushCache()
    outBand.SetNoDataValue(-99)
    
    # georeference the image and set the projection
    outDs.SetGeoTransform(inDs.GetGeoTransform())
    outDs.SetProjection(inDs.GetProjection())
    
    del outData
    
    
  
    
    ## >> Notes
    ## https://ecco-v4-python-tutorial.readthedocs.io/ECCO_v4_Saving_Datasets_and_DataArrays_to_NetCDF.html
    ## https://xarray.pydata.org/en/stable/user-guide/reshaping.html
    

