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


## >> Specify max function
def max_nc(scenario, t, infile1, infile2, infile3, infile4, infile5, outdir):

    ## Read in .nc files
    A = xr.open_dataset(infile1)
    B = xr.open_dataset(infile2)
    C = xr.open_dataset(infile3)
    D = xr.open_dataset(infile4)
    E = xr.open_dataset(infile5)


    ## Get variable names
    a = list(A.keys())
    b = list(B.keys())
    c = list(C.keys())
    d = list(D.keys())
    e = list(E.keys())
    
    
    ## Rename variables
    A = A.rename({a[0]: 'crop'})
    B = B.rename({b[0]: 'crop'})
    C = C.rename({c[0]: 'crop'})
    D = D.rename({d[0]: 'crop'})
    E = E.rename({e[0]: 'crop'})
    
    
    ## Stack arrays
    stacked_arrays = xr.concat([A, B, C, D, E], 'band')


    ## Get max values
    max_of_stack = stacked_arrays.max('band')


    ## Write output .nc file
    outfile = [outdir, i, "_crop", "_t", j,  ".nc"]
    outfile = "".join(outfile)
    max_of_stack.to_netcdf(outfile)



## Specify parameters
ssp = ['ssp1-rcp26', 'ssp3-rcp70', 'ssp5-rcp85']
years = ['2015', '2020', '2025', '2030', '2035', '2040', '2045', '2050', '2055', '2060', '2065', '2070']
inpath = '/tempdata/research-cifs/6300-payalb/uom_data/gsdms_data/outputs/hurtt_future'
outdir = '/tempdata/research-cifs/6300-payalb/uom_data/gsdms_data/outputs/hurtt_future/'

## Get input files
infiles = [os.path.join(inpath, f) for f in os.listdir(inpath) if
os.path.isfile(os.path.join(inpath, f))]

## Loop over parameters
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


    ## >> Run function
    max_nc(i, j, temp[0], temp[1], temp[2], temp[3], temp[4], outdir)





# ## >> Debugging
# ssp = ['ssp5-rcp85']
# years = ['2070']
# inpath = '/tempdata/research-cifs/6300-payalb/uom_data/gsdms_data/outputs/hurtt_future'
# outdir = '/tempdata/research-cifs/6300-payalb/uom_data/gsdms_data/outputs/hurtt_future/'
# 
# infiles = [os.path.join(inpath, f) for f in os.listdir(inpath) if
# os.path.isfile(os.path.join(inpath, f))]
# 
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
#     ## View input files
#     print("View an input file : \n", temp[0])
#     print("View an input file : \n", temp[1])
#     print("View an input file : \n", temp[2])
#     print("View an input file : \n", temp[3])
#     print("View an input file : \n", temp[4])
# 
#     ## Open nc files
#     A = xr.open_dataset(temp[0])
#     B = xr.open_dataset(temp[1])
#     C = xr.open_dataset(temp[2])
#     D = xr.open_dataset(temp[3])
#     E = xr.open_dataset(temp[4])
#     print("\n\n data type : \n", type(A))
#     print("\n", A)
#     print("\n", B)
#     print("\n", C)
#     print("\n", D)
#     print("\n", E)
#   
# 
#     ## Stack arrays
#     stacked_arrays = xr.merge([A, B, C, D, E])
#     ## same as: stacked_arrays = xr.concat([A, B, C, D, E], 'band', data_vars = "different")
#     print("\n\n", stacked_arrays)
# 
#     
#     ## Get variable names
#     a = list(A.keys())
#     b = list(B.keys())
#     c = list(C.keys())
#     d = list(D.keys())
#     e = list(E.keys())
#     print(type(a))
#     print(a[0])
#     
#     
#     ## Rename variables
#     ## a[0] gets first element of list because input required here is a string
#     A = A.rename({a[0]: 'crop'})
#     B = B.rename({b[0]: 'crop'})
#     C = C.rename({c[0]: 'crop'})
#     D = D.rename({d[0]: 'crop'})
#     E = E.rename({e[0]: 'crop'})
#     print("\n\n", A)
#     print("\n", B)
#     print("\n", C)
#     print("\n", D)
#     print("\n", E)
# 
#     ## Stack arrays
#     stacked_arrays = xr.concat([A, B, C, D, E], 'band')
#     print("\n\n", stacked_arrays)
# 
#     ## Get max values
#     max_of_stack = stacked_arrays.max('band')
#     print("\n\n", max_of_stack)
# 
#     ## Write output .nc file
#     ## https://stackoverflow.com/questions/55956270/convert-a-numpy-dataset-to-netcdf/55973923#55973923
#     outfile = [outdir, i, "_crop", "_t", j,  ".nc"]
#     outfile = "".join(outfile)
#     max_of_stack.to_netcdf(outfile)
    
