## ref: https://gis.stackexchange.com/questions/298230/change-no-data-value-geotif-file-with-qgis-or-gdal

#!/bin/bash
 
basename=$(echo "$1" | cut -f 1 -d '.')

## Output filenames
mask=${basename}_temp.tif
output=${basename}_edt.tif

## First it gets the current no data value, then it creates a mask where the data is set as 1 
nodata=$(gdalinfo $1 | grep "NoData" | cut -d "=" -f 2)
gdal_calc.py --NoDataValue=$2 --calc="A!=${nodata}" --outfile="$mask" -A $1

## and then multiplies that mask to get new image where no data is set with a parameter.
gdal_calc.py --NoDataValue=$2 --calc="A*B" --outfile="$output" -A $1 -B $mask
