## ref: https://gis.stackexchange.com/questions/298230/change-no-data-value-geotif-file-with-qgis-or-gdal

#!/bin/bash
 
basename=$(echo "$1" | cut -f 1 -d '.')
mask=${basename}_mask.tif
output=${basename}_edt.tif
nodata=$(gdalinfo $1 | grep "NoData" | cut -d "=" -f 2)
gdal_calc.py --NoDataValue=$2 --calc="A!=${nodata}" --outfile="$mask" -A $1
gdal_calc.py --NoDataValue=$2 --calc="A*B" --outfile="$output" -A $1 -B $mask
