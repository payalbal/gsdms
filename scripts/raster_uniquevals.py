from osgeo import gdal
import numpy as np

def raster_uniquevals(file):
  ds = gdal.Open(file)
  values = np.unique(np.array(ds.GetRasterBand(1).ReadAsArray()))
  return values
