#Start
import os
import numpy as np  
from osgeo import gdal
#from gdalconst import GA_ReadOnly
from osgeo import ogr
from osgeo import osr
import rasterio
import operator
import time
import ntpath
import sys
import subprocess

from rasterstats import zonal_stats

import argparse
parser= argparse.ArgumentParser(description='Program to generate australian future bioclim from global future bioclim')
parser.add_argument('-i','--input_file_path', type=str, metavar='', help='specify the folder where the tif files to be averaged are')
parser.add_argument('-m','--mask_path', type=str, metavar='', help='specify the folder where the tif files to be averaged are')
parser.add_argument('-o','--output_directory', type=str, metavar='', help='specify the folder where the resulting file will be stored')
parser.add_argument('-r','--resolution', type=int, metavar='', help='desired resolution in meters')

args = parser.parse_args()


# Defining input, mask and output directories
input_file_path=args.input_file_path #"/home/ubuntu/mnt/Alex/gsdms_alex/outputs/out_ssp126_2021-2040"
mask_path=args.mask_path
output_dir= args.output_directory #"/home/ubuntu/mnt/Alex/aus-ppms_alex/outputs"


## String for the wgs projection
equalearth_crs = '+proj=eqearth +ellips=WGS84 +wktext'
wgs_crs = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0" 

new_extent = "-180 -60 180 90"


res=args.resolution

new_crs = "EPSG:4326" 

gsdms_dir = "/home/ubuntu/mnt3" # "/tempdata/research-cifs/6300-payalb/uom_data/gsdms_data" in Payal's drive
aus_dir = "/home/ubuntu/mnt2" # "/tempdata/research-cifs/6300-payalb/uom_data/aus-ppms_data" in Payal's drive
data_dir = "/home/ubuntu/mnt/aus-ppms_alex/data"
temp_dir="/home/ubuntu/mnt/Alex/aus-ppms_alex/temp" #Where I'll store temporal files produced by intermediate steps


#
def get_output_name(output_path):
    
    print(f'out is file? {output_path}')
    if str(output_path).endswith('/') or str(output_path).endswith('"\"'):
        print("is NOT file")
        output_folder=os.path.dirname(output_path)
        output_name=None
    else:
        print("is file")
        output_folder=os.path.dirname(output_path)
        output_name=ntpath.basename(output_path)
    print(output_folder)
    print(output_name)
    return output_folder , output_name

def add_suffix_step(infile_path, suffix):
    '''Adds the suffix specified depending on the extension of the input
    '''
    filename=ntpath.basename(infile_path)
    str_infile_path=str(infile_path)
    str_infile_path=str_infile_path.lower()
    if str_infile_path.endswith('.adf'):
        filename=filename.replace('.adf',suffix)
    elif str_infile_path.endswith('.dat'):
        filename=filename.replace('.dat',suffix)
    elif str_infile_path.endswith('.tif'):
        filename=filename.replace('.tif',suffix)

    return filename

# TODO refactor this change_nodata_value to take any data value in the arguments
def change_nodata_value(infile_path, outfolder):
    #filename=ntpath.basename(infile_path)
    #out_filename=filename.replace(".tif","nodata.tif") #TODO Handle this assumption later

    out_filename=add_suffix_step(infile_path,"nodata.tif")
    outfile_path=os.path.join(outfolder,out_filename)

    os.system(f"gdalwarp -dstnodata -9999 {infile_path} {outfile_path} ")

    return outfile_path

def mask_ones(infile_path,outfolder):

    #filename=ntpath.basename(infile_path)
    #out_filename=filename.replace(".tif","mask_ones.tif") #TODO Handle this assumption later
    out_filename=filename.add_suffix_step(infile_path,"mask_ones.tif")
    
    outfile_path=os.path.join(outfolder,out_filename)
    print(f"inside mask_ones")
    call=f'gdal_calc.py -A {infile_path} --outfile {outfile_path} --calc="numpy.where(A!=-9999,1,-9999)" --NoDataValue=-9999'  
    os.system(call)
    print(f"finished call {call}")
    return outfile_path


def get_shapefile(infile_path):
    
    outfile_path = str(infile_path).replace(".tif",".shp")
    
    out_filename = ntpath.basename(outfile_path)

    
    call = f"gdaltindex {outfile_path} {infile_path}"

    os.system(call)

    print(f"\nGetting index for {infile_path}\n")

    return outfile_path


def get_stats_raster(infile_path, silent=True):
    shapefile_path=get_shapefile(infile_path)
    if not silent:
       print(f'Getting stats for {ntpath.basename(infile_path)}')

    stats = zonal_stats(shapefile_path, infile_path, stats=['min', 'max', 'median', 'majority', 'range' , 
                                                            'nodata', 'percentile_25.0', 'percentile_75.0'])
    
    
    # Deletes .shp, .dbf, .prj and .shx files
    target_dir=os.path.dirname(infile_path)
    junk_extensions=[".dbf",".prj",".shx",".shp"]
    for j in junk_extensions:
        command=f'rm *{j}'
        p=subprocess.Popen(command, cwd=target_dir, shell=True)
        p.wait()

    print(stats[0])
    raster_stats=stats[0]

    return raster_stats

def get_stats_directory(input_dir, output_dir=None):
    #Creating summary of stats for all the tif files in the output directory
    
    infiles=[]
    for file in os.listdir(input_dir):
        # check only text files
        if file.endswith('.tif'):
            infiles.append(file)
    
    if not output_dir:
        output_dir=input_dir

    stats_out_path=os.path.join(output_dir,'report.txt')
    with open(stats_out_path,'w') as f:
        for filename in infiles:
            file_to_stats=os.path.join(output_dir,filename)
            f.write(f'\nStats for {filename} :\n')
            stats=get_stats_raster(file_to_stats)
            for key, value in stats.items():
                f.write(f'{key} : ')
                f.write(f'{value}')
                f.write('\n')

def translate(crs, infile_path, outfolder, out_name=None):
    
    #Get the infile name from the path
    #filename=ntpath.basename(infile_path)
    if out_name is None:
        #out_filename=filename.replace(".tif","_wgs.tif")
        out_filename=add_suffix_step(infile_path,"_wgs.tif")
    else:
        out_filename=out_name

    outfile_path=os.path.join(outfolder,out_filename)

    #Makes system call to gdal
    os.system(f"gdal_translate -ot Float32 -a_srs '{crs}' {infile_path} {outfile_path} ")

    return outfile_path

def set_extent(extent,infile_path, outfolder, out_name=None):
    #Get the infile name from the path
    #filename=ntpath.basename(infile_path)
    if out_name is None:
        out_filename=add_suffix_step(infile_path,"_clip.tif")
    else:
        out_filename=out_name

    outfile_path=os.path.join(outfolder,out_filename)
    
    if isinstance(extent, str):
        print(f"gdalwarp -overwrite -ot Float32 -te {extent} {infile_path} {outfile_path}")
        os.system(f"gdalwarp -overwrite -ot Float32 -te {extent} {infile_path} {outfile_path}")
    elif isinstance(extent,list):
        os.system(f"gdalwarp -overwrite -ot Float32 -te {extent[0][0]} {extent[0][1]} \
                                                {extent[1][0]} {extent[1][1]} \
                                                {infile_path} {outfile_path}")
    return outfile_path

def get_extent(infile_path):

    """ Return list of corner coordinates from a gdal Dataset
        Taken from: https://gis.stackexchange.com/questions/57834/how-to-get-raster-corner-coordinates-using-python-gdal-bindings
    """
    ds=gdal.Open(infile_path)

    xmin, xpixel, _, ymax, _, ypixel = ds.GetGeoTransform()
    width, height = ds.RasterXSize, ds.RasterYSize
    xmax = xmin + width * xpixel
    ymin = ymax + height * ypixel

    all_corners=[(xmin, ymax), (xmax, ymax), (xmax, ymin), (xmin, ymin)] #maybe useful later
    min_coord, max_coord =(xmin, ymin), (xmax, ymax)   
    
    return [min_coord , max_coord]

def reproject(crs, resolutionX, resolutionY, infile_path, outfolder, out_name=None):
    #Get the infile name from the path
    #filename=ntpath.basename(infile_path)
    if out_name is None:
        out_filename=add_suffix_step(infile_path,"_reproj.tif")
    else:
        out_filename=out_name

    outfile_path=os.path.join(outfolder,out_filename)

    os.system(f"gdalwarp -overwrite -ot Byte -tr {resolutionX} {resolutionY} -t_srs '{crs}' {infile_path} {outfile_path}")

    return outfile_path

def reproject_with_resampling(resamp_method, resolutionX,resolutionY, s_crs, t_crs, infile_path, outfolder, out_name=None):
    #Get the infile name from the path
    #filename=ntpath.basename(infile_path)
    if out_name is None:
        out_filename=add_suffix_step(infile_path,"_reproj.tif")
    else:
        out_filename=out_name

    
    outfile_path=os.path.join(outfolder,out_filename)

    os.system(f"gdalwarp -overwrite -ot Float32 -r {resamp_method} -tr {resolutionX} {resolutionY}\
                                                 -s_srs '{s_crs}' -t_srs '{t_crs}' {infile_path} {outfile_path}")
    return outfile_path

def apply_mask(infile_path, mask_path, no_data_val, outfolder, out_name=None):
    #Get the infile name from the path
    #filename=ntpath.basename(infile_path)
    if out_name is None:
        out_filename=add_suffix_step(infile_path,"_masked.tif")
    else:
        out_filename=out_name

    outfile_path=os.path.join(outfolder,out_filename)

    #Processing nodataval
    nodatavalue=str(no_data_val)
    
    print(f"no data value is being set to : {nodatavalue}")
    
    os.system(f"gdal_calc.py -A {infile_path} -B {mask_path} --calc='((B==1)*A)+(-9999*(B!=1))' --NoDataValue={nodatavalue} --outfile={outfile_path}")

    return outfile_path


def run_pipeline(infile_path,mask_path,new_crs,new_extent, resolution, s_crs, t_crs, output_path,temp_dir=None):
    '''

    output_path: Path to output. Can be just a directory or a file 
    '''

    output_dir , output_name = get_output_name(output_path)
    
    if not os.path.exists(output_dir):
        print(f'creating output folder: {output_dir}')
        os.makedirs(output_dir)

    if temp_dir is None:
        temp_dir=os.path.join(output_dir,"temp")

    if not os.path.exists(temp_dir):
        os.makedirs(temp_dir)

    # >> Step 1: Check CRS. Set it to EPSG:4326 if it is different
    temp_in=infile_path

    with rasterio.open(temp_in) as src:
        ori_crs=str(src.crs)
        print (f"\nOriginal CRS is:{str(src.crs)}\n")
        if ori_crs != new_crs:
            print(f"Needs to be translated to {new_crs}")
            temp_out = translate(new_crs, temp_in, temp_dir)
            temp_in=temp_out



    ## >> Step two_Reduce layer extent, as specified in WGS 84 ####
    
    print("\n Setting extent \n")
    temp_out=set_extent(new_extent,temp_in,temp_dir)# TODO change this to temp_out after coding step 1


    ## >> Step three_Reproject layer to Equal earth ####
    temp_in=temp_out
    print("\n Reprojecting to equal earth \n")
    temp_out=reproject_with_resampling("bilinear", resolution, resolution, s_crs, t_crs, temp_in, temp_dir, out_name=None)

    ## >> Step four_Clip layer s.t. #cells in layer equal #cells in mask ####
    print("\n Clipping to mask \n")
    temp_in=temp_out
    extent_mask = get_extent(mask_path)
    temp_out=set_extent(extent_mask,temp_in,temp_dir)

    
    ## >> Step five_Change nodata values to -9999 ####
    ## To remove nodata value of -3.4e+38
    ## This is not needed for bio_future layers (check)
    temp_in=temp_out
    print("\n Changing NoDataValue \n")
    temp_out=change_nodata_value(temp_in, temp_dir)

    print("\n Applying mask \n")
    ## >> Step six_Mask ####
    ## https://gitlab.unimelb.edu.au/garberj/gdalutilsaddons/-/blob/master/gdal_calc.R    
    temp_in=temp_out
    final_out=apply_mask(temp_in,mask_path,-9999,output_dir, output_name)
    
    #print(get_stats_raster(final_out))


def run_pipeline_directory(input_dir,mask_path,new_crs,new_extent,resolution,s_crs, t_crs, output_dir,temp_dir=None):
    rasters = []

    # Iterate directory 
    for file in os.listdir(input_dir):
        # check only .tif files
        if file.endswith('.tif'):
            rasters.append(file)
    
    for raster in rasters:
        temp_in=os.path.join(input_dir,raster)
        run_pipeline(temp_in,mask_path,new_crs,new_extent, resolution, s_crs, t_crs, output_dir,temp_dir=None)

    #Now I need the stats of the files I just processed.
    #Therefore, the output_dir becomes the new input_dir for the stats function
    input_dir=output_dir 

    get_stats_directory(input_dir)

def create_mask(infile_path,output_dir,crs,extent,resolution):
    
    temp_dir=os.path.join(output_dir,"temp")

    if not os.path.exists(temp_dir):
        os.makedirs(temp_dir)
    
    temp_in=infile_path
    temp_out= change_nodata_value(temp_in,temp_dir)

    temp_in = temp_out
    temp_out=mask_ones(temp_in,temp_dir)

    temp_in=temp_out

    # Setting crs for the mask
    print(f"\ntranslating {temp_in}\n")
    temp_out=translate(crs, temp_in, temp_dir, "globalmask_wgs_30s.tif" )

    with rasterio.open(temp_out) as src:
        print ("CRS according to rasterio is:",type(src.crs),str(src.crs))

    # Checking that we have the values we want in the mask (1's  for value and -9999 for no data value)
    #print(get_stats_raster(temp_in))

    # Setting extent to the mask
    temp_in=temp_out
    print(f"\nSetting extent {temp_in}\n")
    temp_out= set_extent(new_extent,temp_in,temp_dir, "globalmask_wgs_30s_clip.tif" )

    #Reproject the mask to resolution of 10km
    temp_in=temp_out
    print(f"\nReprojecting {temp_in}\n")
    temp_out=reproject(equalearth_crs, resolution, resolution,temp_in, output_dir, f"globalmask_ee_{resolution/1000}km_nodata.tif")

    #Note: For step 4: NoData value in mask was decided earlier to be -9999

    return temp_out
    temp_in=infile_path


if __name__=="__main__":

    # Handled in the arguments
    #infile_path = "/home/ubuntu/mnt3/soil/71a04c738fde75ee64f57118ed466530/IGBPDIS_SURFPRODS/bulkdens.dat"

    #"/home/ubuntu/mnt3/bio_30s/bio_current_1.tif"
    
    #input_dir="/home/ubuntu/mnt3/bio_30s/"
    #mask_path="/home/ubuntu/mnt/Alex/gsdms_alex/outputs/t2_out/10km/globalmask_ee_10.0km_nodata.tif"
    #output_dir="/home/ubuntu/mnt/Alex/gsdms_alex/outputs/t2_out/10km/bio_current/"
    #if not os.path.exists(output_dir):
    #    os.makedirs(output_dir)

    #Mask creation
    #mask_path=create_mask(infile_path,output_dir, new_crs,new_extent, res)
    #print(get_stats_raster(mask_path))

    #print(get_stats_raster("/home/ubuntu/mnt/Alex/gsdms_alex/temp/mask/globalmask_ee_10km_nodata.tif"))

    
    #######Layers processing
    #Runs all inside a directory
    run_pipeline_directory(input_file_path,mask_path, new_crs, new_extent, res, wgs_crs, equalearth_crs, output_dir, temp_dir=None)
    
    #Runs only one file at a time
    #run_pipeline(input_file_path,mask_path,new_crs, new_extent, res, wgs_crs, equalearth_crs, output_dir, temp_dir=None) #TODO put this in an if
    
    
    
    #print(get_stats_raster("/home/ubuntu/mnt3/outputs/layers_10k/soil/bulkdens_wgs_ee.tif"))
    
    #print(get_stats_raster("/home/ubuntu/mnt/Alex/gsdms_alex/temp/proclay/bulkdens_wgs_clip_reproj_clip2nodata_masked.tif"))
    
    #print(get_stats_raster("/home/ubuntu/mnt3/outputs/byte_layers/layers_10k_byte/bulkdens_wgs_ee.tif"))

    #print(get_stats_raster("/home/ubuntu/mnt3/outputs/layers_10k/soil/bulkdens_wgs_ee.tif"))

    

    




    


    
    


    



