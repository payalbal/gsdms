
import os
import sys
import ntpath
import subprocess
from rasterstats import zonal_stats
import rasterio
import argparse


SUPPORTED_EXTENSIONS=['.tif','.nc','.adf','.dat']

def get_shapefile(infile_path):
    
    supported,ext = is_supported_extension(infile_path)
    if supported:
        outfile_path = str(infile_path).replace(ext,".shp")
        out_filename = ntpath.basename(outfile_path)

        print(f"\nGetting index for {infile_path}\n")
        call = f"gdaltindex {outfile_path} {infile_path}"
        os.system(call)

        return outfile_path


def is_supported_extension(filename):

    for ext in SUPPORTED_EXTENSIONS:
        if filename.endswith(ext):
            #Returns true and extension if the file ends with one of the supported file extensions
            return (True,ext)
        
    return (False,None)


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

    # Now including the nodata value
    with rasterio.open(infile_path) as src:
        stats[0]["nodatavalue"]=src.nodatavals

    print(stats[0])
    raster_stats=stats[0]

    return raster_stats


def handle_directories(input_path, output_dir=None):
    ''' Handles inputs and output paths depending on whether the input path leads to a file or is a directory
    Returns a list of names of input, the path to where those input files are, and also an output folder to store the results
    '''
    #Creating summary of stats for all the tif files in the output directory
    infiles=[] #Contains names of input files (not whole paths)
    if os.path.isdir(input_path): #input path is directory
        for filename in os.listdir(input_path):
            # check only supported files
            if is_supported_extension(filename)[0]:
                infiles.append(filename)
        
        input_dir=input_path
        if not output_dir:
            output_dir=input_dir #outputdir is the same directory as the inpupt file is stored

    else: #input path leads to a file
        infiles=[ntpath.basename(input_path)] #only the name of the file (basename)
        input_dir=os.path.dirname(input_path) # the directory where the file is stored
        
        if not output_dir:
            output_dir=input_dir #outputdir is the same directory as the inpupt file is stored

    if output_dir: #Check if output dir exists only if it is specified by the user
        if not os.path.exists(output_dir):
            print(f'creating folder for stats: {output_dir}')
            os.makedirs(output_dir)

    return infiles , input_dir, output_dir



def get_stats_directory(input_path, output_dir=None, report_name=None):
    
    infiles ,input_dir, output_dir = handle_directories(input_path, output_dir)
    
    if report_name:
        stats_out_path=os.path.join(output_dir,report_name)
    else:
        stats_out_path=os.path.join(output_dir,'report.txt')

    ###This section is where the stats are actually generated and written to a .txt file
    with open(stats_out_path,'w') as f:
        for filename in infiles:
            file_to_stats=os.path.join(input_dir,filename)
            f.write(f'\nStats for {filename} :\n')
            try:
                stats=get_stats_raster(file_to_stats)
                for key, value in stats.items():
                    f.write(f'{key} : ')
                    f.write(f'{value}')
                    f.write('\n')
            except Exception as e:
                print(f"Exception when doing stats for:{filename}")
                print(e)
            


if __name__=="__main__":

    parser= argparse.ArgumentParser(description='Program to generate australian future bioclim from global future bioclim')
    parser.add_argument('-i','--input_file_path', type=str, metavar='', help='specify the folder where the tif files to be averaged are')
    parser.add_argument('-o','--output_directory', type=str, metavar='', required=False ,help='specify the folder where the resulting file will be stored')
    parser.add_argument('-r','--report_name', type=str, metavar='', required=False ,help='specify name of the folder')
    args = parser.parse_args()

    get_stats_directory(args.input_file_path, output_dir=args.output_directory, report_name=args.report_name)

    

