## Python extras

## Basename
import os.path
os.path.basename(infile))

py_run_string("import os.path")
py_run_string( "import re" )
py_run_string("out = [outdir + re.sub('_2015-2100', '', os.path.basename('/tempdata/research-cifs/uom_data/gsdms_data/landuse/Hurtt/ssp1-rcp26_2015-2100.nc')) + '_lu_var, '_ttime_slice + '.nc']")

py_run_string( "import numpy as np" )
py_run_string("out = np.array(np.arange(2015, 2101, 1).tolist())")

py_run_string("for item in out: print(item)")


## Calling python function using reticulate
use_python(Sys.which("python"))
## create virtual env
virtualenv_create("hurtt_processing")
use_virtualenv("hurtt_processing")

## install packages
py_install(c("xarray", "scipy"), envname = "hurtt_processing")

## Import packages
xr <- import("xarray")
scipy <- import("scipy")

## source function
source_python("/tempdata/workdir/gsdms/scripts/xarray_operations.py")

## run function
infile = '/tempdata/research-cifs/uom_data/gsdms_data/landuse/Hurtt/ssp1-rcp26_2015-2100.nc'
outfile = '/tempdata/research-cifs/uom_data/gsdms_data/outputs/hurtt/ssp1_rcp26_primf_t2015.nc'
lu_var = 'primf'
tsteps = 0.1

extract_nc(infile = infile, outfile = outfile, lu_var = lu_var, time_slice = tsteps)

## OR...
## Uncomment L1 and l23..
system(paste("python /tempdata/workdir/gsdms/scripts/xarray_operations.py", infile, outfile, lu_var, tsteps))
## same error...