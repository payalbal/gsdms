## Activate python environment
system(". /home/ubuntu/environments/gsdms_env/bin/activate")
system("bash -c 'source /home/ubuntu/environments/gsdms_env/bin/activate'")
# system("source /home/ubuntu/environments/gsdms_env/bin/activate") 
## source does not work. see https://stackoverflow.com/questions/13702425/source-command-not-found-in-sh-shell


## Create subset tables
system("sudo bash -c /home/ubuntu/pyscripts/gbif/make_subset_gbif.sh")

file.exists()
library(reticulate)
use_condaenv(condaenv = '/home/ubuntu/environments/pyscripts', required = TRUE)
py_run_string('import umap')