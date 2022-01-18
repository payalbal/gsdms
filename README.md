# Species distribution models using bespoke ppm code
scripts for running point process models for many species

# Biodiversity data
Citation: 

# Covariate data sources


# Workflow
1. Download & clean GBIF data: https://github.com/payalbal/gbifprocessing

2. Download covariate data: data_downloads

3. Covariate data processing: data_processing

4. Model fitting: ppm_fitting

5. Model predictions: ppm_predictions

6. Mapping: 

Functions:

Ancillary script files





Notes:
4.Move species occurrence points falling off the mask to nearest 'land' cells
We lose data again i.e. number of unique locations is reduced. This can be problematic for ppms...

nrow(unique(outside_pts))
nrow(unique(land))
  
  
5. Extract covariates for presence points
For landuse: Take the raster value with lowest distance to point AND non-NA value in the raster
This takes very logn to run for global layers; find alternative...


6. Catch errors in models
Think of how to reduce these. Possible quasiseperation issues due to spatially restricted...could resolved when biuilding models on reduced area?
At the moment, the script finds species with errors and reruns the models. If errors persist, species with errors are discarded.


7. LOCATE ERRORS AND RERUN ANALYSIS FOR SPECIES WITH ERROR
To be automated if error problem is not solved by species grouping
At the moment, it appears that error might be when species data is spatially restricted.
