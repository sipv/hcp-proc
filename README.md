

# HCP-proc

Preprocessing of Human Connectome Project data for use in whole-brain modeling.



### Goals

  - Build the structural connectome 
  - Extract the parcellated resting state time series

Parcellation used is Desikan-Killiany (`aparc+aseg` from FreeSurfer), and both the connectome and the resting state time series are built for all cortical and subcortical regions (`res_all/` directory in the generated data) and for cortical regions only (`res_ctx/` directory in the generated data).
The fMRI time series are extracted both without further processing (apart from the ICA-FIX done in HCP pipeline), and also processed by DiCER [1].


### Prerequisities

1. HCP data (available from https://db.humanconnectome.org)

   Data download is the first step of the Snakemake worflow. To successfully execute it, you will need AWS S3 client, and AWS credentials (which you can get through HCP website). 


2. Software
   - Connectome workbench
   - MRtrix3
   - Python 3.7 with standard (neuro-)scientific packages
   - DiCER package (if processing fMRI data; https://github.com/BMHLab/DiCER).


### Usage

1. Get the data and the software.
2. Specify path to DiCER with environment variable `DICERPATH`.
3. From the working directory call `snakemake sc fmri`. Test it with a dry-run first (`-n` switch), and make sure you have sufficient computational resources for the tractography.
   The IDs of the desired subjects should be placed in a text file (one ID per line), passed to Snakemake via env variable `SUBJECTS_FILE`.
   Note that some parts of the workflow require considerable amount of memory, so it is useful to specify also the amount of memory available with Snakemake option `--resources mem_mb=X`.


### References

[1] Aquino et al. (2019). Identifying and removing widespread signal deflections from fMRI data: Rethinking the global signal regression problem. NeuroImage 212: 116614. https://doi.org/10.1016/j.neuroimage.2020.116614. 
