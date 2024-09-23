## run this from the directory where container itself is

#!/bin/bash

# local cache for apptainer
mkdir -p ${HOME}/tmp/.apptainer
export APPTAINER_CACHEDIR=${HOME}/tmp/.apptainer

# otherwise permission error see https://stackoverflow.com/a/60821998
export JUPYTER_RUNTIME_DIR=${HOME}/tmp

# mount extra directory in case you need to install python/R packages so they will not screw up you system installaiton 
localdir=".local_seurat"
mkdir -p $localdir

# run container starting notebooks from this directory
singularity exec -B ${localdir}:${HOME}/.local -B /media /home/${USER}/containers/scanpy_seurat.sif  jupyter-lab 

