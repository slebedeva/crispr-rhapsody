# build a container which contains both tools for single cell analysis: scanpy and seurat

#!/bin/bash

# local cache for apptainer
mkdir -p ${HOME}/tmp/.apptainer
export APPTAINER_CACHEDIR=${HOME}/tmp/.apptainer

# otherwise permission error see https://stackoverflow.com/a/60821998
export JUPYTER_RUNTIME_DIR=${HOME}/tmp

# mount extra directory in case you need to install python/R packages so they will not screw up you system installaiton 
localdir=".local_seurat"
mkdir -p $localdir

# build container preserving environmental variables (-E)
sudo -E apptainer build scanpy_seurat.sif scanpy_seurat.def
