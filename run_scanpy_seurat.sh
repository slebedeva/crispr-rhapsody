## run this from the directory where container itself is
## bash run_container.sh mycontainer.sif

#!/bin/bash


mkdir -p ${HOME}/tmp/.apptainer
export APPTAINER_CACHEDIR=${HOME}/tmp/.apptainer

# otherwise permission erro see https://stackoverflow.com/a/60821998
export JUPYTER_RUNTIME_DIR=${HOME}/tmp

# path is broken because of the whitespace
#container=/media/${USER}/One\ Touch/work_backup_2023/containers/scanpy_seurat.sif

## mount extra directory in case you need to install python/R packages so they will not screw up you system installaiton 
localdir=".local_seurat"
mkdir -p $localdir

singularity exec -B ${localdir}:${HOME}/.local -B /media /media/${USER}/One\ Touch/work_backup_2023/containers/scanpy_seurat.sif  jupyter-lab --notebook-dir /media/${USER}/One\ Touch/work_backup_2023 
