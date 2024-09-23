#TODO
- return to default BD input file names schema
- replace local container by public Seurat image
- interactive shiny app showing expression of marker genes over umap and some statistics
- expand to whole transcriptome (WTA) assay : only plot expression for the targets of the guides, not for all genes.
- fake / synthetic test data
- images in the readme
- (snakemake report)
- R : separate effective and non-effective guideRNAs
- R : maybe try manually input variable genes, if possible, for a better umap. Of course if we cheat and keep guide genes in the same matrix, they will be dirivng the clustering, but this carries no information but convenience.

Snakemake is in pyenv environment.

Container and all the data are in the backup drive One Touch, be sure to mount it before running .sh
Also copied container to ~/containers.

About mount paths:
$HOME, $pwd are mounted automatically. See documentation:
https://apptainer.org/docs/user/latest/bind_paths_and_mounts.html#system-defined-bind-paths

Resources:
- https://github.com/hbctraining/scRNA-seq_online

