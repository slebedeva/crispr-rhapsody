# run a Jupyter notebook inside container
run:
	bash run_scanpy_seurat.sh

# run snakemake pipeline locally
pipeline:
	snakemake --use-apptainer

