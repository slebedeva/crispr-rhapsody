# A minimal workflow to process single-cell CRISPR assay from BD Rhapsody

This workflow will take count matrices coming from BD Rhapsody type of targeted single-cell transcriptome assay. This kind of assay contain around 10,000 cells and 500 panel genes. It will perform filtering, normalization, visualization of gene expression and dimensionality reduction (UMAP).

### Prerequisites

- [snakemake](https://snakemake.readthedocs.io/en/stable/) 
- [apptainer](https://apptainer.org/)  (former singularity)

### Quick start

0. Build an apptainer image (requires root priviledges):
    ```
    make build
    ```

1. Edit `config.yaml` to enter path to your data. Default path will be `data` folder in this directory.

    ```
    # config.yaml
    data_dir: "data"
    ```

2. Make sure your input data files are called like `[SAMPLE]_DBEC_MolsPerCell_correct_gene_names.csv`, where [SAMPLE] is your sample name. Thus, an example input data will look like:
    ```
    ├── data
    │   ├── [SAMPLE1]_DBEC_MolsPerCell_correct_gene_names.csv
    │   ├── [SAMPLE2]_DBEC_MolsPerCell_correct_gene_names.csv

    ```

3. (Optional) check that paths are alright by performing a snakemake dry run:
    ```
    snakemake -n
    ```
    
4. Run the workflow:
    ```
    make pipeline
    ```

### Data input

We start with the count matrix which is produced by the [BD internal mapping pipeline](https://www.bdbiosciences.com/en-us/products/software/rhapsody-sequence-analysis-pipeline). This pipeline takes care of alignment to reference, quality filtering and error correction. It has cell ID in rows and gene ID in columns:

Cell_Index | Gene1 | Gene2 | Gene3
---|---|---|---
7836734 | 0 | 68 | 0 
4277806 | 0 | 25 | 0

In our specific case, there is also a "guide RNA gene" included in the panel, which allows us to identify the guide(s) expressed in each cell.

### Pipeline steps

All the steps are currently performed in R using [Seurat](https://satijalab.org/seurat/) package. 
- First, the count matrix is separated in two: for gene and guide expression.
- We filter both matrices to contain the same cells, as well as require the cell to have at least 10 guide reads and 1000 gene reads to pass the filter.
- We assign cell identity based on the highest expressed guide RNA for each cell.
- We plot expression levels of guide targets to visually assess guide effectiveness.
- We perform differential gene expression analysis: comparing cells that carry guides against your target gene of interest to cells carrying negative control (non-targeting) guides.
- We perform dimensionality reduction and visualize on UMAP if cells with a specific gene targeted form a separate cluster.

### Expected output

Folder `results` will be created, with the following structure:

```
├── results
│   ├── [SAMPLE]_differential_expression.csv
│   ├── [SAMPLE]_genes.h5Seurat
│   ├── [SAMPLE]_guides.h5Seurat
│   ├── plots
│   │   ├── [SAMPLE]_expression_plot.png
│   │   ├── [SAMPLE]_umap.png
```

- `differential_expression.csv` : genes which are differentially expressed in cells carrying guides to a specific gene compared to cell carrying non-targeting guides. Included are genes with p-value cutoff 0.1.

- `genes.h5Seurat` and `guides.h5Seurat` : raw Seurat objects for gene and guide expression, respectively.

- Plots:
  - `expression_plot.png`: a dot plot showing expression of each target gene in cells grouped by target gene.
  - `umap.png` : dimensionality reduction visualization of all cells, colored by their targeted gene.

Additionally, all the log files will be written to a `logs` folder.
  
### Additional commands

Start a jupyter-lab inside the container with R and python kernels available to interactively use Seurat or [scanpy](https://scverse.org/):

```
make run
```
