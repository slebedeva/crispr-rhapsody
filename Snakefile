configfile: "config.yaml"

DATADIR = config["data_dir"]
#SEURAT_IMAGE = config["seurat_image"]

# discover all samples in the input directory
SAMPLE, = glob_wildcards(DATADIR + "/{sample}_DBEC_MolsPerCell_correct_gene_names.csv")

rule all:
  input:
    expand("results/{sample}_differential_expression.csv", sample=SAMPLE),

rule run_seurat:
    input: 
        DATADIR + "/{sample}_DBEC_MolsPerCell_correct_gene_names.csv"
    output:
        "results/{sample}_differential_expression.csv"
    log:
        "logs/run_seurat.{sample}.log",
        stderr = "logs/run_seurat.{sample}_err.log",
        stdout = "logs/run_seurat.{sample}_out.log"
    container:
        "docker://satijalab/seurat:4.3.0"
    shell:
        "Rscript code/short_minimal_pipeline.R {input}"



