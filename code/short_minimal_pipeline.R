##### minimal processing of BD Rhapsody targeted panel with CRISPR assay #####


# libraries

library("Seurat")
library("SeuratDisk")
library("ggplot2")


# functions to separate guides and genes
# TODO: write a package if this section expands too much

get_genes=function(obj){rownames(obj) |> (\(.) .[!grepl('guide',.)])()}
get_guides=function(obj){rownames(obj) |> (\(.) .[grepl('guide',.)])()}


# input file as commandline argument
args = commandArgs(trailingOnly = TRUE)
if( length(args)!=1 ){
  stop("Please provide the input file name as the first argument.")
} 

mysample = args[1]
    
message("processing file... ", mysample)


# Hardcode output directory for now
# TODO as config param
outdir = "results"

if(!dir.exists(outdir)){dir.create(outdir)}

# subdir for plots
plotdir = file.path(outdir, "plots")
if(!dir.exists(plotdir)){dir.create(plotdir)}

# 0. clean/harmonize gene names in the BD sevenbridges pipeline output

# Project name (MeCP2 or VPR) and output directory
myproject = sub("_DBEC_MolsPerCell_correct_gene_names.csv", "",
                sub(".+/", "", mysample))
message("Processing project: ", myproject)

# 1. import (already clean) count matrix and make Seurat object
sample = CreateSeuratObject(
    # remember to transpose the matrix
    # why read.csv2: I prefer to use base R if it doesn't make a difference in speed
    counts=read.csv2(
        mysample
        ,row.names=1, comment.char = "#", header=T,sep=",",check.names=F)|>t()
    ,assay='RNA', project=myproject, min.cells=0, min.features=1, check.matrix=T)

# 2. separate into guides and genes matrix

genes = subset(sample, features = get_genes(sample))
guides = subset(sample, features = get_guides(sample))

# 3. filter out cells with less than 10 guide reads and less than 1000 gene reads

message("removing cells with <10 guide reads and <1000 gene reads...")
guides = subset(guides, subset = nCount_RNA > 10)
genes = subset(genes, subset = nCount_RNA > 1000)

# 4. intersect the matrices and assign top expressed guide

## filter the objects again for common cells

## common cells between the two:
mycells=intersect(Cells(genes), Cells(guides))
genes=subset(genes, cells=mycells)
guides=subset(guides, cells=mycells)


## assign top expressed guide
message("Assingning top expressed guide to cell...")
# prepare guide and target gene metadata
idx = apply(guides@assays$RNA@counts, 2, which.max) # idx of column where max value is == guide identity index
guide_identity = rownames(guides)[idx] # retrieve guide identity
target_gene_identity=gsub('-guide[0-9]+','',guide_identity) # retrieve target gene identity

#add guide metadata to filtered genes seurat object
genes = AddMetaData(object = genes, metadata=guide_identity, col.name='guide_identity')
genes = AddMetaData(object = genes, metadata=target_gene_identity, col.name='target_gene')

# assign top expressed guide as cell identity
Idents(genes)=factor(target_gene_identity, levels=sort(unique(target_gene_identity)))

# TODO: statistics about how many cells/guide
#genes@meta.data |> dplyr::group_by(target_gene) |> dplyr::summarize(count = dplyr::n()) 

# 5. Normalize data

genes = genes |> NormalizeData()

# Optional: plotting (after normalizaiton)
target_genes = unique(target_gene_identity) |> (\(.) .[!grepl('non-targeting',.)])() |> sort()
expression_plot = DotPlot(genes, features = target_genes) + RotatedAxis()
ggplot2::ggsave(file.path("results", "plots", paste0(myproject, "_expression_plot.png")), dpi = 150, width = 15, height = 15, bg = "white")

# TODO: optional: regress out cell cycle effects

# 6. Default differential expression - between targeted cells and non-targeting

message("Calling differentially expressed genes...")

# TODO set as config parameter
p_val_cutoff = 0.1

# make a list of DE genes for each guide target 

DE_list = list()

for(target_gene in target_genes){
    DE_list[[target_gene]] = FindMarkers(genes, ident.1 = target_gene, ident.2 = "non-targeting") |>  (\(.) .[.$p_val_adj<p_val_cutoff,])()}

# unite into one data.frame
DE_df = purrr::map_df(DE_list, rbind, .id="guide_target")

# 7. visualize and save umap
message("Performing dimensionality reduction...")
genes = genes |> FindVariableFeatures() |> ScaleData()
genes = RunPCA(object = genes)
genes = RunUMAP(object = genes,  dims = 1:40)
#genes_p = DimPlot(genes, group.by='target_gene', split.by='target_gene', ncol=4)
genes_p = DimPlot(genes, reduction = "umap", pt.size = 1) + theme(legend.position = "bottom")
ggsave(plot = genes_p, filename = file.path("results", "plots", paste0(myproject, "_umap.png")), dpi = 150, width = 10, height = 10, bg = "white")

# save data table 
#saveRDS(DE_df, file = paste0(myproject, "_differential_expression.Rds"))
# remember to convert row names to column
data.table::fwrite(tibble::rownames_to_column(DE_df, 'gene_name'), file = file.path(outdir ,paste0(myproject, "_differential_expression.csv")))

# save Seurat objects
SaveH5Seurat(object=guides, filename = file.path(outdir,paste0(myproject,"_guides.h5Seurat")), overwrite = TRUE)
SaveH5Seurat(object=genes, filename = file.path(outdir,paste0(myproject,"_genes.h5Seurat")), overwrite = TRUE)
