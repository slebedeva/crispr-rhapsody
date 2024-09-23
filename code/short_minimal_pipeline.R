##### minimal processing of BD Rhapsody targeted panel with CRISPR assay #####


# libraries

library("Seurat")
library("SeuratDisk")


# functions to separate guides and genes
# TODO: write a package if this section expands too much

get_genes=function(obj){rownames(obj) |> (\(.) .[!grepl('guide',.)])()}
get_guides=function(obj){rownames(obj) |> (\(.) .[grepl('guide',.)])()}




# input file as commandline argument
args = commandArgs(trailingOnly = TRUE)
if( length(args)!=1 ){
  #stop("Please provide the input file name as the first argument.")
    message("Using default file")
    mysample = "../results/pipeline_sevenbridges/20230303_united_guide_mRNA/MeCP2/MeCP2_DBEC_MolsPerCell_correct_gene_names.csv"
} else {
    mysample = args[1]
}

message("processing file... ", mysample)

# load functions 
source("short_pipeline_helper_functions.R")

# 0. clean/harmonize gene names in the BD sevenbridges pipeline output

# project name (MeCP2 or VPR)
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

guides = subset(guides, subset = nCount_RNA > 10)
genes = subset(genes, subset = nCount_RNA > 1000)

# 4. intersect the matrices and assign top expressed guide

## filter the objects again for common cells

## common cells between the two:
mycells=intersect(Cells(genes), Cells(guides))
genes=subset(genes, cells=mycells)
guides=subset(guides, cells=mycells)


## assign top expressed guide
# prepare guide and target gene metadata
idx = apply(guides@assays$RNA@counts, 2, which.max) # idx of column where max value is == guide identity index
guide_identity = rownames(guides)[idx] # retrieve guide identity
target_gene_identity=gsub('-guide[0-9]+','',guide_identity) # retrieve target gene identity

#add guide metadata to filtered genes seurat object
genes = AddMetaData(object = genes, metadata=guide_identity, col.name='guide_identity')
genes = AddMetaData(object = genes, metadata=target_gene_identity, col.name='target_gene')

# assign top expressed guide as cell identity
Idents(genes)=factor(target_gene_identity, levels=sort(unique(target_gene_identity)))


# 5. Normalize data

genes = genes |> NormalizeData()

# Optional: plotting (after normalizaiton)
target_genes = unique(target_gene_identity) |> (\(.) .[!grepl('non-targeting',.)])() |> sort()
DotPlot(genes, features = target_genes) + RotatedAxis()

# TODO: optional: regress out cell cycle effects

# 6. Default differential expression - between targeted cells and non-targeting

p_val_cutoff = 0.1

# make a list of DE genes for each guide target 

DE_list = list()

for(target_gene in target_genes){
    DE_list[[target_gene]] = FindMarkers(genes, ident.1 = target_gene, ident.2 = "non-targeting") |>  (\(.) .[.$p_val_adj<p_val_cutoff,])()}

# unite into one data.frame
DE_df = purrr::map_df(DE_list, rbind, .id="guide_target")

# 7. visualize a umap
#guides = guides |> FindVariableFeatures() |> ScaleData()
#guides = RunPCA(object = guides)
#guides = RunUMAP(object = guides,  dims = 1:40)
#guides_p = DimPlot(guides, group.by='target_gene', split.by='target_gene', ncol=4, reduction = "umap")

# save data (plots, table) #TODO proper path for results in config
saveRDS(DE_df, file = "../results_mini_pipeline/differential_expression.Rds")
saveRDS(guides_p, file = "../results_mini_pipeline/umap.Rds")

# save Seurat objects
SaveH5Seurat(object=guides, filename = "../results_mini_pipeline/guides.h5Seurat")
SaveH5Seurat(object=genes, filename = "../results_mini_pipeline/genes.h5Seurat")
