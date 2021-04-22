########The analysis steps have already been done for the R object. However, if you would like to run through them you are more than welcome to!


# Load libraries
library(Seurat)
library(Matrix)
# library(magrittr)
# library(dplyr)
# require(scales)

# Load the dataset
print('--- Loading count data...')
dscCountsMatrix <- read.table("~/Code/srivastava/from-julian-2021-01-28/sc.counts_J.matrix", sep = "\t")

write.csv(rownames(dscCountsMatrix), snakemake@output[["transcripts"]])

dscSeurat <- CreateSeuratObject(counts = dscCountsMatrix, min.cells = 5, min.features = 200)

# dscSeurat <- subset(dscSeurat, subset.names = "nGene", low.thresholds = 500)
dscSeurat <- subset(dscSeurat, subset = nFeature_RNA > 500)

print('--- Normalizing and scaling data...')
dscSeurat <- NormalizeData(object = dscSeurat, normalization.method = "LogNormalize", scale.factor = 10000)
dscSeurat <- FindVariableFeatures(object = dscSeurat, do.plot = F)
dscSeurat <- ScaleData(object = dscSeurat, display.progress = F)
dscSeurat <- RunPCA(object = dscSeurat, pc.genes = dscSeurat@var.genes, do.print = TRUE)

print('--- Running UMAP clustering...')
dscSeurat <- FindNeighbors(object = dscSeurat, dims = 1:20)
dscSeurat <- FindClusters(dscSeurat, resolution = 0.5, print.output = 0, save.SNN = T)
dscSeurat <- RunUMAP(dscSeurat, reduction = "pca", dims= 1:20, min.dist = 0.6, n.neighbors = 50)

print('--- Outputting embeddings and metadata...')
write.csv(Embeddings(dscSeurat, reduction = "umap"), snakemake@output[["embeddings"]])
write.csv(dscSeurat@meta.data, snakemake@output[["metadata"]])

print('--- Outputting normalized expression matrix...')
writeMM(GetAssayData(dscSeurat, slot="data"), snakemake@output[["expressions"]])


##########The code below is useful for exploring the dataset!

#identifying cluster markers (you will get an output file that contains the genes that are up and down regulated)
# cluster1.markers <- FindMarkers(object = juv.data, ident.1 = "1" ,min.pct = 0.25)
