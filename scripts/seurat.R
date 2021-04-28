library(Seurat)
library(Matrix)


print(paste("--- Loading count data from", snakemake@input[[1]]))
dscCountsMatrix <- read.table(snakemake@input[[1]], sep = "\t")

write.csv(rownames(dscCountsMatrix), snakemake@output[["transcripts"]])

dscSeurat <- CreateSeuratObject(counts = dscCountsMatrix, min.cells = 5, min.features = 200)

# dscSeurat <- subset(dscSeurat, subset = nFeature_RNA > 500)

print('--- Normalizing and scaling data...')
dscSeurat <- NormalizeData(object = dscSeurat, normalization.method = "LogNormalize", scale.factor = 10000)
dscSeurat <- FindVariableFeatures(object = dscSeurat)
dscSeurat <- ScaleData(object = dscSeurat)
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
