library(Seurat)
library(future)
library(parallel)

options(future.seed=TRUE)
num.cpus <- detectCores() - 1

if (num.cpus == 0) {
	num.cpus <- 1
}

print(paste("--- Using ", num.cpus, " cores"))

plan('multiprocess', workers=num.cpus)


print(paste("--- Loading count data from", snakemake@input[[1]]))
dscCountsMatrix <- read.table(snakemake@input[[1]], sep = "\t")

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

print('--- Saving RDS object...')
saveRDS(dscSeurat, file=snakemake@output[[1]])
