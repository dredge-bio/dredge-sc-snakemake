library(Seurat)
library(Matrix)
library(jsonlite)
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

print('--- Outputting differential expressions... (This will take a while)')
differential.expressions = list()
for (i in levels(dscSeurat)) {
	print(paste('    ...finding markers in cluster ', i))
	differential.expressions[[i]] = FindMarkers(dscSeurat, ident.1=i, min.pct=.25, verbose=TRUE)
}
write_json(differential.expressions, snakemake@output[["differential_expressions"]])
