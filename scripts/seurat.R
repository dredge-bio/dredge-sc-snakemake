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


print(paste("--- Loading seurat object from", snakemake@input[[1]]))
dscSeurat <- readRDS(file=snakemake@input[[1]])

# Create a metadata column called `__dredge_clusters` from the current ident so
# that we can retrieve it later
dscSeurat[["X__dredge_cluster"]] <- Idents(dscSeurat)


print('--- Outputting embeddings and metadata...')
write.csv(Embeddings(dscSeurat, reduction = "umap"), snakemake@output[["embeddings"]])
write.csv(dscSeurat@meta.data, snakemake@output[["metadata"]])
write.csv(rownames(dscSeurat[[dscSeurat@active.assay]]), snakemake@output[["transcripts"]])


print('--- Outputting normalized expression matrix...')
writeMM(GetAssayData(dscSeurat, slot="data"), snakemake@output[["expressions"]])


print('--- Outputting differential expressions... (This will take a while)')
differential.expressions = list()
for (i in levels(dscSeurat)) {
	print(paste('    ...finding markers in cluster ', i))
	differential.expressions[[i]] = FindMarkers(dscSeurat, ident.1=i, min.pct=.25, verbose=TRUE)
}
write_json(differential.expressions, snakemake@output[["differential_expressions"]])
