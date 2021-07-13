# Setup

This pipeline requires [Snakemake](https://snakemake.readthedocs.io/) to run.

# Configuration

All configuration for generating the files required by DrEdGE is placed in a file named `config.yaml`. Create that file in this directory.

There are two methods for loading your dataset which you can configure in `config.yaml`. The first takes precedence over the second.

## 1. Using a pre-configured Seurat object

If you know how to configure your dataset with Seurat and you want to manually tune clustering, UMAP generation, and cluster labeling yourself, you can load in a pre-configured Seurat object via an `.rds` file. To do so, add the following line to your configuration file:

```
seurat_object: /path/to/your/object.rds
```

This Seurat object should already have been normalized, scaled, and had `FindNeighbors`, `FindClusters`, and `RunUMAP` performed on it.


## 2. Using an expression count matrix

If you do not have experience with Seurat, or you want to let this pipeline run the [typical Seurat workflow](https://satijalab.org/seurat/articles/essential_commands.html) for you, then you need to provide an expression count matrix to be read. To do so, add the following line to your configuration file:

```
count_matrix: /path/to/your/counts/file.matrx
```

# Running

Once you have configured your dataset, run `snakemake -j 1` in this directory. Generated files will be placed in the `out` directory.
