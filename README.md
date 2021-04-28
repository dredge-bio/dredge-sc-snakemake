To create all files required for DrEdGE configuration, first install [Snakemake](https://snakemake.readthedocs.io/). Then create a file in this directory called `config.yaml`. That file will inform this pipeline where it should look for the expression count matrix for you dataset. It must look like this:

```
count_matrix: /path/to/your/counts/file.matrx
```

Next, run `snakemake -j 1` in this directory. Generated files will be placed in the `out` directory.
