configfile: "config.yaml"

rule all:
	input:
		"out/transcripts.csv",
		"out/embeddings.csv",
		"out/metadata.csv",
		"out/expressions.bin.gz",
		"out/cluster_dge.json"

rule seurat:
	input:
		config["count_matrix"]
	output:
		transcripts="out/transcripts_raw.csv",
		embeddings="out/embeddings.csv",
		metadata="out/metadata.csv",
		expressions="out/expressions.mtx",
		differential_expressions="out/seurat_cluster_dge.json"
	script:
		"scripts/seurat.R"

rule pack_matrix:
	input:
		"out/expressions.mtx"
	output:
		"out/expressions.bin"
	shell:
		"./scripts/compress_matrix.py {input} {output}"

rule compress_matrix:
	input:
		"out/expressions.bin"
	output:
		"out/expressions.bin.gz"
	shell:
		"gzip -k {input}"

rule parse_transcripts:
	input:
		"out/transcripts_raw.csv"
	output:
		"out/transcripts.csv"
	shell:
		"./scripts/parse_transcripts.py {input} > {output}"

rule merge_dges:
	input:
		"out/seurat_cluster_dge.json"
	output:
		"out/cluster_dge.json"
	shell:
		"./scripts/merge_dges.py {input} > {output}"
