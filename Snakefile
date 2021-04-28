rule all:
	input:
		"out/transcripts.csv",
		"out/embeddings.csv",
		"out/metadata.csv",
		"out/expressions.bin.gz"

rule seurat:
	input:
		"../from-julian-2021-01-28/sc.counts_J.matrix"
	output:
		transcripts="out/transcripts_raw.csv",
		embeddings="out/embeddings.csv",
		metadata="out/metadata.csv",
		expressions="out/expressions.mtx"
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
