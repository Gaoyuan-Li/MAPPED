# Prokaryotic Gene Count Matrix Pipeline

This Nextflow DSL2 pipeline processes raw FASTQ files to produce a gene-level count matrix for prokaryotic genomes. Docker is enabled by default, so you do NOT need to supply `-with-docker` when launching.

## Requirements
- Nextflow (>=20.10.0)
- Docker daemon running

## Inputs
- `--workdir` : Path to a directory containing:
  - `samplesheet/samplesheet.csv` with columns `sample`, `run_accession`, `fastq_1`, `fastq_2` pointing to FASTQ files.
  - `seqFiles/fastq/` containing the FASTQ files.
  - `seqFiles/ref_genome/` containing the reference genome files.
- `--genome_fa` : Path to the reference FASTA (`.fna`) file.
- `--annotation_gff` : Path to the genome GFF annotation file.
- `--ref_genome` : Path to the reference FASTA (`.fna`) file.
- `--ref_gff` : Path to the genome GFF annotation file.

## Example
```bash
cd 3_generate_count_matrix
nextflow run main.nf \
  --workdir ~/PhD/iMM_2/test_results \
  --ref_genome ~/PhD/iMM_2/test_results/seqFiles/ref_genome/GCF_000046845.1_ASM4684v1_genomic.fna \
  --ref_gff ~/PhD/iMM_2/test_results/seqFiles/ref_genome/genomic.gff
```

All outputs will be written into subdirectories of `--workdir`:
- `fastqc/` : FastQC reports (raw & trimmed)
- `trimmed/` : Trimmed FASTQ files
- `bbsplit/` : BBSplit output (microbial reads)
- `salmon/` : Salmon quantification results
- `multiqc/` : MultiQC report
- `count_matrix.tsv` in the top level of `--workdir`
