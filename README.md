# MAPPED (Modular Analysis Pipeline for Prokaryotic Expression Data)

## Overview
MAPPED is a comprehensive Nextflow-based workflow for analyzing prokaryotic RNA-seq data. The pipeline automates the entire process from data acquisition to expression quantification through four main steps:

1. Download metadata from SRA
2. Download FASTQ files
3. Download reference genome
4. Generate count matrices

## Requirements
- Nextflow (>=22.10.0)
- Docker

## Usage

```bash
# Run the complete workflow
nextflow run main.nf --organism 'Acinetobacter baylyi' --outdir results --library_layout paired

# Resume a previous run
nextflow run main.nf --organism 'Acinetobacter baylyi' --outdir results --library_layout paired -resume
```

## Parameters

| Parameter | Description | Default |
|-----------|-------------|---------|
| `--organism` | Organism name (e.g. 'Acinetobacter baylyi') | (required) |
| `--outdir` | Output directory | (required) |
| `--library_layout` | Library layout (paired/single/both) | paired |
| `--skip_fastq_download` | Skip FASTQ download step | false |
| `--download_method` | Download method (ftp/aspera/sratools) | ftp |

## Output Structure

```
outdir/
├── metadata/              # Metadata files
├── samplesheet/           # Sample sheets for downstream analysis
├── seqFiles/              # Reference genome files
│   └── ref_genome/
├── fastqc/                # FastQC reports
├── trimmed/               # Trimmed FASTQ files
├── salmon/                # Salmon quantification results
├── expression_matrices/   # TPM and count matrices
└── multiqc/               # MultiQC report
```

## Individual Steps
Each step can also be run independently:

```bash
# Step 1: Download metadata
cd 1_download_metadata_efetch
nextflow run main.nf --organism 'Acinetobacter baylyi' --outdir ../results --library_layout paired

# Step 2: Download FASTQ files
cd 2_download_fastq
nextflow run main.nf --outdir ../results

# Step 3: Download reference genome
cd 3_download_reference_genome
nextflow run main.nf --organism 'Acinetobacter baylyi' --outdir ../results

# Step 4: Generate count matrix
cd 4_generate_count_matrix
nextflow run main.nf --outdir ../results
```