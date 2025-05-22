# Step 4: Generate Gene Count Matrix Pipeline

This Nextflow DSL2 pipeline processes raw FASTQ files to produce a gene-level count matrix for prokaryotic genomes.

## Prerequisites

- Nextflow
- Docker

## Usage

```bash
nextflow run main.nf --outdir ../test_results/  -resume
```

## Parameters

- `--outdir <Path>`: Directory containing:
  - `samplesheet/samplesheet.csv` with columns `sample`, `run_accession`, `fastq_1`, `fastq_2`.
  - `seqFiles/fastq/` containing the FASTQ files.
  - `seqFiles/ref_genome/` containing exactly one FASTA (`.fna`/`.fa`) and one GFF (`.gff`) file, which will be auto-detected.

## Outputs

- `fastqc/`: FastQC reports.
- `trimmed/`: Trimmed FASTQ files.
- `salmon/`: Salmon quantification results.
- `multiqc/`: MultiQC report.
- `expression_matrix/`: TPM and counts matrices