# Step 4: Generate Gene Count Matrix Pipeline

This Nextflow DSL2 pipeline processes raw FASTQ files to produce a gene-level count matrix for prokaryotic genomes.

## Prerequisites

- Nextflow
- Docker

## Usage

```bash
nextflow run main.nf \
  --workdir ~/PhD/iMM_2/test_results \
  --ref_genome ~/PhD/iMM_2/test_results/seqFiles/ref_genome/GCF_000046845.1_ASM4684v1_genomic.fna \
  --ref_gff ~/PhD/iMM_2/test_results/seqFiles/ref_genome/genomic.gff
```

## Parameters

- `--workdir <Path>`: Directory containing:
  - `samplesheet/samplesheet.csv` with columns `sample`, `run_accession`, `fastq_1`, `fastq_2`.
  - `seqFiles/fastq/` containing the FASTQ files.
  - `seqFiles/ref_genome/` containing the reference genome files.
- `--ref_genome <Path>`: Path to the reference FASTA (`.fna`) file.
- `--ref_gff <Path>`: Path to the genome GFF annotation file.

## Outputs

- `fastqc/`: FastQC reports.
- `trimmed/`: Trimmed FASTQ files.
- `salmon/`: Salmon quantification results.
- `multiqc/`: MultiQC report.
- `expression_matrix/`: TPM and counts matrices