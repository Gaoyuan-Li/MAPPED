# MAPPED - Modular Automated Pipeline for Public Expression Data

MAPPED (Modular Automated Pipeline for Public Expression Data) is a comprehensive Nextflow-based workflow designed to analyze public RNA-seq data from NCBI SRA. It automates the entire process from metadata retrieval to expression matrix generation, making large-scale transcriptomics analysis accessible and reproducible.

## Overview

MAPPED consists of four integrated modules that work together to process public expression data:

1. **Metadata Download**: Retrieves and formats metadata from NCBI SRA based on organism name
2. **Reference Genome Download**: Obtains reference genome sequences and annotations
3. **FASTQ Download**: Efficiently downloads sequencing data using optimized protocols
4. **Expression Quantification**: Performs quality control, trimming, and gene expression quantification

The pipeline is designed to handle large-scale datasets with built-in error handling, resume capabilities, and resource optimization.

## Batch Processing Version

This directory contains a **batch processing version** of MAPPED that processes RNA-seq samples in batches of 500 to manage disk space efficiently. This is particularly useful for organisms with thousands of samples where processing all at once would require excessive disk space.

### Key Features of Batch Processing Version

- **Batch Processing**: Automatically splits large sample sets into batches of 500 samples
- **Disk Space Management**: Cleans up raw and trimmed sequencing files after each batch
- **Smart Batching**: For remaining samples (N < 500):
  - If N < 250: merges with the previous batch
  - If N ≥ 250: creates a separate batch
- **Reference Genome Efficiency**: Downloads reference genome once and reuses for all batches
- **Merged Output**: Automatically merges expression matrices and samplesheets from all batches

## Features (Common to Both Versions)

- **Automated end-to-end workflow**: From organism name to expression matrices in a single command
- **Flexible reference genome selection**: Use default reference strains or specify custom genome accessions
- **Robust error handling**: Automatic retries and graceful failure management
- **Resume capability**: Continue from any interruption point without re-processing
- **Resource optimization**: Configurable CPU allocation and efficient storage management
- **Docker integration**: No manual dependency installation required
- **Comprehensive quality control**: FastQC and MultiQC reports included

## Prerequisites

- **[Nextflow](https://www.nextflow.io/)** (version 21.04.0 or later)
- **[Docker](https://www.docker.com/)** (version 20.10 or later)

## Installation

1. Clone the MAPPED repository:
```bash
git clone https://github.com/your-org/MAPPED.git
cd MAPPED
```

2. Copy to MAPPED-batch directory (or clone directly):
```bash
cp -r MAPPED MAPPED-batch
cd MAPPED-batch
```

3. Ensure the wrapper script is executable:
```bash
chmod +x run_MAPPED_batch.sh
```

4. Verify Nextflow and Docker are installed:
```bash
nextflow -version
docker --version
```

## Quick Start

Process RNA-seq data for an organism using the batch processing version:

```bash
./run_MAPPED_batch.sh \
    --organism "Escherichia coli" \
    --outdir ./results \
    --workdir ./work \
    --library_layout paired \
    --cpu 48
```

## Usage

### Basic Usage

The `run_MAPPED_batch.sh` wrapper script orchestrates all pipeline modules with automatic batching:

```bash
./run_MAPPED_batch.sh [OPTIONS]
```

### Using a Specific Reference Genome

To use a specific genome assembly instead of the default reference strain:

```bash
./run_MAPPED_batch.sh \
    --organism "Streptomyces coelicolor" \
    --ref-accession GCA_008931305.1 \
    --outdir ./results \
    --workdir ./work \
    --library_layout paired \
    --cpu 24
```

### Example with Large Dataset

For organisms with thousands of samples:

```bash
./run_MAPPED_batch.sh \
    --organism "Salmonella enterica" \
    --outdir ./salmonella_results \
    --workdir ./salmonella_work \
    --library_layout both \
    --cpu 32 \
    --max_concurrent_downloads 10
```

## Pipeline Modules

### 1. Download Metadata (Module 1)
- Queries NCBI SRA for RNA-seq experiments matching the specified organism
- Filters samples based on library layout (single-end, paired-end, or both)
- Generates formatted metadata files for downstream processing

### 2. Download Reference Genome (Module 2)
- Downloads reference genome assemblies from NCBI
- Retrieves genome sequence (FASTA), annotations (GFF), and protein sequences (FAA)
- Supports two modes:
  - **Default mode**: Automatically selects the largest reference genome for the organism
  - **Accession mode**: Downloads a specific genome assembly using its accession number
- **In batch mode**: Downloaded once and reused for all batches

### 3. Download FASTQ (Module 3) - Processed Per Batch
- Downloads raw sequencing data for each batch of 500 samples
- Validates downloaded files
- Creates a samplesheet for downstream analysis
- **In batch mode**: Processes only the samples for the current batch

### 4. Generate Count Matrix (Module 4) - Processed Per Batch
- Performs quality control on raw reads (FastQC)
- Trims adapters and low-quality bases (TrimGalore)
- Quantifies gene expression using Salmon
- Generates normalized expression matrices (TPM and raw counts)
- Creates comprehensive quality reports (MultiQC)
- **In batch mode**: Results are merged across all batches

## Parameters

### Required Parameters

| Parameter | Description | Example |
|-----------|-------------|---------|
| `--organism` | Full taxonomic name of the target organism | `"Escherichia coli"` |
| `--outdir` | Output directory for all results | `/path/to/results` |
| `--workdir` | Nextflow work directory for temporary files | `/path/to/work` |
| `--library_layout` | Sequencing library type: `single`, `paired`, or `both` | `paired` |

### Optional Parameters

| Parameter | Description | Default | Example |
|-----------|-------------|---------|---------|
| `--cpu` | Number of CPUs to allocate per process | System dependent | `16` |
| `--ref-accession` | Specific reference genome accession | Auto-selected | `GCA_008931305.1` |
| `--max_concurrent_downloads` | Maximum number of concurrent FASTQ downloads | `20` | `10` |
| `-h, --help` | Display help message | - | (flag) |

**Note**: The batch version automatically cleans intermediate files after each batch, so there is no `--clean-mode` parameter.

## Output Structure

The batch processing pipeline creates a streamlined output directory:

```
${outdir}/
├── metadata/                    # Downloaded and formatted metadata
│   ├── <organism>_metadata.tsv  # Complete metadata for all samples
│   └── sample_id.csv            # List of SRA accessions
├── samplesheet/                 # Sample information for processing
│   └── samplesheet.csv          # Merged metadata for samples that passed QC
├── ref_genome/                  # Reference genome files (downloaded once)
│   ├── *.fna                    # Genome sequence (FASTA)
│   ├── *.gff                    # Gene annotations (GFF3)
│   ├── *.faa                    # Protein sequences (FASTA)
│   └── datasets_summary.json
├── expression_matrices/         # Final merged expression data
│   ├── counts.csv               # Raw count matrix
│   ├── tpm.csv                  # TPM normalized matrix
│   ├── log_tpm.csv              # Log-transformed TPM
│   └── log_tpm_norm.csv         # Log-transformed normalized TPM
└── batch_logs/                  # Processing logs for each batch
    ├── batch_1.log
    ├── batch_2.log
    └── ...
```

### Batch Processing Workflow

1. **Metadata Download**: Downloads metadata for all samples
2. **Reference Genome**: Downloads reference genome once (reused for all batches)
3. **Batch Processing**: For each batch of 500 samples:
   - Creates batch-specific samplesheet
   - Downloads FASTQ files for the batch
   - Processes through QC, trimming, and quantification
   - Saves expression matrices and samplesheet
   - **Automatically cleans up raw and trimmed files to save space**
4. **Merge Results**: Combines all batch outputs into final expression matrices
5. **Normalization**: Generates normalized log TPM matrix from merged data

### Differences from Standard MAPPED

- **No intermediate files**: Raw FASTQ, trimmed files, and individual QC reports are deleted after each batch
- **Batch logs**: Additional `batch_logs/` directory contains processing logs for each batch
- **Automatic cleanup**: No need for `--clean-mode` flag as cleanup is mandatory
- **Reference genome location**: Stored at top level instead of under `seqFiles/`

## Testing

A test script is provided to verify the pipeline with a small dataset:

```bash
# Activate conda environment first
conda activate nf_core

# Run test
./test_batch_pipeline.sh
```

## Performance Considerations

- Processing time scales linearly with the number of batches
- Disk space usage is kept minimal by cleaning up after each batch
- Memory usage is controlled by processing smaller sample sets
- Network bandwidth is utilized efficiently with concurrent downloads

## Troubleshooting

- Check `batch_logs/` directory for individual batch processing logs
- If a batch fails, you can manually re-run from that batch by modifying the script
- Ensure sufficient disk space for at least one batch (500 samples) plus final outputs