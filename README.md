# MAPPED - Modular Automated Pipeline for Public Expression Data

MAPPED (Modular Automated Pipeline for Public Expression Data) is a comprehensive Nextflow-based workflow designed to streamline the analysis of public RNA-seq data from NCBI SRA. It automates the entire process from metadata retrieval to expression matrix generation, making large-scale transcriptomics analysis accessible and reproducible.

## Table of Contents
- [Overview](#overview)
- [Features](#features)
- [Prerequisites](#prerequisites)
- [Installation](#installation)
- [Quick Start](#quick-start)
- [Usage](#usage)
- [Pipeline Modules](#pipeline-modules)
- [Parameters](#parameters)
- [Output Structure](#output-structure)
- [Examples](#examples)
- [Troubleshooting](#troubleshooting)
- [Contributing](#contributing)
- [License](#license)

## Overview

MAPPED consists of four integrated modules that work together to process public expression data:

1. **Metadata Download**: Retrieves and formats metadata from NCBI SRA based on organism name
2. **FASTQ Download**: Efficiently downloads sequencing data using optimized protocols
3. **Reference Genome Download**: Obtains reference genome sequences and annotations
4. **Expression Quantification**: Performs quality control, trimming, and gene expression quantification

The pipeline is designed to handle large-scale datasets with built-in error handling, resume capabilities, and resource optimization.

## Features

- **Automated end-to-end workflow**: From organism name to expression matrices in a single command
- **Flexible reference genome selection**: Use default reference strains or specify custom genome accessions
- **Robust error handling**: Automatic retries and graceful failure management
- **Resume capability**: Continue from any interruption point without re-processing
- **Resource optimization**: Configurable CPU allocation and efficient storage management
- **Clean mode**: Automatic cleanup of intermediate files to save disk space
- **Docker integration**: No manual dependency installation required
- **Comprehensive quality control**: FastQC and MultiQC reports included

## Prerequisites

- **[Nextflow](https://www.nextflow.io/)** (version 21.04.0 or later)
- **[Docker](https://www.docker.com/)** (version 20.10 or later)
- **Disk space**: Approximately 5-10 GB per sample (varies by organism and sequencing depth)
- **Memory**: Minimum 8 GB RAM (16 GB or more recommended for large datasets)
- **Internet connection**: Required for downloading data from NCBI

## Installation

1. Clone the MAPPED repository:
```bash
git clone https://github.com/your-org/MAPPED.git
cd MAPPED
```

2. Ensure the wrapper script is executable:
```bash
chmod +x run_MAPPED.sh
```

3. Verify Nextflow and Docker are installed:
```bash
nextflow -version
docker --version
```

## Quick Start

Process RNA-seq data for an organism using the default reference genome:

```bash
./run_MAPPED.sh \
    --organism "Escherichia coli" \
    --outdir ./results \
    --workdir ./work \
    --library_layout paired \
    --cpu 8
```

## Usage

### Basic Usage

The `run_MAPPED.sh` wrapper script orchestrates all pipeline modules:

```bash
./run_MAPPED.sh [OPTIONS]
```

### Using a Specific Reference Genome

To use a specific genome assembly instead of the default reference strain:

```bash
./run_MAPPED.sh \
    --organism "Acinetobacter baylyi" \
    --ref-accession GCA_008931305.1 \
    --outdir ./results \
    --workdir ./work \
    --library_layout paired \
    --cpu 16
```

### Clean Mode

To automatically clean up intermediate files after successful completion:

```bash
./run_MAPPED.sh \
    --organism "Pseudomonas putida" \
    --outdir ./results \
    --workdir ./work \
    --library_layout paired \
    --cpu 16 \
    --clean-mode
```

## Pipeline Modules

### 1. Download Metadata (Module 1)
- Queries NCBI SRA for RNA-seq experiments matching the specified organism
- Filters samples based on library layout (single-end, paired-end, or both)
- Generates formatted metadata files for downstream processing

### 2. Download FASTQ (Module 2)
- Downloads raw sequencing data using optimized protocols (ENA when available, SRA fallback)
- Validates downloaded files
- Creates a samplesheet for downstream analysis

### 3. Download Reference Genome (Module 3)
- Downloads reference genome assemblies from NCBI
- Retrieves genome sequence (FASTA), annotations (GFF), and protein sequences (FAA)
- Supports two modes:
  - **Default mode**: Automatically selects the largest reference genome for the organism
  - **Accession mode**: Downloads a specific genome assembly using its accession number

### 4. Generate Count Matrix (Module 4)
- Performs quality control on raw reads (FastQC)
- Trims adapters and low-quality bases (TrimGalore)
- Quantifies gene expression using Salmon
- Generates normalized expression matrices (TPM and raw counts)
- Creates comprehensive quality reports (MultiQC)

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
| `--clean-mode` | Remove intermediate files after completion | `false` | (flag) |
| `-h, --help` | Display help message | - | (flag) |

## Output Structure

The pipeline creates a well-organized output directory structure:

```
${outdir}/
├── metadata/                 # Downloaded and formatted metadata
│   ├── all_samples.csv      # Complete metadata for all samples
│   └── sample_id.csv        # List of SRA accessions
├── samplesheet/             # Sample information for processing
│   └── samplesheet.csv      # Formatted samplesheet with file paths
├── seqFiles/                # Reference genome files
│   └── ref_genome/
│       ├── *.fna            # Genome sequence (FASTA)
│       ├── *.gff            # Gene annotations (GFF3)
│       ├── *.faa            # Protein sequences (FASTA)
│       └── datasets_summary.json
├── fastqc/                  # Quality control reports
│   ├── *_fastqc.html        # Per-sample QC reports
│   └── *_fastqc.zip         # QC data files
├── trimmed/                 # Adapter-trimmed FASTQ files
│   ├── *_trimmed.fq.gz      # Trimmed sequences
│   └── *_trimming_report.txt
├── salmon/                  # Expression quantification
│   └── <sample_id>/
│       └── quant.sf         # Quantification results
├── expression_matrices/     # Final expression data
│   ├── counts.csv           # Raw count matrix
│   ├── tpm.csv              # TPM normalized matrix
│   ├── log_tpm.csv          # Log-transformed TPM
│   └── log_tpm_norm.csv     # Log-transformed normalized TPM
└── multiqc/                 # Aggregated quality reports
    ├── multiqc_report.html  # Interactive report
    └── multiqc_data/        # Raw MultiQC data
```

### Clean Mode Output

When using `--clean-mode`, only essential outputs are retained:
```
${outdir}/
├── expression_matrices/     # Final expression matrices
├── samplesheet/            # Sample metadata
└── ref_genome/             # Reference genome files
```

## Examples

### Example 1: Basic E. coli Analysis
```bash
./run_MAPPED.sh \
    --organism "Escherichia coli" \
    --outdir ./ecoli_results \
    --workdir ./ecoli_work \
    --library_layout paired \
    --cpu 8
```

### Example 2: Specific Strain Analysis
```bash
# Use a specific P. aeruginosa strain
./run_MAPPED.sh \
    --organism "Pseudomonas aeruginosa" \
    --ref-accession GCA_000006765.1 \
    --outdir ./pa_pa01_results \
    --workdir ./pa_work \
    --library_layout paired \
    --cpu 16
```

### Example 3: Large-Scale Analysis with Cleanup
```bash
# Process all available samples and clean up afterwards
./run_MAPPED.sh \
    --organism "Salmonella enterica" \
    --outdir ./salmonella_results \
    --workdir ./salmonella_work \
    --library_layout both \
    --cpu 32 \
    --clean-mode
```

### Example 4: Resume an Interrupted Run
```bash
# Simply re-run the same command - Nextflow will resume automatically
./run_MAPPED.sh \
    --organism "Mycobacterium tuberculosis" \
    --outdir ./mtb_results \
    --workdir ./mtb_work \
    --library_layout paired \
    --cpu 16
```

## Troubleshooting

### Common Issues

1. **"No reference genomes found" error**
   - Verify the organism name matches NCBI taxonomy exactly
   - Check if the organism has reference genomes available in NCBI

2. **Docker permission errors**
   - Ensure your user is in the docker group: `sudo usermod -aG docker $USER`
   - Log out and back in for changes to take effect

3. **Disk space issues**
   - Use `--clean-mode` to reduce storage requirements
   - Ensure sufficient space in both `--outdir` and `--workdir`

4. **Memory errors during Salmon quantification**
   - Reduce the number of CPUs with `--cpu` parameter
   - Close other memory-intensive applications

5. **Network timeouts during download**
   - The pipeline will automatically retry failed downloads
   - For persistent issues, check your internet connection and firewall settings

### Debug Mode

For detailed debugging information, run individual modules:
```bash
# Run only the metadata download module
cd 1_download_metadata_efetch
nextflow run main.nf -with-trace -with-report
```

## Contributing

We welcome contributions! Please see our [Contributing Guidelines](CONTRIBUTING.md) for details on:
- Reporting bugs
- Suggesting enhancements
- Submitting pull requests
- Code style guidelines

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Citation

If you use MAPPED in your research, please cite:
```
[Citation information to be added]
```

## Contact

For questions, issues, or suggestions:
- Open an issue on [GitHub](https://github.com/your-org/MAPPED/issues)
- Email: [contact information]

---

**Note**: MAPPED is under active development. Please check for updates regularly and report any issues you encounter.