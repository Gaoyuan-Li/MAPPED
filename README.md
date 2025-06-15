# MAPPED - Modular Automated Pipeline for Public Expression Data

MAPPED is a nextflow-based workflow orchestrates four modules for processing public expression data:

1. **download_metadata_efetch**: Download and format metadata from NCBI SRA.
2. **download_fastq**: Download FASTQ files for samples.
3. **download_reference_genome**: Download reference genome from NCBI.
4. **generate_count_matrix**: Perform QC, trimming, quantification, and generate count/tpm matrices.

## Prerequisites

- [Nextflow](https://www.nextflow.io/)
- [Docker](https://www.docker.com/)

## Pipeline Wrapper Script

To run the entire MAPPED pipeline with a single command, use the provided `run_MAPPED.sh` script:

```bash
./run_MAPPED.sh --cpu 20 --organism "Acinetobacter baylyi" --outdir /path/to/output --library_layout paired
```

**Parameters:**

- `--cpu`: Number of threads to be used in the workflow
- `--organism`: Full taxonomic name of the target organism (e.g., "Acinetobacter baylyi") used by metadata and reference genome download modules.
- `--outdir`: Path to the output directory where all results will be saved. If a relative path is provided, it will be converted to an absolute path. The output directory is organized as follows:

```bash
${outdir}/
├── metadata              # Formatted metadata files
├── samplesheet           # CSV samplesheet for FASTQ download
├── seqFiles              # Reference genome FASTA and GFF
├── fastqc                # FastQC reports on raw reads
├── trimmed               # Trimmed FASTQ files
├── salmon                # Salmon quantification results
├── expression_matrices   # TPM and counts matrices
└── multiqc               # MultiQC report
```

- `--library_layout`: Sequencing layout type; use `paired` for paired-end samples only or `single` for single-end samples only, `both` for no filtration. (Note that only `paired` is supported for now)
- `--workdir`: Path to the work directory for the nextflow pipeline.
- `--clean-mode`: Clean up intermediate files and caches after pipeline completion.
