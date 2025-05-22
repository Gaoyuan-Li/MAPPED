# Download FASTQ Files with Nextflow

This pipeline downloads FASTQ files for SRA accessions specified in a CSV file and generates a consolidated sample sheet.

## Prerequisites

- Nextflow (version >= 21.04)
- Docker

## Input

Place your `sample_id.csv` file in the `metadata` subdirectory of your work directory. The CSV should have a header and a single column of SRA accession IDs. For example:

    metadata/sample_id.csv
    ------------------------
    id
    SRR123456
    SRR234567

## Usage

Run the pipeline from this directory:

```bash
nextflow run main.nf \
  --workdir /path/to/your/workdir
```

By default, output files (fastq, MD5 checksums, runinfo TSVs, and the final samplesheet) will be stored in `/path/to/your/workdir/metadata`.

## Parameters

- `--workdir`: Path to your working directory containing `metadata/sample_id.csv`.
- `--outdir`: (Optional) Override the default output directory (defaults to `<workdir>/metadata`).
- `--ena_metadata_fields`: (Optional) Additional ENA metadata fields to fetch.
- `--skip_fastq_download`: (Optional) Set to `true` to fetch metadata only without downloading FASTQ files.

## Outputs

- `metadata/runinfo_ftp.tsv`: FTP links for FASTQ files.
- `metadata/fastq/`: Downloaded FASTQ files.
- `metadata/samplesheet/samplesheet.csv`: Consolidated samplesheet for downstream nf-core pipelines.
- `metadata/samplesheet/id_mappings.csv`: Sample ID mappings file.
