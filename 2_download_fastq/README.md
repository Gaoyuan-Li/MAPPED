# Step 2: Download FASTQ Files Pipeline

This Nextflow DSL2 pipeline downloads FASTQ files for specified SRA accessions and generates a consolidated samplesheet.

## Prerequisites

- Nextflow
- Docker

## Usage

```bash
nextflow run main.nf \
  --workdir /path/to/your/workdir
```

## Parameters

- `--workdir <Path>`: Directory containing `metadata/sample_id.csv`.
- `--outdir <Path>`: (Optional) Override default output directory (defaults to `<workdir>/metadata`).
- `--ena_metadata_fields <Fields>`: (Optional) Additional ENA metadata fields to fetch.
- `--skip_fastq_download <true|false>`: (Optional) Fetch metadata only without downloading FASTQ files if set to true.

## Inputs

- `metadata/sample_id.csv`: CSV file with header `id` and a column of SRA accession IDs. Example:

  ```csv
  id
  SRR123456
  SRR234567
  ```

## Outputs

- `metadata/runinfo_ftp.tsv`: FTP links for FASTQ files.
- `metadata/fastq/`: Downloaded FASTQ files.
- `metadata/samplesheet/samplesheet.csv`: Samplesheet for downstream pipelines.
- `metadata/samplesheet/id_mappings.csv`: Sample ID mapping file.

## Example

```bash
nextflow run main.nf --outdir ../test_results -resume
```
