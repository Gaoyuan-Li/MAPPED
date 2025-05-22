# Step 1: Download RNA-seq Metadata Pipeline

This Nextflow DSL2 pipeline downloads and cleans RNA-seq metadata for a specified organism from NCBI SRA, using Entrez Direct and a Pandas-based formatting script.

## Prerequisites

- Nextflow
- Docker

## Usage

```bash
nextflow run main.nf \
  --organism 'Bacillus subtilis' \
  --outdir ../../test_results
```

## Parameters

- `--organism <String>`: Scientific name of the organism (in quotes if containing spaces).
- `--outdir <Path>`: Directory where output will be saved (created if it does not exist).
- `--library-layout <paired|single|both>`: Filter runs by library layout (default: both).

## Outputs

- `{organism_name}_metadata.tsv`: Cleaned metadata file with spaces replaced by underscores.

## Pipeline Steps

1. **fetch_metadata**: Uses `esearch` and `efetch` to retrieve SRA run info.
2. **format_metadata**: Cleans and formats the metadata via `clean_metadata_file.py`.
3. **clean_metadata_tmp**: Removes temporary directories and old logs, keeping the two most recent runs.

## Example

```bash
nextflow run main.nf --organism 'Acinetobacter baylyi' --outdir ../test_results --library_layout paired -resume
# Output: ../test_results/Acinetobacter_baylyi_metadata.tsv and sample_id.csv
```

## Notes

- Only scripts in this folder are used.
- Containers:
  - `quay.io/biocontainers/entrez-direct:22.4--he881be0_0`
  - `felixlohmeier/pandas:1.3.3`
