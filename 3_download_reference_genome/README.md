# Download RNAseq metadata (Nextflow pipeline)

This folder now contains a Nextflow pipeline to download and clean all RNA-seq metadata for an organism from NCBI SRA, using only the scripts in this folder and two Biocontainers images: `quay.io/biocontainers/entrez-direct:22.4--he881be0_0` (for `fetch_metadata`) and `felixlohmeier/pandas:1.3.3` (for `format_metadata`).

## Requirements
- [Nextflow](https://www.nextflow.io/)
- Docker or Singularity/Apptainer (to use the Biocontainers image)

## Usage

```bash
nextflow run main.nf --organism 'Bacillus subtilis' --workdir ../../test_results
```

- `--organism`: Name of the organism (in quotes if it contains spaces)
- `--workdir`: Path to the working/output folder (will be created if it does not exist)
- `--library-layout`: One of `paired`, `single`, or `both` (default), filters metadata to runs with the specified library layout

The output file will be named `{organism_name}_metadata.tsv` (spaces replaced with underscores) and placed in the specified folder.

## Pipeline steps
1. **fetch_metadata**: Uses `esearch` and `efetch` to download SRA run info for the organism.
2. **format_metadata**: Cleans and aggregates the metadata using the provided Python script (`clean_metadata_file.py`).
3. **clean_metadata_tmp**: Prunes the `tmp` directory and old work subdirectories, retaining only the two most recent runs and removing rotated logs.

## Example

```bash
nextflow run main.nf --organism 'Escherichia coli' --workdir ../../test_results
# Output: ../../test_results/Escherichia_coli_metadata.tsv
```

## Notes
- Only the code in this folder is used.
- The output is a cleaned metadata file suitable for downstream analysis.
- The pipeline uses two containers:
  - `quay.io/biocontainers/entrez-direct:22.4--he881be0_0` for metadata fetching.
  - `felixlohmeier/pandas:1.3.3` for metadata formatting.
- A temporary `tmp` directory and older work subdirectories/rotated logs are pruned by `clean_metadata_tmp`.
