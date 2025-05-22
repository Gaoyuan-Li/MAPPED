# Step 3: Download Reference Genome

This Nextflow workflow downloads the reference genome and annotation files for a specified organism from NCBI using the `ncbi-datasets` command-line tool, then publishes the results to the desired output directory.

## Prerequisites

- Nextflow
- Docker
- Internet connection for accessing NCBI datasets

## Configuration

The pipeline uses DSL2 and is configured to run inside the `staphb/ncbi-datasets:18.0.2` Docker container by default (see `nextflow.config`).

## Parameters

- `--organism <String>` (required)
  - The scientific name or taxon ID of the organism to download (e.g., `Escherichia coli`).
- `--outdir <Path>` (required)
  - Path to the directory where output files will be published. The workflow will create a subdirectory named `seqFiles` under this path.

## Usage

```bash
nextflow run main.nf \
  --organism 'Escherichia coli' \
  --outdir ../test_results
```

## Outputs

Upon successful completion, the following files will be available under `<outdir>/seqFiles/ref_genome/`:

- `*.fna` &nbsp;&nbsp;– FASTA file(s) containing the reference genome sequence.
- `*.gff` &nbsp;&nbsp;– GFF3 file(s) containing genome annotations.

All files are copied in world-readable mode to ensure they can be accessed by downstream steps.

## Cleaning Up

An `onComplete` handler automatically removes any rotated Nextflow log files (`.nextflow.log.<n>`). To clean up temporary files created during the download process, simply delete any remaining `tmp/` directories or `.zip` archives if present.
