# Step 3: Download Reference Genome Pipeline

This Nextflow DSL2 pipeline downloads reference genome and annotation files for a specified organism from NCBI using the `ncbi-datasets` tool.

## Prerequisites

- Nextflow
- Docker

## Usage

```bash
nextflow run main.nf \
  --organism 'Escherichia coli' \
  --outdir ../test_results
```

## Parameters

- `--organism <String>`: Scientific name or taxon ID of the organism (e.g., `Escherichia coli`).
- `--outdir <Path>`: Directory where outputs will be published (workflow creates `<outdir>/seqFiles`).

## Outputs

- `<outdir>/seqFiles/ref_genome/*.fna`: Reference genome FASTA files.
- `<outdir>/seqFiles/ref_genome/*.gff`: Genome annotation GFF3 files.

## Example

```bash
nextflow run main.nf --organism 'Acinetobacter baylyi' --outdir ../test_results
```

## Notes

- Default Docker container: `staphb/ncbi-datasets:18.0.2`.
- Rotated Nextflow log files (`.nextflow.log.<n>`) are removed automatically.
