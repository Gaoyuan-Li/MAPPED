#!/usr/bin/env bash
set -euo pipefail

function usage() {
  cat <<EOF
Usage: $0 --organism ORGANISM --outdir OUTDIR --library_layout LIB_LAYOUT

Options:
  --organism        Organism name (e.g., "Acinetobacter baylyi")
  --outdir          Output directory for pipeline results
  --library_layout  Library layout (e.g., paired or single)
  -h, --help        Show this help message and exit
EOF
}

# Parse arguments
ORGANISM=""
OUTDIR=""
LIB_LAYOUT=""

while [[ $# -gt 0 ]]; do
  key="$1"
  case $key in
    --organism)
      ORGANISM="$2"
      shift 2
      ;;
    --outdir)
      OUTDIR="$2"
      shift 2
      ;;
    --library_layout)
      LIB_LAYOUT="$2"
      shift 2
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    *)
      echo "Unknown option: $1"
      usage
      exit 1
      ;;
  esac
done

# Check required arguments
if [[ -z "$ORGANISM" || -z "$OUTDIR" || -z "$LIB_LAYOUT" ]]; then
  echo "Error: Missing required arguments."
  usage
  exit 1
fi

# Convert OUTDIR to an absolute path and ensure it exists
if [[ "$OUTDIR" != /* ]]; then
  OUTDIR="$(pwd)/$OUTDIR"
fi
mkdir -p "$OUTDIR"

# Step 1: Download metadata
echo "=== Step 1: Download metadata ==="
pushd 1_download_metadata_efetch > /dev/null
nextflow run main.nf --organism "$ORGANISM" --outdir "$OUTDIR" --library_layout "$LIB_LAYOUT" -resume
popd > /dev/null

# Step 2: Download FASTQ
echo "=== Step 2: Download FASTQ ==="
pushd 2_download_fastq > /dev/null
nextflow run main.nf --outdir "$OUTDIR" -resume
popd > /dev/null

# Step 3: Download reference genome
echo "=== Step 3: Download reference genome ==="
pushd 3_download_reference_genome > /dev/null
nextflow run main.nf --organism "$ORGANISM" --outdir "$OUTDIR" -resume
popd > /dev/null

# Step 4: Generate count/tpm matrix
echo "=== Step 4: Generate count/tpm matrix ==="
pushd 4_generate_count_matrix > /dev/null
nextflow run main.nf --outdir "$OUTDIR" -resume
popd > /dev/null

echo "\nAll steps completed successfully!" 