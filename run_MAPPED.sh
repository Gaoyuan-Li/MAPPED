#!/usr/bin/env bash
set -euo pipefail

function usage() {
  cat <<EOF
Usage: $0 --organism ORGANISM --outdir OUTDIR --library_layout LIB_LAYOUT --workdir WORKDIR --clean-mode CLEAN_MODE --cpu CPU

Options:
  --organism        Organism name (e.g., "Acinetobacter baylyi")
  --outdir          Output directory for pipeline results
  --workdir         Work directory for Nextflow 'work' files
  --library_layout  Library layout (e.g., paired or single)
  --clean-mode      Clean up intermediate files and caches after pipeline completion.
  --cpu             Number of CPUs to allocate per process
  -h, --help        Show this help message and exit
EOF
}

# Parse arguments
ORGANISM=""
OUTDIR=""
LIB_LAYOUT=""
CLEAN_MODE="false"
CPU=""
WORKDIR=""

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
    --workdir)
      WORKDIR="$2"
      shift 2
      ;;
    --clean-mode)
      CLEAN_MODE="true"
      shift
      ;;
    --cpu)
      CPU="$2"
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
if [[ -z "$ORGANISM" || -z "$OUTDIR" || -z "$LIB_LAYOUT" || -z "$WORKDIR" ]]; then
  echo "Error: Missing required arguments."
  usage
  exit 1
fi

# Convert OUTDIR to an absolute path and ensure it exists
if [[ "$OUTDIR" != /* ]]; then
  OUTDIR="$(pwd)/$OUTDIR"
fi
mkdir -p "$OUTDIR"

# Convert WORKDIR to an absolute path and ensure it exists
if [[ "$WORKDIR" != /* ]]; then
  WORKDIR="$(pwd)/$WORKDIR"
fi
mkdir -p "$WORKDIR"

# Step 1: Download metadata
echo "=== Step 1: Download metadata ==="
pushd 1_download_metadata_efetch > /dev/null
nextflow run main.nf -work-dir "$WORKDIR" --organism "$ORGANISM" --outdir "$OUTDIR" --library_layout "$LIB_LAYOUT" -resume
popd > /dev/null

# Step 2: Download FASTQ
echo "=== Step 2: Download FASTQ ==="
pushd 2_download_fastq > /dev/null
nextflow run main.nf -work-dir "$WORKDIR" --outdir "$OUTDIR" -resume
popd > /dev/null

# Step 3: Download reference genome
echo "=== Step 3: Download reference genome ==="
pushd 3_download_reference_genome > /dev/null
nextflow run main.nf -work-dir "$WORKDIR" --organism "$ORGANISM" --outdir "$OUTDIR" ${CPU:+--cpu $CPU} -resume
popd > /dev/null

# Step 4: Generate count/tpm matrix
echo "=== Step 4: Generate count/tpm matrix ==="
pushd 4_generate_count_matrix > /dev/null
nextflow run main.nf -work-dir "$WORKDIR" --outdir "$OUTDIR" ${CPU:+--cpu $CPU} -resume
popd > /dev/null

echo "All steps completed successfully!"

if [[ "$CLEAN_MODE" == "true" ]]; then
  echo "=== Clean mode enabled: cleaning intermediate files ==="
  
  # Preserve ref_genome folder by moving it to a temporary location
  if [[ -d "$OUTDIR/seqFiles/ref_genome" ]]; then
    echo "Preserving ref_genome folder..."
    mv "$OUTDIR/seqFiles/ref_genome" "$OUTDIR/ref_genome_temp"
  fi
  
  # Delete everything in OUTDIR except expression_matrices and samplesheet
  find "$OUTDIR" -mindepth 1 -maxdepth 1 ! -name expression_matrices ! -name samplesheet ! -name ref_genome_temp -exec rm -rf {} +
  
  # Move ref_genome back to the same level as expression_matrices and samplesheet
  if [[ -d "$OUTDIR/ref_genome_temp" ]]; then
    echo "Moving ref_genome to final location..."
    mv "$OUTDIR/ref_genome_temp" "$OUTDIR/ref_genome"
  fi
  
  # Delete work, .nextflow, and .nextflow.log in module directories
  for sub in 1_download_metadata_efetch 2_download_fastq 3_download_reference_genome 4_generate_count_matrix; do
    rm -rf "$sub/work" "$sub/.nextflow" "$sub/.nextflow.log"
  done
  
  # Clean the global Nextflow work directory
  rm -rf "$WORKDIR"
fi 