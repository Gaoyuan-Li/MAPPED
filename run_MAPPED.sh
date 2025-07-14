#!/usr/bin/env bash
set -euo pipefail

function usage() {
  cat <<EOF
Usage: $0 --organism ORGANISM --outdir OUTDIR --library_layout LIB_LAYOUT --workdir WORKDIR --clean-mode CLEAN_MODE --cpu CPU [--ref-accession REF_ACCESSION] [--max_concurrent_downloads N]

Options:
  --organism        Organism name (e.g., "Acinetobacter baylyi") - required for metadata download
  --outdir          Output directory for pipeline results
  --workdir         Work directory for Nextflow 'work' files
  --library_layout  Library layout: 'single', 'paired', or 'both'
  --clean-mode      Clean up intermediate files and caches after pipeline completion.
  --cpu             Number of CPUs to allocate per process
  --ref-accession   Optional: specific reference genome accession (e.g., "GCA_008931305.1"). 
                    If not provided, automatically selects the reference strain for the organism.
  --max_concurrent_downloads  Optional: Maximum number of concurrent downloads (default: 20)
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
REF_ACCESSION=""
MAX_CONCURRENT_DOWNLOADS=""

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
    --ref-accession)
      REF_ACCESSION="$2"
      shift 2
      ;;
    --max_concurrent_downloads)
      MAX_CONCURRENT_DOWNLOADS="$2"
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

# Validate library_layout parameter
if [[ "$LIB_LAYOUT" != "single" && "$LIB_LAYOUT" != "paired" && "$LIB_LAYOUT" != "both" ]]; then
  echo "Error: Invalid library_layout value: $LIB_LAYOUT"
  echo "Valid values are: single, paired, both"
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
pushd 1_download_metadata_efetch > /dev/null 2>&1
nextflow run main.nf -work-dir "$WORKDIR" --organism "$ORGANISM" --outdir "$OUTDIR" --library_layout "$LIB_LAYOUT" -resume
popd > /dev/null 2>&1

# Step 2: Download FASTQ
echo "=== Step 2: Download FASTQ ==="
pushd 2_download_fastq > /dev/null 2>&1
nextflow run main.nf -work-dir "$WORKDIR" --outdir "$OUTDIR" ${MAX_CONCURRENT_DOWNLOADS:+--max_concurrent_downloads $MAX_CONCURRENT_DOWNLOADS} -resume
popd > /dev/null 2>&1

# Step 3: Download reference genome
echo "=== Step 3: Download reference genome ==="
pushd 3_download_reference_genome > /dev/null 2>&1
if [[ -n "$REF_ACCESSION" ]]; then
  nextflow run main.nf -work-dir "$WORKDIR" --ref_accession "$REF_ACCESSION" --outdir "$OUTDIR" ${CPU:+--cpu $CPU} -resume
else
  nextflow run main.nf -work-dir "$WORKDIR" --organism "$ORGANISM" --outdir "$OUTDIR" ${CPU:+--cpu $CPU} -resume
fi
popd > /dev/null 2>&1

# Step 4: Generate count/tpm matrix
echo "=== Step 4: Generate count/tpm matrix ==="
pushd 4_generate_count_matrix > /dev/null 2>&1
nextflow run main.nf -work-dir "$WORKDIR" --outdir "$OUTDIR" ${CPU:+--cpu $CPU} -resume
popd > /dev/null 2>&1

# Print sample counts after Step 4
echo "=== Sample Count Summary ==="
if [[ -f "$OUTDIR/samplesheet/samplesheet_download.csv" ]]; then
  download_count=$(tail -n +2 "$OUTDIR/samplesheet/samplesheet_download.csv" | grep -c '^')
  echo "Downloaded samples (samplesheet_download.csv): $download_count"
else
  echo "samplesheet_download.csv not found"
fi

if [[ -f "$OUTDIR/samplesheet/samplesheet.csv" ]]; then
  filtered_count=$(tail -n +2 "$OUTDIR/samplesheet/samplesheet.csv" | grep -c '^')
  echo "Samples passing filtration (samplesheet.csv): $filtered_count"
else
  echo "samplesheet.csv not found"
fi
echo "============================="

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