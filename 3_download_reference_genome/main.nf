#!/usr/bin/env nextflow

// Workflow to download a reference genome from NCBI and save to outdir

workflow {
    if (!params.organism) {
        error "Missing required parameter --organism"
    }

    DOWNLOAD_REFERENCE(params.organism)
}

process DOWNLOAD_REFERENCE {
    publishDir "${params.outdir}/seqFiles", mode: 'copy', overwrite: true

    input:
    val organism

    output:
    path 'ref_genome/*.fna'
    path 'ref_genome/*.gff'
    path 'ref_genome/datasets_summary.json'

    script:
    """
    # First get the summary to find GenBank accession ID
    datasets summary genome taxon '${organism}' --reference --assembly-source refseq > summary.json
    
    # Parse the JSON to get GCA accession(s)
    # Extract all GCA accessions with their sizes
    gca_info=\$(python3 -c "
import json
import sys

with open('summary.json', 'r') as f:
    data = json.load(f)

if 'reports' not in data or len(data['reports']) == 0:
    print('ERROR: No reference genomes found', file=sys.stderr)
    sys.exit(1)

# Collect GenBank accessions with their total sequence lengths
gca_accessions = []
for report in data['reports']:
    if 'paired_accession' in report and report['paired_accession'].startswith('GCA_'):
        gca = report['paired_accession']
        size = int(report.get('assembly_stats', {}).get('total_sequence_length', 0))
        gca_accessions.append((gca, size))

if not gca_accessions:
    print('ERROR: No GenBank (GCA) accessions found', file=sys.stderr)
    sys.exit(1)

# Sort by size and pick the largest
gca_accessions.sort(key=lambda x: x[1], reverse=True)
selected_gca = gca_accessions[0][0]

print(selected_gca)
")
    
    if [ "\$?" -ne 0 ]; then
        echo "Error parsing datasets summary"
        exit 1
    fi
    
    echo "Selected GenBank accession: \$gca_info"
    
    # Download the specific GenBank accession
    datasets download genome accession "\$gca_info" --include gff3,genome --filename ref.zip
    
    # Extract and organize files
    unzip ref.zip -d tmp
    
    # Find the GCA directory
    gca_dir=\$(find tmp/ncbi_dataset/data -mindepth 1 -maxdepth 1 -type d -name "GCA_*" | head -n1)
    
    if [ -z "\$gca_dir" ]; then
        echo "Error: GenBank assembly directory not found after download"
        exit 1
    fi
    
    mkdir -p ref_genome
    
    # Copy fna files
    for fna in "\$gca_dir"/*_genomic.fna "\$gca_dir"/*.fna; do
        if [ -f "\$fna" ]; then
            cp "\$fna" "ref_genome/\$(basename "\$fna")"
            break
        fi
    done
    
    # Copy gff files
    for gff in "\$gca_dir"/*genomic.gff "\$gca_dir"/*.gff; do
        if [ -f "\$gff" ]; then
            cp "\$gff" "ref_genome/\$(basename "\$gff")"
            break
        fi
    done
    
    # Save the datasets summary
    cp summary.json ref_genome/datasets_summary.json
    
    # Ensure output files are world-readable for publishDir
    chmod a+r ref_genome/*
    
    # Cleanup
    rm -rf tmp ref.zip summary.json
    """
}

// Add an onComplete event handler to always delete rotated Nextflow log files
workflow.onComplete {
    def logPattern = ~/\.nextflow\.log\.\d+/  
    new File('.').listFiles().findAll { it.name ==~ logPattern }.each { it.delete() }
}