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
    
    # Check if we got any results
    total_count=\$(grep '"total_count"' summary.json | sed 's/.*"total_count"[[:space:]]*:[[:space:]]*\\([0-9]*\\).*/\\1/')
    
    if [ -z "\$total_count" ] || [ "\$total_count" -eq 0 ]; then
        echo "ERROR: No reference genomes found for ${organism}"
        exit 1
    fi
    
    echo "Found \$total_count reference genome(s)"
    
    # Extract all GCA accessions with their sizes
    # Create a temporary file with GCA accession and size pairs
    > gca_list.txt
    
    # Extract paired_accession and total_sequence_length from each report
    # This handles multiple reports in the JSON
    grep -B50 -A50 '"paired_accession"' summary.json | \
    awk '
        /"paired_accession".*"GCA_/ {
            match(\$0, /"GCA_[^"]+/)
            gca = substr(\$0, RSTART+1, RLENGTH-1)
        }
        /"total_sequence_length"/ && gca {
            match(\$0, /[0-9]+/)
            size = substr(\$0, RSTART, RLENGTH)
            print gca, size
            gca = ""
        }
    ' >> gca_list.txt
    
    # Check if we found any GCA accessions
    if [ ! -s gca_list.txt ]; then
        echo "ERROR: No GenBank (GCA) accessions found in the reference genomes"
        exit 1
    fi
    
    # Sort by size (second column) and get the largest
    selected_gca=\$(sort -k2 -nr gca_list.txt | head -n1 | awk '{print \$1}')
    
    if [ -z "\$selected_gca" ]; then
        echo "ERROR: Failed to select a GenBank accession"
        exit 1
    fi
    
    echo "Selected GenBank accession: \$selected_gca (largest genome)"
    
    # Download the specific GenBank accession
    datasets download genome accession "\$selected_gca" --include gff3,genome --filename ref.zip
    
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
    rm -rf tmp ref.zip summary.json gca_list.txt
    """
}

// Add an onComplete event handler to always delete rotated Nextflow log files
workflow.onComplete {
    def logPattern = ~/\.nextflow\.log\.\d+/  
    new File('.').listFiles().findAll { it.name ==~ logPattern }.each { it.delete() }
}