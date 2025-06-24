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

    script:
    """
    datasets download genome taxon '${organism}' --reference --include genome,gff3 --filename ref.zip
    unzip ref.zip -d tmp
    # find largest GenBank assembly (GCA) subdirectory in data
    largest=\$(find tmp/ncbi_dataset/data -mindepth 1 -maxdepth 1 -type d -name "GCA_*" -exec du -s {} + | sort -nr | head -n1 | awk '{print \$2}')
    
    if [ -z "\$largest" ]; then
        echo "Error: No GenBank assembly (GCA) found in downloaded data"
        exit 1
    fi
    
    echo "Selected GenBank assembly: \$(basename \$largest)"
    
    mkdir -p ref_genome
    # Copy fna files - these should begin with GCA
    for fna in "\$largest"/*_genomic.fna "\$largest"/*.fna; do
        if [ -f "\$fna" ]; then
            basename=\$(basename "\$fna")
            # Prefer GCA-prefixed files
            if [[ "\$basename" =~ ^GCA_ ]]; then
                cp "\$fna" "ref_genome/\$basename"
                break  # Use the first GCA-prefixed fna file found
            fi
        fi
    done
    
    # If no GCA-prefixed fna found, copy the first available fna with GCA prefix added
    if [ ! -f ref_genome/*.fna ]; then
        for fna in "\$largest"/*.fna; do
            if [ -f "\$fna" ]; then
                assembly=\$(basename "\$largest")
                basename=\$(basename "\$fna")
                cp "\$fna" "ref_genome/\${assembly}_\${basename}"
                break
            fi
        done
    fi
    
    # Copy genomic.gff file (doesn't need GCA prefix)
    for gff in "\$largest"/*genomic.gff "\$largest"/*.gff; do
        if [ -f "\$gff" ]; then
            basename=\$(basename "\$gff")
            cp "\$gff" "ref_genome/\$basename"
            break  # Use the first gff file found (preferably genomic.gff)
        fi
    done
    
    # ensure output files are world-readable for publishDir
    chmod a+r ref_genome/*.fna ref_genome/*.gff
    # cleanup
    rm -rf tmp ref.zip
    """
}

// Add an onComplete event handler to always delete rotated Nextflow log files
workflow.onComplete {
    def logPattern = ~/\.nextflow\.log\.\d+/  
    new File('.').listFiles().findAll { it.name ==~ logPattern }.each { it.delete() }
}