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
    # find largest subdirectory in data
    largest=\$(find tmp/ncbi_dataset/data -mindepth 1 -maxdepth 1 -type d -exec du -s {} + | sort -nr | head -n1 | awk '{print \$2}')
    mkdir -p ref_genome
    cp "\$largest"/*.fna ref_genome/
    cp "\$largest"/*.gff ref_genome/
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