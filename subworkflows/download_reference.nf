#!/usr/bin/env nextflow

// Subworkflow to download reference genome
workflow DOWNLOAD_REFERENCE {
    take:
        organism
        
    main:
        // Download reference genome from NCBI
        reference_files = DOWNLOAD_REFERENCE_GENOME(organism)
        
    emit:
        reference = reference_files
}

// Process to download reference genome
process DOWNLOAD_REFERENCE_GENOME {
    publishDir "${params.outdir}/seqFiles", mode: 'copy', overwrite: true
    container 'staphb/ncbi-datasets:18.0.2'

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