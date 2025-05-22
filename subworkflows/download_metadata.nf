#!/usr/bin/env nextflow

// Subworkflow to download and format metadata
workflow DOWNLOAD_METADATA {
    take:
        organism
        library_layout

    main:
        // Fetch metadata from SRA
        raw_metadata = FETCH_METADATA(organism)
        
        // Clean and format metadata
        clean_script = file("${projectDir}/subworkflows/bin/clean_metadata_file.py")
        (cleaned_metadata, sample_ids) = FORMAT_METADATA(
            raw_metadata,
            clean_script,
            organism,
            library_layout
        )

    emit:
        metadata = cleaned_metadata
        samples = sample_ids
}

// Process to fetch metadata from SRA using efetch
process FETCH_METADATA {
    publishDir "${params.outdir}/tmp", mode: 'copy'
    container 'quay.io/biocontainers/entrez-direct:22.4--he881be0_0'

    input:
        val organism
    
    output:
        path 'tmp_metadata.tsv'

    script:
        def query = '"' + organism + '"[Organism] AND "rna seq"[Strategy] AND "transcriptomic"[Source]'
        """
        esearch -db sra -query '${query}' | efetch -db sra -format runinfo > tmp_metadata.tsv
        """
}

// Process to format and clean metadata
process FORMAT_METADATA {
    publishDir "${params.outdir}/metadata", mode: 'copy'
    container 'felixlohmeier/pandas:1.3.3'

    input:
        path raw_tsv
        path clean_script
        val  organism
        val  library_layout

    output:
        path "*_metadata.tsv", emit: metadata
        path "sample_id.csv", emit: sample_ids

    script:
        def safe_name = organism.replaceAll(/\s+/, '_')
        def outfile   = "${safe_name}_metadata.tsv"

        """
        python3 ${clean_script} -i ${raw_tsv} -o ${outfile} -l ${library_layout}
        """
} 