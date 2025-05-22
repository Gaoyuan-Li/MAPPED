#!/usr/bin/env nextflow

def runStartTime = System.currentTimeMillis()

params.organism = null
params.workdir = null

process fetch_metadata {

    publishDir 'tmp', mode: 'copy'

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

process format_metadata {

    publishDir "${params.workdir}/metadata", mode: 'copy'   // ⬅ copy results into metadata subfolder of workdir

    container 'felixlohmeier/pandas:1.3.3'

    input:
        path raw_tsv
        path clean_script
        val  organism
        val  library_layout

    output:
        path "*_metadata.tsv"                 // ⬅ any metadata file
        path "sample_id.csv"                  // ⬅ sample IDs for downstream use

    script:
        def safe_name = organism.replaceAll(/\s+/, '_')
        def outfile   = "${safe_name}_metadata.tsv"

        """
        python3 ${clean_script} -i ${raw_tsv} -o ${outfile} -l ${library_layout}
        """
}

// Add a dedicated cleanup process to prune tmp and old work dirs
process clean_metadata_tmp {
    cache false
    input:
        path metadata_tsv
    script:
    """
    rm -rf ${projectDir}/tmp
    # remove old Nextflow log rotations, keep only .nextflow.log
    rm -f ${projectDir}/.nextflow.log.[0-9]* || true
    """
}

workflow {
    if ( !params.organism || !params.workdir ) {
        error "You must provide both --organism and --workdir parameters."
    }

    raw_metadata = fetch_metadata( params.organism )

    clean_script = file( 'bin/clean_metadata_file.py' )

    ( cleaned_metadata, sample_ids ) = format_metadata(
        raw_metadata,
        clean_script,
        params.organism,
        params.library_layout
    )

    clean_metadata_tmp(cleaned_metadata)
}
