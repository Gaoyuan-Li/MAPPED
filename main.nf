#!/usr/bin/env nextflow

// Enable DSL2
nextflow.enable.dsl = 2

// Define parameters
params.organism = null
params.outdir = null
params.library_layout = 'paired'

// Check required parameters
if (!params.organism) error "Please provide --organism parameter"
if (!params.outdir) error "Please provide --outdir parameter"

// Include sub-workflows
include { VERIFY_DEPENDENCIES } from './subworkflows/verify_dependencies'
include { DOWNLOAD_METADATA } from './subworkflows/download_metadata'
include { DOWNLOAD_FASTQ } from './subworkflows/download_fastq'
include { DOWNLOAD_REFERENCE } from './subworkflows/download_reference'
include { GENERATE_COUNT_MATRIX } from './subworkflows/generate_count_matrix'

// Main workflow
workflow {
    // First, verify that all dependencies are present
    VERIFY_DEPENDENCIES()
    
    log.info """
    =======================================
    MAPPED WORKFLOW
    =======================================
    Organism:        ${params.organism}
    Output dir:      ${params.outdir}
    Library layout:  ${params.library_layout}
    =======================================
    """
    
    // Step 1: Download metadata
    DOWNLOAD_METADATA(params.organism, params.library_layout)
    
    // Step 2: Download FastQ files
    DOWNLOAD_FASTQ()
    
    // Step 3: Download reference genome
    DOWNLOAD_REFERENCE(params.organism)
    
    // Step 4: Generate count matrix
    GENERATE_COUNT_MATRIX()
}

// Clean up rotated log files on completion
workflow.onComplete {
    log.info "Pipeline completed at: ${new Date()}"
    log.info "Results have been saved to: ${params.outdir}"
    def logPattern = ~/\.nextflow\.log\.\d+/
    new File('.').listFiles().findAll { it.name ==~ logPattern }.each { it.delete() }
}
