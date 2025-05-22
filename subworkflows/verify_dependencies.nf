#!/usr/bin/env nextflow

// Workflow to verify that all necessary files and dependencies are present
workflow VERIFY_DEPENDENCIES {
    main:
        // 1. Check that necessary Python scripts exist
        def required_scripts = [
            "${projectDir}/subworkflows/bin/clean_metadata_file.py",
            "${projectDir}/subworkflows/bin/sra_ids_to_runinfo.py", 
            "${projectDir}/subworkflows/bin/sra_runinfo_to_ftp.py"
        ]
        
        for (script in required_scripts) {
            if (!file(script).exists()) {
                error "Required script not found: ${script}"
            }
        }
        
        // 2. Check that required modules exist
        def required_modules = [
            "${projectDir}/2_download_fastq/modules/sra_fastq_ftp",
            "${projectDir}/2_download_fastq/modules/sra_ids_to_runinfo", 
            "${projectDir}/2_download_fastq/modules/sra_runinfo_to_ftp",
            "${projectDir}/2_download_fastq/modules/sra_to_samplesheet"
        ]
        
        for (module in required_modules) {
            if (!file(module).exists()) {
                error "Required module not found: ${module}"
            }
        }
        
        // 3. Check that required parameters are set
        if (!params.organism) {
            error "Missing required parameter: --organism"
        }
        if (!params.outdir) {
            error "Missing required parameter: --outdir"
        }
        
        log.info "All dependencies verified successfully"
} 