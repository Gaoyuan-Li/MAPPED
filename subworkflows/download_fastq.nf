#!/usr/bin/env nextflow

// Subworkflow to download FastQ files from SRA
workflow DOWNLOAD_FASTQ {
    main:
        // Make Python scripts in bin folder executable
        "chmod +x ${projectDir}/subworkflows/bin/sra_ids_to_runinfo.py".execute().waitFor()
        "chmod +x ${projectDir}/subworkflows/bin/sra_runinfo_to_ftp.py".execute().waitFor()

        // Define input channel for sample IDs from metadata CSV
        sample_ids = Channel
            .fromPath("${params.outdir}/metadata/sample_id.csv")
            .splitCsv(header:true)
            .map { row -> row.values().first() }

        // Get SRA run information for public database ids
        runinfo = SRA_IDS_TO_RUNINFO(
            sample_ids,
            params.ena_metadata_fields
        )

        // Parse SRA run information, create file containing FTP links
        ftp_links = SRA_RUNINFO_TO_FTP(runinfo.tsv)

        // Process metadata
        sra_metadata = ftp_links.tsv
            .splitCsv(header:true, sep:'\t')
            .map { meta ->
                def meta_clone = meta.clone()
                meta_clone.single_end = meta_clone.single_end.toBoolean()
                return meta_clone
            }
            .unique()

        // Check if we should skip FastQ download
        if (!params.skip_fastq_download) {
            // Branch by download method
            download_reads = sra_metadata
                .branch { meta ->
                    def download_method = 'ftp'
                    if (meta.fastq_aspera && params.download_method == 'aspera') {
                        download_method = 'aspera'
                    }
                    if ((!meta.fastq_aspera && !meta.fastq_1) || params.download_method == 'sratools') {
                        download_method = 'sratools'
                    }

                    ftp: download_method == 'ftp'
                        return [ meta, [ meta.fastq_1, meta.fastq_2 ] ]
                    aspera: download_method == 'aspera'
                        return [ meta, meta.fastq_aspera.tokenize(';').take(2) ]
                    sratools: download_method == 'sratools'
                        return [ meta, meta.run_accession ]
                }

            // Download FastQ files via FTP
            fastq_files = SRA_FASTQ_FTP(download_reads.ftp)

            // Create samplesheet for downloaded files
            SRA_TO_SAMPLESHEET(
                sra_metadata,
                params.nf_core_pipeline,
                params.nf_core_rnaseq_strandedness,
                params.sample_mapping_fields
            )

            // Merge samplesheets across all samples
            samplesheet = SRA_TO_SAMPLESHEET.out.samplesheet
                .map { it[1] }
                .collectFile(name:'tmp_samplesheet.csv', newLine: true, keepHeader: true, sort: { it.baseName })
                .map { it.text.tokenize('\n').join('\n') }
                .collectFile(name:'samplesheet.csv', storeDir: "${params.outdir}/samplesheet")
        }
    
    emit:
        sra_metadata = sra_metadata
}

// Import modules from the 2_download_fastq directory
include { SRA_FASTQ_FTP } from '../2_download_fastq/modules/sra_fastq_ftp/main'
include { SRA_IDS_TO_RUNINFO } from '../2_download_fastq/modules/sra_ids_to_runinfo/main'
include { SRA_RUNINFO_TO_FTP } from '../2_download_fastq/modules/sra_runinfo_to_ftp/main'
include { SRA_TO_SAMPLESHEET } from '../2_download_fastq/modules/sra_to_samplesheet/main' 