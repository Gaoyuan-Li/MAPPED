//
// Required parameters
//
params.workdir    = null
params.ref_genome = null
params.ref_gff    = null

if ( ! params.workdir )    error "Please provide --workdir"
if ( ! params.ref_genome ) error "Please provide --ref_genome"
if ( ! params.ref_gff )    error "Please provide --ref_gff"

// Process: extract CDS
//
process EXTRACT_CDS {
    tag 'extract_cds'
    container 'quay.io/biocontainers/gffread:0.12.7--h077b44d_6'

    input:
      path genome
      path annotation

    output:
      path 'cds.fa'

    script:
    """
    gffread ${annotation} -g ${genome} -x cds.fa
    """
}

//
// Process: build Salmon index
//
process SALMON_INDEX {
    tag 'salmon_index'
    container 'quay.io/biocontainers/salmon:1.10.3--h45fbf2d_4'

    input:
      path cds_fa

    output:
      path 'salmon_index'

    script:
    """
    salmon index -t ${cds_fa} -i salmon_index --gencode
    """
}

//
// Process: FastQC on raw reads
//
process FASTQC_RAW {
    tag '$sample'
    container 'quay.io/biocontainers/fastqc:0.12.1--hdfd78af_0'
    publishDir "${params.workdir}/fastqc", mode: 'copy'

    input:
      tuple val(sample), path(fq1), path(fq2)

    output:
      path "*.{zip,html}", emit: fastqc_files

    script:
    """
    fastqc -o . ${fq1} ${fq2}
    """
}

//
// Process: TrimGalore
//
process TRIMGALORE {
    tag '$sample'
    container 'quay.io/biocontainers/trim-galore:0.6.9--hdfd78af_0'
    publishDir "${params.workdir}/trimmed", mode: 'copy'

    input:
      tuple val(sample), path(fq1), path(fq2)

    output:
      tuple val(sample),
            path("${sample}_val_1.fq.gz"),
            path("${sample}_val_2.fq.gz")

    script:
    """
    trim_galore --paired --basename ${sample} --output_dir . ${fq1} ${fq2}
    """
}

//
// Process: Salmon quantification
//
process SALMON_QUANT {
    tag '$sample'
    container 'quay.io/biocontainers/salmon:1.10.3--h45fbf2d_4'
    publishDir "${params.workdir}/salmon", mode: 'copy'

    input:
      tuple val(sample), path(fq1), path(fq2)
      path index

    output:
      path "${sample}_quant"

    script:
    """
    salmon quant \
      -i ${index} -l A \
      -1 ${fq1} -2 ${fq2} \
      -o ${sample}_quant \
      --validateMappings
    """
}

//
// Process: merge count matrices
//
process MERGE_COUNTS {
    publishDir "${params.workdir}/expression_matrices", mode: 'copy'

    input:
      path quant_dirs

    output:
      path 'tpm.csv'
      path 'counts.csv'

    script:
    """
    # Extract gene IDs from first quant directory
    firstDir=\$(ls -d *_quant | head -n 1)
    cut -f1 \$firstDir/quant.sf | tail -n +2 > gene_ids.txt

    # Generate TPM matrix
    for d in *_quant; do
      sample=\${d%_quant}
      cut -f4 \$d/quant.sf | tail -n +2 > \${sample}.tpm
    done
    paste gene_ids.txt *.tpm > tmp_tpm.tsv
    {
      printf 'GeneID'
      for d in *_quant; do
        var=\${d%_quant}
        sample=\${var%%_*}
        printf '\t%s' "\${sample}"
      done
      printf '\n'
    } > tpm.csv
    cat tmp_tpm.tsv >> tpm.csv

    # Generate counts matrix
    for d in *_quant; do
      sample=\${d%_quant}
      cut -f5 \$d/quant.sf | tail -n +2 > \${sample}.counts
    done
    paste gene_ids.txt *.counts > tmp_counts.tsv
    {
      printf 'GeneID'
      for d in *_quant; do
        var=\${d%_quant}
        sample=\${var%%_*}
        printf '\t%s' "\${sample}"
      done
      printf '\n'
    } > counts.csv
    cat tmp_counts.tsv >> counts.csv
    """
}

//
// Process: MultiQC report
//
process MULTIQC {
    container 'quay.io/biocontainers/multiqc:1.12--pyhdfd78af_0'
    publishDir "${params.workdir}/multiqc", mode: 'copy'

    input:
      path qc_files

    output:
      path 'multiqc_report.html'

    script:
    """
    multiqc . -o .
    """
}

//
// Main workflow
//
workflow {
    // load samples
    samples_ch = Channel
        .fromPath( file("${params.workdir}/samplesheet/samplesheet.csv") )
        .ifEmpty { error "Sample sheet not found at: ${params.workdir}/samplesheet/samplesheet.csv" }
        .splitCsv(header: true, sep: ',')
        // remove surrounding quotes and normalize keys
        .map { row ->
            row.collectEntries { k, v ->
                def key   = k.trim().toLowerCase().replaceAll(/^"|"$/, '')
                def value = v?.toString()?.trim()?.replaceAll(/^"|"$/, '')
                [ key, value ]
            }
        }
        // keep only rows that have all required fields
        .filter { row ->
            row.sample && row.fastq_1 && row.fastq_2 && row.run_accession
        }
        .ifEmpty {
            error """
    No valid samples found after filtering.
    Make sure your CSV has columns exactly named:
    sample, fastq_1, fastq_2, run_accession.
    """
        }
        // build the tuple for each sample
        .map { row ->
            tuple(
                "${row.sample}_${row.run_accession}",
                file(row.fastq_1),
                file(row.fastq_2)
            )
        }

    // build index
    cds_fa_ch       = EXTRACT_CDS( file(params.ref_genome), file(params.ref_gff) )
    salmon_index_ch = SALMON_INDEX(cds_fa_ch)

    // QC raw
    qc1_ch = FASTQC_RAW(samples_ch)

    // trim
    trimmed_ch = TRIMGALORE(samples_ch)

    // quant
    quant_ch  = SALMON_QUANT(trimmed_ch, salmon_index_ch)

    // merge + report
    count_matrix_ch = MERGE_COUNTS( quant_ch.collect() )
    
    // Collect all QC files for MultiQC
    multiqc_files = qc1_ch.mix(count_matrix_ch).collect()
    MULTIQC( multiqc_files )

}

// Add an onComplete event handler to always delete rotated Nextflow log files
workflow.onComplete {
    def logPattern = ~/\.nextflow\.log\.\d+/  
    new File('.').listFiles().findAll { it.name ==~ logPattern }.each { it.delete() }
}