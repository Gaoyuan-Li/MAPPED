//
// Required parameters
//
params.outdir    = null

if ( ! params.outdir )    error "Please provide --outdir"

// Resource configuration: threads per task and max parallel runs
params.cpu = params.cpu ?: 20
def threads_per_task = 4
def max_parallel = (params.cpu / threads_per_task) as Integer

// Auto-detect reference genome and GFF in seqFiles/ref_genome under outdir
def refDir = new File("${params.outdir}/seqFiles/ref_genome")
if ( ! refDir.exists() ) error "Reference genome directory not found: ${refDir}"
def fastaFiles = refDir.list().findAll { it.endsWith('.fna') || it.endsWith('.fa') }
if ( fastaFiles.size() == 0 ) error "No FASTA (.fna/.fa) file found in ${refDir}"
if ( fastaFiles.size() > 1 ) error "Multiple FASTA files found in ${refDir}: ${fastaFiles}"
params.ref_genome = "${refDir}/${fastaFiles[0]}"
def gffFiles = refDir.list().findAll { it.endsWith('.gff') }
if ( gffFiles.size() == 0 ) error "No GFF (.gff) file found in ${refDir}"
if ( gffFiles.size() > 1 ) error "Multiple GFF files found in ${refDir}: ${gffFiles}"
params.ref_gff = "${refDir}/${gffFiles[0]}"

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
process FASTQC {
    tag '$sample'
    container 'quay.io/biocontainers/fastqc:0.12.1--hdfd78af_0'
    publishDir "${params.outdir}/fastqc", mode: 'copy'
    cpus threads_per_task
    maxForks max_parallel

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
    publishDir "${params.outdir}/trimmed", mode: 'copy'
    cpus threads_per_task
    maxForks max_parallel

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
    publishDir "${params.outdir}/salmon", mode: 'copy'
    cpus threads_per_task
    maxForks max_parallel

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
    publishDir "${params.outdir}/expression_matrices", mode: 'copy'

    input:
      path quant_dirs

    output:
      path 'tpm.tsv'
      path 'counts.tsv'

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
    } > tpm.tsv
    cat tmp_tpm.tsv >> tpm.tsv

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
    } > counts.tsv
    cat tmp_counts.tsv >> counts.tsv
    """
}

//
// Process: MultiQC report
//
process MULTIQC {
    container 'quay.io/biocontainers/multiqc:1.12--pyhdfd78af_0'
    publishDir "${params.outdir}/multiqc", mode: 'copy'

    input:
      path qc_files

    output:
      path 'multiqc_report.html', emit: html
      path 'multiqc_data.json', emit: json

    script:
    """
    multiqc . -o . --json multiqc_data.json
    """
}

//
// Process: parse MultiQC JSON to extract samples passing key metrics
process PARSE_QC {
    tag 'parse_multiqc'
    container 'python:3.9-slim'

    input:
      path multiqc_json

    output:
      path 'passed_samples.txt', emit: passlist

    script:
    """
    python3 - << 'EOF'
    import json
    from pathlib import Path
    data = json.load(open('multiqc_data.json'))
    # Navigate to FastQC module results
    fastqc = data.get('modules', {}).get('fastqc', {})
    results = fastqc.get('plot_data', fastqc.get('results', fastqc.get('stats', {})))
    passed = []
    for sample, metrics in results.items():
        ok = True
        for m in ['per_base_sequence_quality', 'per_sequence_quality_scores', 'per_base_n_content']:
            if metrics.get(m, {}).get('status') != 'pass':
                ok = False
                break
        if ok:
            passed.append(sample)
    open('passed_samples.txt', 'w').write('\n'.join(passed))
    EOF
    """
}

// Process: filter sample sheet CSV based on passed samples
process FILTER_SAMPLESHEET {
    tag 'filter_samplesheet'
    container 'alpine:latest'
    publishDir "${params.outdir}/samplesheet", mode: 'copy', overwrite: true

    input:
      path samplesheet
      path passedlist

    output:
      path 'samplesheet.csv'

    script:
    """
    cp ${samplesheet} tmp.csv
    head -n1 tmp.csv > samplesheet.csv
    grep -F -f ${passedlist} tmp.csv >> samplesheet.csv
    """
}

//
// Main workflow
//
workflow {
    // load samples
    samples_ch = Channel
        .fromPath( file("${params.outdir}/samplesheet/samplesheet.csv") )
        .ifEmpty { error "Sample sheet not found at: ${params.outdir}/samplesheet/samplesheet.csv" }
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
                file("${params.outdir}/${row.fastq_1}"),
                file("${params.outdir}/${row.fastq_2}")
            )
        }

    // build index
    cds_fa_ch       = EXTRACT_CDS( file(params.ref_genome), file(params.ref_gff) )
    salmon_index_ch = SALMON_INDEX(cds_fa_ch)

    // trim
    trimmed_ch = TRIMGALORE(samples_ch)

    // QC trimmed reads
    qc_ch = FASTQC(trimmed_ch)

    // MultiQC on trimmed QC results
    multiqc = MULTIQC( qc_ch.collect() )
    multiqc_json_ch = multiqc.json
    // Parse MultiQC JSON for passed samples
    passed_ch = PARSE_QC( multiqc_json_ch )

    // Filter original sample sheet based on passed samples
    FILTER_SAMPLESHEET( file("${params.outdir}/samplesheet/samplesheet.csv"), passed_ch )

    // Filter trimmed reads and perform Salmon quantification only on passed samples
    pass_tuple_ch      = passed_ch.map { id -> tuple(id, id) }
    matched_trim_ch    = trimmed_ch.join(pass_tuple_ch)
    filtered_trimmed_ch = matched_trim_ch.map { key, trim_val, _ -> trim_val }
    quant_ch           = SALMON_QUANT(filtered_trimmed_ch, salmon_index_ch)

    // merge count matrices
    count_matrix_ch = MERGE_COUNTS( quant_ch.collect() )
}

// Add an onComplete event handler to always delete rotated Nextflow log files
workflow.onComplete {
    def logPattern = ~/\.nextflow\.log\.\d+/  
    new File('.').listFiles().findAll { it.name ==~ logPattern }.each { it.delete() }
}