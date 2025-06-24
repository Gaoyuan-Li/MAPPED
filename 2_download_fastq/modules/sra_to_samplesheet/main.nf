process SRA_TO_SAMPLESHEET {
    tag "$meta.id"
    errorStrategy 'ignore'

    executor 'local'
    memory 100.MB

    input:
    val meta
    val pipeline
    val strandedness
    val mapping_fields

    output:
    tuple val(meta), path("*samplesheet.csv"), emit: samplesheet, optional: true

    script:
    // Build only the samplesheet with local fastq paths
    def mclone = meta.clone()
    ['fastq_1','fastq_2','md5_1','md5_2','single_end'].each { mclone.remove(it) }
    def sampleId = meta.id.split('_')[0..-2].join('_')
    def fastqDir = "${params.workdir}/seqFiles/fastq"
    def baseMap = [ sample: sampleId ]
    if (meta.single_end.toString().toBoolean()) {
        baseMap.fastq_1 = "${fastqDir}/${meta.id}.fastq.gz"
        baseMap.fastq_2 = ""
    } else {
        baseMap.fastq_1 = "${fastqDir}/${meta.id}_1.fastq.gz"
        baseMap.fastq_2 = "${fastqDir}/${meta.id}_2.fastq.gz"
    }
    def pipeline_map = baseMap + mclone
    def header = pipeline_map.keySet().collect{ "\"${it}\"" }.join(',')
    def values = pipeline_map.values().collect{ "\"${it}\"" }.join(',')
    return """#!/usr/bin/env bash
echo '${header}' > ${meta.id}.samplesheet.csv
echo '${values}' >> ${meta.id}.samplesheet.csv
"""
}
