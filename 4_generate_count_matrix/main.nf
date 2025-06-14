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
    cache true
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
    fastqc --threads 4 -o . ${fq1} ${fq2}
    """
}

//
// Process: TrimGalore
//
process TRIMGALORE {
    cache true
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
    trim_galore --cores 4 --paired --basename ${sample} --output_dir . ${fq1} ${fq2}
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
      --validateMappings \
      -p 4
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
    echo "Found QC files:"
    ls -la
    
    echo "Running MultiQC..."
    multiqc . -o . --data-format json -v
    
    echo "MultiQC output files:"
    ls -la multiqc*
    
    # Check for JSON in both possible locations
    if [ -f multiqc_data.json ]; then
        echo "MultiQC JSON found in current directory"
        echo "MultiQC JSON size: \$(wc -c < multiqc_data.json) bytes"
        echo "First few lines of JSON:"
        head -n 5 multiqc_data.json
    elif [ -f multiqc_data/multiqc_data.json ]; then
        echo "MultiQC JSON found in multiqc_data directory"
        echo "MultiQC JSON size: \$(wc -c < multiqc_data/multiqc_data.json) bytes"
        echo "First few lines of JSON:"
        head -n 5 multiqc_data/multiqc_data.json
        # Move it to the expected location for the next process
        cp multiqc_data/multiqc_data.json ./multiqc_data.json
        echo "Copied multiqc_data.json to current directory"
    else
        echo "ERROR: multiqc_data.json was not found in current directory or multiqc_data/ subdirectory!"
        echo "Available files:"
        find . -name "*.json" -type f
        exit 1
    fi
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
      path 'qc_summary.csv', emit: qc_summary

    script:
    """
python3 - << 'EOF'
import json
from pathlib import Path

try:
    data = json.load(open('multiqc_data.json'))
    print("Successfully loaded MultiQC JSON data")
    
    # Find the FastQC raw data with direct pass/fail/warn status
    fastqc_data = None
    
    if ('report_saved_raw_data' in data and 
        'multiqc_fastqc' in data['report_saved_raw_data']):
        fastqc_data = data['report_saved_raw_data']['multiqc_fastqc']
        print("Found FastQC raw data in report_saved_raw_data/multiqc_fastqc")
    
    if not fastqc_data:
        print("Warning: Could not find FastQC raw data in MultiQC JSON")
        print("Available keys:", list(data.keys()))
        if 'report_saved_raw_data' in data:
            print("Available raw data keys:", list(data['report_saved_raw_data'].keys()))
        # Create empty files for both outputs
        open('passed_samples.txt', 'w').write('')
        # Create empty CSV with headers
        import csv
        with open('qc_summary.csv', 'w', newline='') as csvfile:
            writer = csv.writer(csvfile)
            writer.writerow(['sample', 'per_base_sequence_quality', 'per_sequence_quality_scores', 'per_base_n_content', 'overall_status'])
        print("Created empty output files (no FastQC raw data found)")
        exit(0)
    
    print(f"Found FastQC data for {len(fastqc_data)} samples")
    
    # Target metrics we want to extract
    target_metrics = ['per_base_sequence_quality', 'per_sequence_quality_scores', 'per_base_n_content']
    
    # Process each sample
    passed = []
    failed_samples = []
    qc_results = []
    
    for sample_name, sample_data in fastqc_data.items():
        print(f"Processing sample: {sample_name}")
        ok = True
        failed_metrics = []
        sample_qc = {'sample': sample_name}
        
        # Extract the three critical metrics directly
        for metric in target_metrics:
            status = sample_data.get(metric, 'unknown')
            sample_qc[metric] = status
            
            print(f"  {metric}: {status}")
            
            # If status is not 'pass', mark sample as failed
            if status != 'pass':
                ok = False
                failed_metrics.append(metric)
        
        # Add overall status
        sample_qc['overall_status'] = 'PASS' if ok else 'FAIL'
        qc_results.append(sample_qc)
        
        if ok:
            passed.append(sample_name)
            print(f"  -> PASSED")
        else:
            failed_samples.append((sample_name, failed_metrics))
            print(f"  -> FAILED ({', '.join(failed_metrics)})")
    
    print(f"Total samples processed: {len(fastqc_data)}")
    print(f"Samples passed: {len(passed)}")
    print(f"Samples failed: {len(failed_samples)}")
    
    if failed_samples:
        print("Failed samples:")
        for sample, metrics in failed_samples:
            print(f"  {sample}: {', '.join(metrics)}")
    
    # Write passed samples to file
    with open('passed_samples.txt', 'w') as f:
        f.write('\\n'.join(passed))
        if passed:
            f.write('\\n')  # Add final newline
    
    print(f"Written {len(passed)} passed samples to passed_samples.txt")
    
    # Write QC summary CSV
    import csv
    with open('qc_summary.csv', 'w', newline='') as csvfile:
        fieldnames = ['sample', 'per_base_sequence_quality', 'per_sequence_quality_scores', 'per_base_n_content', 'overall_status']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        
        if qc_results:
            for row in qc_results:
                writer.writerow(row)
            print(f"Written QC summary CSV with {len(qc_results)} samples")
        else:
            print("Created empty QC summary CSV (no samples found)")
    
    print("QC summary CSV created successfully")
    
except Exception as e:
    print(f"Error processing MultiQC data: {e}")
    import traceback
    traceback.print_exc()
    # Create empty files to avoid pipeline failure
    open('passed_samples.txt', 'w').write('')
    # Create empty CSV with headers
    import csv
    with open('qc_summary.csv', 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(['sample', 'per_base_sequence_quality', 'per_sequence_quality_scores', 'per_base_n_content', 'overall_status'])
    print("Created empty output files due to error")
EOF
    """
}

// Process: filter sample sheet CSV based on passed samples
process FILTER_SAMPLESHEET {
    tag 'filter_samplesheet'
    container 'ubuntu:22.04'
    publishDir "${params.outdir}/samplesheet", mode: 'copy', overwrite: true

    input:
      path samplesheet
      path passedlist

    output:
      path 'samplesheet.csv'

    script:
    """
    # Copy original samplesheet
    cp ${samplesheet} tmp.csv
    
    # Check if passed samples file is empty
    if [ ! -s ${passedlist} ]; then
        echo "WARNING: No samples passed QC filters!"
        # Create samplesheet with only header
        head -n1 tmp.csv > samplesheet.csv
    else
        # Copy header
        head -n1 tmp.csv > samplesheet.csv
        
        # Filter rows based on passed samples
        # Use a more robust approach that handles sample names that may be part of larger strings
        while IFS= read -r sample_id; do
            if [ -n "\$sample_id" ]; then
                # Look for lines where the sample_id appears at the beginning of a field
                grep "^\$sample_id\\|,\$sample_id" tmp.csv >> samplesheet.csv || true
            fi
        done < ${passedlist}
        
        # Remove duplicates while preserving order
        awk '!seen[\$0]++' samplesheet.csv > samplesheet_dedup.csv
        mv samplesheet_dedup.csv samplesheet.csv
    fi
    
    # Report results
    original_count=\$(tail -n +2 tmp.csv | wc -l)
    filtered_count=\$(tail -n +2 samplesheet.csv | wc -l)
    echo "Original samples: \$original_count"
    echo "Filtered samples: \$filtered_count"
    echo "Samples removed: \$((\$original_count - \$filtered_count))"
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
    parse_qc_result = PARSE_QC( multiqc_json_ch )
    passed_ch = parse_qc_result.passlist
    qc_summary_ch = parse_qc_result.qc_summary

    // Copy QC summary to multiqc folder
    qc_summary_ch.subscribe { qc_file ->
        def target_dir = file("${params.outdir}/multiqc")
        target_dir.mkdirs()
        qc_file.copyTo(target_dir.resolve("qc_summary.csv"))
    }

    // Filter original sample sheet based on passed samples
    FILTER_SAMPLESHEET( file("${params.outdir}/samplesheet/samplesheet.csv"), passed_ch )

    // Filter trimmed reads and perform Salmon quantification only on passed samples
    // Read the passed samples list and convert to channel
    passed_samples_ch = passed_ch
        .splitText()
        .map { it.trim() }
        .filter { it }  // Remove empty lines
    
    // Create a set of passed sample IDs for filtering
    passed_set_ch = passed_samples_ch.collect().map { it.toSet() }
    
    // Filter trimmed channel to only include samples that passed QC
    filtered_trimmed_ch = trimmed_ch
        .combine(passed_set_ch)
        .filter { sample_tuple, passed_set ->
            def sample_id = sample_tuple[0]
            return passed_set.contains(sample_id)
        }
        .map { sample_tuple, passed_set -> sample_tuple }
    
    quant_ch = SALMON_QUANT(filtered_trimmed_ch, salmon_index_ch)

    // merge count matrices
    count_matrix_ch = MERGE_COUNTS( quant_ch.collect() )
}

// Add an onComplete event handler to always delete rotated Nextflow log files
workflow.onComplete {
    def logPattern = ~/\.nextflow\.log\.\d+/  
    new File('.').listFiles().findAll { it.name ==~ logPattern }.each { it.delete() }
}