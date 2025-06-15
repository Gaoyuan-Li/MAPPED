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
    container 'quay.io/biocontainers/gffread:0.9.12--0'

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
    cache true
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
// Process: merge count matrices by experiment ID
//
process MERGE_COUNTS {
    publishDir "${params.outdir}/expression_matrices", mode: 'copy'
    container 'felixlohmeier/pandas:1.3.3'

    input:
      path quant_dirs
      path passed_samples_file

    output:
      path 'tpm.tsv'
      path 'log_tpm.tsv'
      path 'counts.tsv'

    script:
    """
    python3 - << 'EOF'
import os
import pandas as pd
import numpy as np
from collections import defaultdict

# Read passed sample IDs
passed_sample_ids = set()
if os.path.exists('${passed_samples_file}'):
    with open('${passed_samples_file}', 'r') as f:
        for line in f:
            line = line.strip()
            if line:
                passed_sample_ids.add(line)

print(f"Found {len(passed_sample_ids)} passed sample IDs")

# Find all quant directories
quant_dirs = [d for d in os.listdir('.') if d.endswith('_quant')]
print(f"Found {len(quant_dirs)} quantification directories")

# Group by experiment ID and collect data
experiment_data = defaultdict(list)
gene_ids = None

for quant_dir in quant_dirs:
    sample_name = quant_dir.replace('_quant', '')
    
    # Extract experiment ID (everything before the first underscore)
    # This handles SRX_SRR, DRX_DRR, and ERX_ERR patterns
    experiment_id = sample_name.split('_')[0]
    
    # Get base sample name (remove _val_1/_val_2 suffixes)
    base_sample_name = sample_name.replace('_val_1', '').replace('_val_2', '')
    
    # Only process samples from passed experiments
    if base_sample_name not in passed_sample_ids:
        print(f"Skipping {sample_name} (base sample {base_sample_name} not in passed samples)")
        continue
    
    # Read quantification file
    quant_file = os.path.join(quant_dir, 'quant.sf')
    if not os.path.exists(quant_file):
        print(f"Warning: {quant_file} not found")
        continue
    
    df = pd.read_csv(quant_file, sep='\\t')
    
    # Store gene IDs from first file
    if gene_ids is None:
        gene_ids = df['Name'].tolist()
    
    # Store counts and TPM data for this experiment
    experiment_data[experiment_id].append({
        'sample': sample_name,
        'base_sample': base_sample_name,
        'counts': df['NumReads'].tolist(),
        'tpm': df['TPM'].tolist(),
        'length': df['Length'].tolist()
    })

print(f"Grouped into {len(experiment_data)} experiments")

# Merge data by experiment
final_counts = {}
final_tpm = {}

for experiment_id, runs in experiment_data.items():
    if len(runs) == 1:
        # Single run - use data directly, use experiment ID as column header
        final_counts[experiment_id] = runs[0]['counts']
        final_tpm[experiment_id] = runs[0]['tpm']
    else:
        # Multiple runs - sum counts and recalculate TPM
        print(f"Merging {len(runs)} runs for experiment {experiment_id}")
        
        # Sum counts across runs
        summed_counts = [0] * len(gene_ids)
        avg_lengths = [0] * len(gene_ids)
        
        for run in runs:
            for i in range(len(gene_ids)):
                summed_counts[i] += run['counts'][i]
                avg_lengths[i] = run['length'][i]  # Length should be same across runs
        
        # Calculate TPM from summed counts
        # TPM = (counts / length) * 1e6 / sum(counts / length)
        rpk = [summed_counts[i] / avg_lengths[i] if avg_lengths[i] > 0 else 0 for i in range(len(gene_ids))]
        scaling_factor = sum(rpk) / 1e6 if sum(rpk) > 0 else 1
        recalculated_tpm = [rpk[i] / scaling_factor if scaling_factor > 0 else 0 for i in range(len(gene_ids))]
        
        # Use experiment ID as column header
        final_counts[experiment_id] = summed_counts
        final_tpm[experiment_id] = recalculated_tpm

# Create output matrices
if gene_ids and final_counts:
    # Sort experiment IDs for consistent output
    sorted_experiments = sorted(final_counts.keys())
    
    # Write counts matrix
    with open('counts.tsv', 'w') as f:
        f.write('GeneID\\t' + '\\t'.join(sorted_experiments) + '\\n')
        for i, gene_id in enumerate(gene_ids):
            f.write(gene_id)
            for exp_id in sorted_experiments:
                f.write(f'\\t{final_counts[exp_id][i]}')
            f.write('\\n')
    
    # Write TPM matrix  
    with open('tpm.tsv', 'w') as f:
        f.write('GeneID\\t' + '\\t'.join(sorted_experiments) + '\\n')
        for i, gene_id in enumerate(gene_ids):
            f.write(gene_id)
            for exp_id in sorted_experiments:
                f.write(f'\\t{final_tpm[exp_id][i]:.6f}')
            f.write('\\n')
    
    # Write log TPM matrix (log2(TPM + 1))
    with open('log_tpm.tsv', 'w') as f:
        f.write('GeneID\\t' + '\\t'.join(sorted_experiments) + '\\n')
        for i, gene_id in enumerate(gene_ids):
            f.write(gene_id)
            for exp_id in sorted_experiments:
                log_tpm = np.log2(final_tpm[exp_id][i] + 1)
                f.write(f'\\t{log_tpm:.6f}')
            f.write('\\n')
    
    print(f"Generated matrices with {len(gene_ids)} genes and {len(sorted_experiments)} experiments")
else:
    print("No data to process - creating empty files")
    with open('counts.tsv', 'w') as f:
        f.write('GeneID\\n')
    with open('tpm.tsv', 'w') as f:
        f.write('GeneID\\n')
    with open('log_tpm.tsv', 'w') as f:
        f.write('GeneID\\n')

EOF
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
      path 'qc_summary.txt', emit: qc_summary_txt

    script:
    """
python3 - << 'EOF'
import json
from pathlib import Path

try:
    data = json.load(open('multiqc_data.json'))
    
    # Find the FastQC raw data with direct pass/fail/warn status
    fastqc_data = None
    
    if ('report_saved_raw_data' in data and 
        'multiqc_fastqc' in data['report_saved_raw_data']):
        fastqc_data = data['report_saved_raw_data']['multiqc_fastqc']
    
    if not fastqc_data:
        # Create empty files for both outputs
        open('passed_samples.txt', 'w').write('')
        # Create empty CSV with headers
        import csv
        with open('qc_summary.csv', 'w', newline='') as csvfile:
            writer = csv.writer(csvfile)
            writer.writerow(['sample', 'per_base_sequence_quality', 'per_sequence_quality_scores', 'per_base_n_content', 'overall_status'])
        # Create summary text file
        with open('qc_summary.txt', 'w') as f:
            f.write("QC Summary Report\\n")
            f.write("=================\\n")
            f.write("Total samples processed: 0\\n")
            f.write("Samples passed: 0\\n")
            f.write("Samples failed: 0\\n")
            f.write("\\nError: No FastQC raw data found in MultiQC JSON\\n")
        exit(0)
    
    # Target metrics we want to extract
    target_metrics = ['per_base_sequence_quality', 'per_sequence_quality_scores', 'per_base_n_content']
    
    # Process each sample and group by experiment ID (SRX, DRX, or ERX)
    experiment_status = {}  # Track status for each experiment
    experiment_samples = {}  # Track which samples belong to each experiment
    qc_results = []
    
    for sample_name, sample_data in fastqc_data.items():
        # Extract experiment ID (everything before the first underscore)
        # This handles SRX_SRR, DRX_DRR, and ERX_ERR patterns
        experiment_id = sample_name.split('_')[0]
        
        # Remove _val_1/_val_2 suffix to get base sample name
        base_sample_name = sample_name.replace('_val_1', '').replace('_val_2', '')
        
        ok = True
        failed_metrics = []
        sample_qc = {'sample': sample_name}
        
        # Extract the three critical metrics directly
        for metric in target_metrics:
            status = sample_data.get(metric, 'unknown')
            sample_qc[metric] = status
            
            # If status is not 'pass', mark sample as failed
            if status != 'pass':
                ok = False
                failed_metrics.append(metric)
        
        # Add overall status
        sample_qc['overall_status'] = 'PASS' if ok else 'FAIL'
        qc_results.append(sample_qc)
        
        # Track experiment status - if any sample for an experiment fails, the whole experiment fails
        if experiment_id not in experiment_status:
            experiment_status[experiment_id] = True
            experiment_samples[experiment_id] = set()
        
        experiment_samples[experiment_id].add(base_sample_name)
        
        if not ok:
            experiment_status[experiment_id] = False
    
    # Generate passed sample list (full sample IDs, not just experiment IDs)
    passed_samples = set()
    failed_samples = []
    
    for experiment_id, status in experiment_status.items():
        if status:
            # Add all base sample names for this experiment
            passed_samples.update(experiment_samples[experiment_id])
        else:
            # Find which samples failed for this experiment
            for sample_name, sample_data in fastqc_data.items():
                if sample_name.split('_')[0] == experiment_id:
                    base_sample_name = sample_name.replace('_val_1', '').replace('_val_2', '')
                    failed_metrics = []
                    for metric in target_metrics:
                        status = sample_data.get(metric, 'unknown')
                        if status != 'pass':
                            failed_metrics.append(metric)
                    if failed_metrics:  # Only include samples that actually failed
                        failed_samples.append((sample_name, failed_metrics))
    
    # Write passed sample IDs to file (full sample IDs like SRX_SRR or DRX_DRR)
    with open('passed_samples.txt', 'w') as f:
        for sample_id in sorted(passed_samples):
            f.write(sample_id + '\\n')
    
    # Write QC summary CSV
    import csv
    with open('qc_summary.csv', 'w', newline='') as csvfile:
        fieldnames = ['sample', 'per_base_sequence_quality', 'per_sequence_quality_scores', 'per_base_n_content', 'overall_status']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        
        if qc_results:
            for row in qc_results:
                writer.writerow(row)
    
    # Write QC summary text file
    with open('qc_summary.txt', 'w') as f:
        f.write("QC Summary Report\\n")
        f.write("=================\\n")
        f.write(f"Total individual samples processed: {len(fastqc_data)}\\n")
        f.write(f"Total experiments processed: {len(experiment_status)}\\n")
        f.write(f"Experiments passed: {sum(1 for status in experiment_status.values() if status)}\\n")
        f.write(f"Experiments failed: {sum(1 for status in experiment_status.values() if not status)}\\n")
        f.write(f"Individual samples failed: {len(failed_samples)}\\n")
        f.write(f"Unique sample IDs passed: {len(passed_samples)}\\n")
        
        if failed_samples:
            f.write("\\nFailed individual samples:\\n")
            for sample, metrics in failed_samples:
                f.write(f"  {sample}: {', '.join(metrics)}\\n")
    
except Exception as e:
    import traceback
    # Create empty files to avoid pipeline failure
    open('passed_samples.txt', 'w').write('')
    # Create empty CSV with headers
    import csv
    with open('qc_summary.csv', 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(['sample', 'per_base_sequence_quality', 'per_sequence_quality_scores', 'per_base_n_content', 'overall_status'])
    # Create error summary text file
    with open('qc_summary.txt', 'w') as f:
        f.write("QC Summary Report\\n")
        f.write("=================\\n")
        f.write("Total samples processed: 0\\n")
        f.write("Samples passed: 0\\n")
        f.write("Samples failed: 0\\n")
        f.write(f"\\nError processing MultiQC data: {e}\\n")
        f.write(f"Traceback: {traceback.format_exc()}\\n")
EOF
    """
}

//
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
    # Copy original samplesheet and remove empty lines
    grep -v '^[[:space:]]*\$' ${samplesheet} > tmp.csv
    
    # Check if passed samples file is empty
    if [ ! -s ${passedlist} ]; then
        echo "WARNING: No samples passed QC filters!"
        # Create samplesheet with only header
        head -n1 tmp.csv > samplesheet.csv
    else
        # Copy header
        head -n1 tmp.csv > samplesheet.csv
        
        # Filter rows based on passed sample IDs
        # The id column (4th column) contains the full sample ID (SRX_SRR, DRX_DRR, or ERX_ERR format)
        while IFS= read -r sample_id; do
            if [ -n "\$sample_id" ]; then
                # Look for lines where the id column (4th field) matches the sample ID exactly
                # Use awk to properly parse CSV and check the id column
                awk -F',' -v sample="\$sample_id" '
                    NR > 1 && \$4 == "\\"" sample "\\"" { print }
                ' tmp.csv >> samplesheet.csv || true
            fi
        done < ${passedlist}
        
        # Remove duplicates while preserving order
        awk '!seen[\$0]++' samplesheet.csv > samplesheet_dedup.csv
        mv samplesheet_dedup.csv samplesheet.csv
    fi
    
    # Report results
    original_count=\$(tail -n +2 tmp.csv | grep -c '^')
    filtered_count=\$(tail -n +2 samplesheet.csv | grep -c '^')
    echo "Original samples: \$original_count"
    echo "Filtered samples: \$filtered_count"
    echo "Samples removed: \$((\$original_count - \$filtered_count))"
    """
}

//
// Process: filter out samples with >50% zero values from expression matrices and samplesheet
//
process FILTER_LOW_EXPRESSION_SAMPLES {
    tag 'filter_low_expression'
    container 'felixlohmeier/pandas:1.3.3'
    publishDir "${params.outdir}/expression_matrices", mode: 'copy', overwrite: true, pattern: "*.tsv"
    publishDir "${params.outdir}/samplesheet", mode: 'copy', overwrite: true, pattern: "*.csv"

    input:
      path tpm_matrix
      path log_tpm_matrix  
      path counts_matrix
      path samplesheet

    output:
      path 'tpm.tsv', emit: tpm
      path 'log_tpm.tsv', emit: log_tpm
      path 'counts.tsv', emit: counts
      path 'samplesheet.csv', emit: samplesheet

    script:
    """
    python3 - << 'EOF'
import pandas as pd
import numpy as np

print("Reading expression matrices...")
# Read the three expression matrices
tpm_df = pd.read_csv('${tpm_matrix}', sep='\\t', index_col=0)
log_tpm_df = pd.read_csv('${log_tpm_matrix}', sep='\\t', index_col=0)
counts_df = pd.read_csv('${counts_matrix}', sep='\\t', index_col=0)

print(f"Original matrices shape: {tpm_df.shape}")
print(f"Sample columns: {list(tpm_df.columns)}")

# Identify samples with >50% zero values
samples_to_remove = []
total_genes = len(tpm_df)

for sample in tpm_df.columns:
    # Count zeros in the counts matrix (most appropriate for this analysis)
    zero_count = (counts_df[sample] == 0).sum()
    zero_percentage = zero_count / total_genes
    
    print(f"Sample {sample}: {zero_count}/{total_genes} zeros ({zero_percentage:.2%})")
    
    if zero_percentage > 0.5:
        samples_to_remove.append(sample)
        print(f"  -> Will be removed (>50% zeros)")

print(f"\\nSamples to remove ({len(samples_to_remove)}): {samples_to_remove}")

# Filter matrices by removing low-expression samples
if samples_to_remove:
    samples_to_keep = [col for col in tpm_df.columns if col not in samples_to_remove]
    
    tpm_filtered = tpm_df[samples_to_keep]
    log_tpm_filtered = log_tpm_df[samples_to_keep]
    counts_filtered = counts_df[samples_to_keep]
    
    print(f"Filtered matrices shape: {tpm_filtered.shape}")
else:
    print("No samples need to be removed")
    tpm_filtered = tpm_df
    log_tpm_filtered = log_tpm_df  
    counts_filtered = counts_df

# Save filtered matrices (overwrite originals)
tpm_filtered.to_csv('tpm.tsv', sep='\\t')
log_tpm_filtered.to_csv('log_tpm.tsv', sep='\\t')
counts_filtered.to_csv('counts.tsv', sep='\\t')

print("Reading and filtering samplesheet...")
# Read and filter samplesheet
samplesheet_df = pd.read_csv('${samplesheet}')
print(f"Original samplesheet shape: {samplesheet_df.shape}")

if samples_to_remove:
    # Filter samplesheet by removing rows where 'id' column matches samples to remove
    # The sample ID in the matrix corresponds to the 'id' column in samplesheet
    samplesheet_filtered = samplesheet_df[~samplesheet_df['id'].isin(samples_to_remove)]
    print(f"Filtered samplesheet shape: {samplesheet_filtered.shape}")
    print(f"Removed {len(samplesheet_df) - len(samplesheet_filtered)} rows from samplesheet")
else:
    samplesheet_filtered = samplesheet_df

# Save filtered samplesheet (overwrite original)
samplesheet_filtered.to_csv('samplesheet.csv', index=False)

print("\\nFiltering completed successfully!")
print(f"Final matrices have {len(tpm_filtered.columns)} samples and {len(tpm_filtered)} genes")

EOF
    """
}

//
// Main workflow
//
workflow {
    // load samples from the original download samplesheet
    samples_ch = Channel
        .fromPath( file("${params.outdir}/samplesheet/samplesheet_download.csv") )
        .ifEmpty { error "Sample sheet not found at: ${params.outdir}/samplesheet/samplesheet_download.csv" }
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
            row.id && row.fastq_1 && row.fastq_2 && row.run_accession
        }
        .ifEmpty {
            error """
    No valid samples found after filtering.
    Make sure your CSV has columns exactly named:
    id, fastq_1, fastq_2, run_accession.
    """
        }
        // build the tuple for each sample using the id column which has SRX_SRR, DRX_DRR, or ERX_ERR format
        .map { row ->
            tuple(
                row.id,
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
    qc_summary_txt_ch = parse_qc_result.qc_summary_txt

    // Copy QC summary files to multiqc folder
    qc_summary_ch.subscribe { qc_file ->
        def target_dir = file("${params.outdir}/multiqc")
        target_dir.mkdirs()
        qc_file.copyTo(target_dir.resolve("qc_summary.csv"))
    }
    
    qc_summary_txt_ch.subscribe { qc_file ->
        def target_dir = file("${params.outdir}/multiqc")
        target_dir.mkdirs()
        qc_file.copyTo(target_dir.resolve("qc_summary.txt"))
    }

    // Filter original sample sheet based on passed samples
    filtered_samplesheet_ch = FILTER_SAMPLESHEET( file("${params.outdir}/samplesheet/samplesheet_download.csv"), passed_ch )

    // Create a channel of passed sample IDs
    passed_samples_ch = passed_ch
        .splitText()
        .map { it.trim() }
        .filter { it }
        .map { sample_id -> tuple(sample_id, sample_id) }

    // Transform trimmed channel to include sample ID as key
    trimmed_with_sample_ch = trimmed_ch
        .map { sample_tuple ->
            def sample_id = sample_tuple[0]
            // Remove _val_1/_val_2 suffixes to get base sample name
            def base_sample_id = sample_id.replace('_val_1', '').replace('_val_2', '')
            tuple(base_sample_id, sample_tuple)
        }

    // Join to filter only QC-passed samples
    filtered_trimmed_ch = trimmed_with_sample_ch
        .join(passed_samples_ch)
        .map { base_sample_id, sample_tuple, passed_sample -> sample_tuple }

    // Run Salmon quantification only on QC-passed samples
    quant_ch = SALMON_QUANT(filtered_trimmed_ch, salmon_index_ch)

    // merge count matrices - wait for both quantification and samplesheet filtering
    count_matrix_ch = MERGE_COUNTS( 
        quant_ch.collect(), 
        passed_ch
    )

    // Filter out samples with >50% zero values from expression matrices and samplesheet
    filtered_results = FILTER_LOW_EXPRESSION_SAMPLES(
        count_matrix_ch[0],  // tpm.tsv
        count_matrix_ch[1],  // log_tpm.tsv  
        count_matrix_ch[2],  // counts.tsv
        filtered_samplesheet_ch
    )
}

// Add an onComplete event handler to always delete rotated Nextflow log files
workflow.onComplete {
    def logPattern = ~/\.nextflow\.log\.\d+/  
    new File('.').listFiles().findAll { it.name ==~ logPattern }.each { it.delete() }
}