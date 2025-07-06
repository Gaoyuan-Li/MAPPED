#!/usr/bin/env python3
import argparse
import pandas as pd
import csv

def main(metadata_file, sample_ids_file, output_file):
    # Read the batch sample IDs
    with open(sample_ids_file, 'r') as f:
        reader = csv.reader(f)
        next(reader)  # Skip header
        batch_sample_ids = set(row[0] for row in reader if row)
    
    print(f"Processing {len(batch_sample_ids)} samples for this batch")
    
    # Read the full metadata
    metadata_df = pd.read_csv(metadata_file, sep='\t', index_col=0)
    
    # Filter to only include samples in this batch
    batch_metadata = metadata_df[metadata_df.index.isin(batch_sample_ids)]
    
    print(f"Found {len(batch_metadata)} samples in metadata")
    
    # Create samplesheet in the expected format
    samplesheet_data = []
    
    for sample_id, row in batch_metadata.iterrows():
        # Parse the Run column which contains SRR IDs separated by semicolon
        run_ids = row['Run'].split(';') if pd.notna(row['Run']) else []
        
        for i, run_id in enumerate(run_ids):
            run_id = run_id.strip()
            if not run_id:
                continue
                
            # Determine layout
            layout = row.get('LibraryLayout', 'SINGLE').upper()
            
            # Create file paths
            if layout == 'PAIRED':
                fastq_1 = f"seqFiles/{run_id}_1.fastq.gz"
                fastq_2 = f"seqFiles/{run_id}_2.fastq.gz"
            else:
                fastq_1 = f"seqFiles/{run_id}.fastq.gz"
                fastq_2 = ""
            
            # Add to samplesheet
            # Use sample_id_run format for id to match the original pipeline
            samplesheet_data.append({
                'id': f"{sample_id}_{run_id}",
                'run_accession': run_id,
                'experiment_accession': sample_id,
                'fastq_1': fastq_1,
                'fastq_2': fastq_2,
                'layout': layout.lower()
            })
    
    # Write samplesheet
    with open(output_file, 'w', newline='') as f:
        if samplesheet_data:
            fieldnames = ['id', 'run_accession', 'experiment_accession', 'fastq_1', 'fastq_2', 'layout']
            writer = csv.DictWriter(f, fieldnames=fieldnames, quoting=csv.QUOTE_ALL)
            writer.writeheader()
            writer.writerows(samplesheet_data)
            print(f"Created samplesheet with {len(samplesheet_data)} entries")
        else:
            # Write empty samplesheet with header
            f.write('"id","run_accession","experiment_accession","fastq_1","fastq_2","layout"\n')
            print("Warning: No matching samples found for this batch")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Create batch-specific samplesheet from metadata")
    parser.add_argument("--metadata", required=True, help="Path to metadata TSV file")
    parser.add_argument("--sample-ids", required=True, help="Path to batch sample IDs CSV file")
    parser.add_argument("--output", required=True, help="Output samplesheet CSV file")
    args = parser.parse_args()
    
    main(args.metadata, args.sample_ids, args.output)