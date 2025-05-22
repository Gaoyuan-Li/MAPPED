import argparse

import pandas as pd


def main(infile, outfile, layout):

    # Read in metadata with one SRR as row
    DF_SRR = pd.read_csv(infile, header=0)

    # Merge internal tables
    DF_SRR = DF_SRR[~DF_SRR.Run.isin(["", "Run"])]

    # Merge rows so it has one SRX per row
    agg_rules = {col: ";".join if col == "Run" else "last" for col in DF_SRR.columns}
    DF_SRX = DF_SRR.groupby("Experiment").agg(agg_rules)

    DF_SRX.set_index("Experiment", inplace=True)

    # Add R1 and R2 columns
    DF_SRX["R1"] = None
    DF_SRX["R2"] = None

    # Sort rows by ReleaseDate descending
    if 'ReleaseDate' in DF_SRX.columns:
        DF_SRX['ReleaseDate'] = pd.to_datetime(DF_SRX['ReleaseDate'], format='%Y-%m-%d', errors='coerce')
        DF_SRX = DF_SRX.sort_values('ReleaseDate', ascending=False)

    # Filter by LibraryLayout parameter
    if 'LibraryLayout' in DF_SRX.columns:
        lib_vals = DF_SRX['LibraryLayout'].str.upper()
        if layout.lower() == 'paired':
            DF_SRX = DF_SRX[ lib_vals == 'PAIRED' ]
        elif layout.lower() == 'single':
            DF_SRX = DF_SRX[ lib_vals == 'SINGLE' ]
        else:
            DF_SRX = DF_SRX[ lib_vals.isin(['PAIRED','SINGLE']) ]

    # Save to file
    DF_SRX.to_csv(outfile, sep="\t")

    # Save only sample IDs (first column) to CSV without header
    pd.Series(DF_SRX.index).to_csv('sample_id.csv', index=False, header=True)


if __name__ == "__main__":
    # Argument parsing
    p = argparse.ArgumentParser(description="Clean raw SRA metadata")
    p.add_argument("-i", "--input", help="Input filename")
    p.add_argument("-o", "--output", help="Output filename")
    p.add_argument("-l", "--layout", choices=['paired','single','both'], default='both', help="Library layout to include: paired, single, or both")
    args = p.parse_args()

    main(args.input, args.output, args.layout)
