import argparse

import pandas as pd


def main(infile, outfile, layout, strain=None):

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

    # Optional: Filter by strain using ScientificName column
    # - Split ScientificName by spaces and keep rows where any token equals
    #   or contains the provided strain string (case-insensitive)
    if strain:
        col = 'ScientificName'
        if col in DF_SRX.columns:
            s = DF_SRX[col].fillna("").astype(str)
            query = str(strain).strip()
            if query:
                qlower = query.lower()
                def match_tokens(name: str) -> bool:
                    tokens = name.split()
                    for t in tokens:
                        tlower = t.lower()
                        if tlower == qlower or (qlower in tlower):
                            return True
                    return False
                mask = s.apply(match_tokens)
                DF_SRX = DF_SRX[mask]

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
    p.add_argument("--strain", default=None, help="Optional strain filter applied to ScientificName; matches if any space-delimited token equals or contains the provided string (case-insensitive)")
    args = p.parse_args()

    main(args.input, args.output, args.layout, args.strain)
