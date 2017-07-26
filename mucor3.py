import argparse
import os
import sys

# Local imports
import annfield

def consume_vcf(fn):

    return pd.DataFrame()

def greeting():
    print("Mucor 3")

def valediction():
    print("Have a nice day.")

def main():
    greeting()

    parser = argparse.ArgumentParser()
    parser.add_argument('-a', '--annotation', action='append', help="Annotation field to include in output (default: all)")
    parser.add_argument('-i', '--input', action='append', help="Input VCF file")
    parser.add_argument('-o', '--output', help="Output file base name (CSV and XLSX)")
    args = parser.parse_args()

    # Prevalidate so as not to interrupt long-running procedure
    for vcf_file in args.input:
        # Validate file exists

    master_df = pd.DataFrame()

    for vcf_file in args.input:
        vcf_df = consume_vcf(fn)
        # TODO merge to main dataframe

    # TODO subset on columns included in --annotation
    if args.annotation:
        # TODO print diagnostic pre subset
        master_df = master_df[ args.annotation ]
        # TODO print diagnostic post subset

    # TODO pivot and/or aggregate
    # Option 1, aggregate on location (+- X>Y change)
    # Option 2, aggregate on gene symbol from annotations

    # TODO SHould we include any summary metrics (row wise or group wise) ?

    master_df.to_csv( args.output + ".csv" )

    # TODO output summary JSON

    valediction()

if __name__ == "__main__":
    main()
