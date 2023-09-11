import argparse
from snpio import GenotypeData


def main():
    parser = argparse.ArgumentParser(
        description="Generate missingness report from Genotype data."
    )
    parser.add_argument("-p", "--phylip", required=True, help="Input phylip file.")
    parser.add_argument("-m", "--popmap", required=True, help="Input popmap file.")
    parser.add_argument(
        "-o",
        "--output_prefix",
        default="missingness_report",
        help="Prefix for the output file. Default is 'missingness_report'.",
    )

    args = parser.parse_args()

    # Read the alignment, popmap, and tree files
    gd = GenotypeData(
        filename=args.phylip,
        popmapfile=args.popmap,
        force_popmap=True,
        filetype="auto",
    )

    # Make missingness report plots.
    gd.missingness_reports(prefix=args.output_prefix)


if __name__ == "__main__":
    main()
