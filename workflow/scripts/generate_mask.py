import argparse
import copy
from snpio import GenotypeData
from pgsui.data_processing.transformers import SimGenotypeDataTransformer

def parse_args():
    parser = argparse.ArgumentParser(
        description="Generate a random mask for a sequence alignment"
    )

    parser.add_argument("-i", "--input", type=str, required=True, help="Input Phylip file.")
    parser.add_argument("-p", "--popmap", type=str, required=True, help="Popmap file.")
    parser.add_argument("-t", "--tree", type=str, required=True, help="Tree file.")
    parser.add_argument("-m", "--missing", type=float, required=True, help="Missing data proportion.")
    parser.add_argument("-s", "--strategy", type=str, required=True, help="Missing data strategy.")
    parser.add_argument("-o", "--output", type=str, required=True, help="Output file prefix.")

    return parser.parse_args()

def generate_mask(args):
    genotype_data = GenotypeData(
        filename=args.input,
        popmapfile=args.popmap,
        guidetree=args.tree,
        prefix=args.output,
        force_popmap=True,
        plot_format="png",
    )

    transformer = SimGenotypeDataTransformer(
        genotype_data=genotype_data,
        prop_missing=args.missing,
        strategy=args.strategy,
    )

    transformer.fit(genotype_data.genotypes_012(fmt="numpy"))
    simulated_data = copy.deepcopy(genotype_data)
    simulated_data.genotypes_012 = transformer.transform(genotype_data.genotypes_012(fmt="numpy"))
    transformer.write_mask(filename_prefix=args.output)

if __name__ == "__main__":
    args = parse_args()
    generate_mask(args)
