"""Random subset of the regions of a bed file, as `chrombpnet bias pipeline` selects the peaks it
interprets for its QC report (pandas `sample(n, random_state=seed)`; all regions if there are fewer).

Usage: python subsample_bed.py input.bed output.bed [--n 30000] [--seed 1234]
"""
import argparse

import pandas as pd

parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument("input")
parser.add_argument("output")
parser.add_argument("--n", type=int, default=30000)
parser.add_argument("--seed", type=int, default=1234)
args = parser.parse_args()

regions = pd.read_csv(args.input, sep="\t", header=None)
if regions.shape[0] > args.n:
    regions = regions.sample(args.n, random_state=args.seed)
regions.to_csv(args.output, sep="\t", header=False, index=False)
