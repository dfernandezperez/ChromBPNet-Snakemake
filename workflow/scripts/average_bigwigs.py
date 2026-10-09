"""Average bigwig files (e.g. the same track from different ChromBPNet folds) base by base.

Positions covered by at least one input are written with the mean over the inputs that
cover them; uncovered positions are left empty. Chromosomes and their order are taken from
the first bigwig.

Usage: python average_bigwigs.py -o mean.bw fold_0.bw fold_1.bw ...
"""
import argparse

import numpy as np
import pyBigWig


def contiguous_blocks(mask):
    """Return (start, end) pairs of the runs of True values in a boolean array."""
    padded = np.concatenate(([False], mask, [False]))
    edges = np.flatnonzero(np.diff(padded.astype(np.int8)))
    return edges[0::2], edges[1::2]


def average_bigwigs(input_bws, output_bw):
    bws = [pyBigWig.open(path) for path in input_bws]
    chroms = bws[0].chroms() # dict keeps the chromosome order of the file

    out = pyBigWig.open(output_bw, "w")
    out.addHeader(list(chroms.items()), maxZooms=10)

    for chrom, length in chroms.items():
        total = np.zeros(length, dtype=np.float32)
        count = np.zeros(length, dtype=np.uint16)
        for bw in bws:
            if chrom not in bw.chroms() or bw.intervals(chrom) is None:
                continue
            values = np.asarray(bw.values(chrom, 0, length, numpy=True), dtype=np.float32)
            covered = ~np.isnan(values)
            total[covered] += values[covered]
            count[covered] += 1
            del values, covered

        covered = count > 0
        if not covered.any():
            continue
        mean = np.zeros(length, dtype=np.float32)
        mean[covered] = total[covered] / count[covered]
        del total, count

        starts, ends = contiguous_blocks(covered)
        for start, end in zip(starts, ends):
            out.addEntries(chrom, int(start), values=mean[start:end].astype(np.float64), span=1, step=1)

    out.close()
    for bw in bws:
        bw.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-o", "--output", required=True, help="Output bigwig")
    parser.add_argument("inputs", nargs="+", help="Input bigwigs")
    args = parser.parse_args()
    average_bigwigs(args.inputs, args.output)
