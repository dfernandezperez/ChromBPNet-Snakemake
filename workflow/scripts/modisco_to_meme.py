"""Export the motifs (patterns) of a TF-MoDISco results h5 to a MEME motif file, for tools such as
FIMO, Tomtom, HOMER (after conversion) or TOBIAS.

Each pattern's position frequency matrix (PFM) is trimmed like the logos of the TF-MoDISco report:
positions where sum(|CWM|) >= trim_threshold * max, plus `flank` bp on each side. Motif ids are the
TF-MoDISco names (e.g. pos_patterns.pattern_0, as in the report and the Fi-NeMo hits); the alternative
name adds the given prefix (e.g. the sample) so files of several samples can be concatenated.
Negative patterns (sequences with negative contributions) are included unless --positive-only.

Usage: python modisco_to_meme.py -i modisco.h5 -o motifs.meme [--name-prefix sample] [--positive-only]
"""
import argparse

import h5py
import numpy as np


def trimmed_pfm(pattern, trim_threshold, flank):
    pfm = np.array(pattern["sequence"][:], dtype=float)
    score = np.abs(np.array(pattern["contrib_scores"][:])).sum(axis=1)
    passing = np.where(score >= trim_threshold * score.max())[0]
    start, end = max(passing.min() - flank, 0), min(passing.max() + flank + 1, len(score))
    pfm = pfm[start:end]
    totals = pfm.sum(axis=1, keepdims=True)
    return np.where(totals > 0, pfm / np.where(totals > 0, totals, 1), 0.25)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-i", "--modisco-h5", required=True)
    parser.add_argument("-o", "--output", required=True)
    parser.add_argument("--name-prefix", default="", help="Added to the alternative motif names")
    parser.add_argument("--trim-threshold", type=float, default=0.3)
    parser.add_argument("--flank", type=int, default=4)
    parser.add_argument("--positive-only", action="store_true")
    args = parser.parse_args()

    metaclusters = ["pos_patterns"] if args.positive_only else ["pos_patterns", "neg_patterns"]
    prefix = args.name_prefix + "_" if args.name_prefix else ""
    n_motifs = 0
    with h5py.File(args.modisco_h5, "r") as results, open(args.output, "w") as out:
        out.write("MEME version 4\n\nALPHABET= ACGT\n\nstrands: + -\n\n")
        out.write("Background letter frequencies\nA 0.25 C 0.25 G 0.25 T 0.25\n\n")
        for metacluster in metaclusters:
            if metacluster not in results:
                continue
            patterns = sorted(results[metacluster].items(), key=lambda item: int(item[0].split("_")[-1]))
            for name, pattern in patterns:
                pfm = trimmed_pfm(pattern, args.trim_threshold, args.flank)
                n_seqlets = int(pattern["seqlets"]["n_seqlets"][()][0])
                short = metacluster.split("_")[0]
                out.write(f"MOTIF {metacluster}.{name} {prefix}{short}_{name}\n")
                out.write(f"letter-probability matrix: alength= 4 w= {len(pfm)} nsites= {n_seqlets} E= 0\n")
                for row in pfm:
                    out.write(" ".join(f"{p:.6f}" for p in row) + "\n")
                out.write("\n")
                n_motifs += 1
    print(f"Wrote {n_motifs} motifs to {args.output}")


if __name__ == "__main__":
    main()
