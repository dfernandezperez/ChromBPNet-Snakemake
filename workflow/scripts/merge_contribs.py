"""Merge chunks of chrombpnet contribution scores (from parallel_contribs.py) into the files written by
`chrombpnet contribs_bw`: <out-prefix>.interpreted_regions.bed, <out-prefix>.{profile,counts}_scores.h5
and <out-prefix>.{profile,counts}_scores.bw. Chunks are concatenated in the given order.
Must run inside the chrombpnet image.

Usage: python merge_contribs.py --parts chunk_000 chunk_001 ... --out-prefix out/prefix \
           --chrom-sizes genome.chrom.sizes [--heads profile counts]
"""
import argparse
import os
import sys
from concurrent.futures import ProcessPoolExecutor

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from contribs_utils import concat_beds, concat_h5, write_contribs_bigwig  # noqa: E402


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--parts", nargs="+", required=True, help="Output prefixes of the chunks, in order")
    parser.add_argument("--out-prefix", required=True)
    parser.add_argument("--chrom-sizes", required=True)
    parser.add_argument("--heads", nargs="+", default=["counts", "profile"], choices=["counts", "profile"])
    args = parser.parse_args()

    beds = [part + ".interpreted_regions.bed" for part in args.parts]
    out_bed = args.out_prefix + ".interpreted_regions.bed"
    for head in args.heads:
        print(f"Concatenating {len(args.parts)} {head} chunks", flush=True)
        concat_h5([f"{part}.{head}_scores.h5" for part in args.parts], beds, f"{args.out_prefix}.{head}_scores.h5")
    concat_beds(beds, out_bed)

    print("Writing bigwigs", flush=True)
    with ProcessPoolExecutor(len(args.heads)) as pool:
        jobs = [pool.submit(write_contribs_bigwig, f"{args.out_prefix}.{head}_scores.h5", out_bed,
                            args.chrom_sizes, f"{args.out_prefix}.{head}_scores.bw") for head in args.heads]
        for job in jobs:
            job.result()
    print("Done", flush=True)


if __name__ == "__main__":
    main()
