"""Keep the top N peaks by significance (narrowPeak column 7) located on the canonical chromosomes.

Peaks on other chromosomes (chrM, unplaced contigs...) are removed first: chrombpnet only trains on the
chromosomes of the fold files, and e.g. chrM summits would otherwise take some of the N places.
Exactly N peaks are kept (ties broken by file order); the selected lines are written unchanged and in
their original order.
"""
import logging
import re

import pandas as pd

input_peaks_file = snakemake.input.peaks
output_peaks_file = snakemake.output.top_peaks
max_peaks = int(snakemake.params.top_n_peaks)
canonical = re.compile(snakemake.params.canonical_chroms_regex)

logging.basicConfig(filename=snakemake.log[0], level=logging.INFO,
                    format='%(asctime)s %(message)s', datefmt='%Y-%m-%d %H:%M:%S')
logging.info(f"Input file: {input_peaks_file}")
logging.info(f"Output file: {output_peaks_file}")
logging.info(f"Max peaks requested (N): {max_peaks}")

with open(input_peaks_file) as f:
    lines = [line for line in f if line.strip()]

if lines:
    peaks = pd.DataFrame([line.rstrip("\n").split("\t") for line in lines]).iloc[:, [0, 6]]
    peaks.columns = ["chrom", "score"]
    peaks["score"] = peaks["score"].astype(float)

    is_canonical = peaks["chrom"].map(lambda chrom: canonical.match(chrom) is not None)
    removed = peaks.loc[~is_canonical, "chrom"].value_counts()
    logging.info(f"Found {len(peaks)} peaks; removed {(~is_canonical).sum()} on non-canonical chromosomes"
                 + (": " + ", ".join(f"{chrom} ({n})" for chrom, n in removed.items()) if len(removed) else ""))
    peaks = peaks[is_canonical]

    top = peaks.sort_values("score", ascending=False, kind="stable").head(max_peaks)
    keep = sorted(top.index)
    if len(top):
        logging.info(f"Score threshold of the selected peaks: {top['score'].min()}")
else:
    logging.warning("Input file is empty. Output file will be empty.")
    keep = []

with open(output_peaks_file, "w") as f_out:
    f_out.writelines(lines[i] for i in keep)
logging.info(f"Wrote {len(keep)} peaks to {output_peaks_file}.")
