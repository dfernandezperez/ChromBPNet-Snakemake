"""Helpers for chrombpnet contribution score files, used by parallel_contribs.py and merge_contribs.py.

`chrombpnet contribs_bw` (and its interpret module) writes, for an output prefix:
  <prefix>.interpreted_regions.bed   regions that could be interpreted, in input order
  <prefix>.{profile,counts}_scores.h5  raw/seq (one-hot), shap/seq (hypothetical contributions) and
                                       projected_shap/seq (raw * shap), shape (N, 4, input_len)
  <prefix>.{profile,counts}_scores.bw  projected contributions as a bigwig
Splitting the regions, interpreting the parts and concatenating them in order gives the same files.
Must run inside the chrombpnet image (bigwig writing uses chrombpnet's own helper).
"""
import os
import shutil

import h5py
import hdf5plugin  # registers the blosc filter used by chrombpnet (deepdish) to write the h5 files
import numpy as np
import pyBigWig

KEYS = ["raw/seq", "shap/seq", "projected_shap/seq"]
BLOCK_ROWS = 1000
# Same filter as the files written by chrombpnet (deepdish), so any of its readers can open them
COMPRESSION = hdf5plugin.Blosc(cname="lz4", clevel=5, shuffle=hdf5plugin.Blosc.SHUFFLE)


def count_lines(path):
    with open(path) as f:
        return sum(1 for line in f if line.strip())


def split_lines(path, n_parts, out_paths):
    """Split a file into `n_parts` contiguous blocks of lines of (almost) equal size."""
    with open(path) as f:
        lines = [line for line in f if line.strip()]
    for part, out_path in zip(np.array_split(np.arange(len(lines)), n_parts), out_paths):
        with open(out_path, "w") as f:
            f.writelines(lines[i] for i in part)


def concat_beds(beds, out_bed):
    tmp = out_bed + ".tmp"
    with open(tmp, "w") as out:
        for bed in beds:
            with open(bed) as f:
                shutil.copyfileobj(f, out)
    os.replace(tmp, out_bed)


def concat_h5(h5s, beds, out_h5):
    """Concatenate contribution h5 files row-wise, in the given order.

    `beds` are the interpreted_regions files of each part, used to check that every h5 has one row
    per interpreted region.
    """
    inputs = [h5py.File(path, "r") for path in h5s]
    try:
        shape = inputs[0]["raw/seq"].shape[1:]
        dtypes = {key: inputs[0][key].dtype for key in KEYS}
        n_rows = 0
        for path, bed, f in zip(h5s, beds, inputs):
            n_bed = count_lines(bed)
            for key in KEYS:
                if f[key].shape[1:] != shape or f[key].shape[0] != n_bed:
                    raise ValueError(f"{path}:{key} has shape {f[key].shape}, expected ({n_bed}, {shape})")
            n_rows += n_bed

        tmp = out_h5 + ".tmp"
        with h5py.File(tmp, "w") as out:
            for key in KEYS:
                if n_rows == 0:
                    out.create_dataset(key, shape=(0,) + shape, dtype=dtypes[key])
                    continue
                dset = out.create_dataset(key, shape=(n_rows,) + shape, dtype=dtypes[key],
                                          chunks=(min(BLOCK_ROWS, n_rows),) + shape, **COMPRESSION)
                start = 0
                for f in inputs:
                    n = f[key].shape[0]
                    for offset in range(0, n, BLOCK_ROWS):
                        block = f[key][offset:offset + BLOCK_ROWS]
                        dset[start + offset:start + offset + len(block)] = block
                    start += n
        os.replace(tmp, out_h5)
    finally:
        for f in inputs:
            f.close()


def write_contribs_bigwig(h5, bed, chrom_sizes, out_bw):
    """Bigwig of the projected contributions, written as `chrombpnet contribs_bw` does
    (chrombpnet.evaluation.make_bigwigs.importance_hdf5_to_bigwig)."""
    import chrombpnet.evaluation.make_bigwigs.bigwig_helper as bigwig_helper

    gs = bigwig_helper.read_chrom_sizes(chrom_sizes)
    with h5py.File(h5, "r") as f:
        dset = f["projected_shap/seq"]
        n_rows, _, seqlen = dset.shape
        data = np.empty((n_rows, seqlen), dtype=dset.dtype)
        for start in range(0, n_rows, BLOCK_ROWS):
            data[start:start + BLOCK_ROWS] = dset[start:start + BLOCK_ROWS].sum(1)

    tmp = out_bw + ".tmp"
    if n_rows == 0: # write_bigwig fails without regions
        bw = pyBigWig.open(tmp, "w")
        bw.addHeader(gs)
        bw.close()
    else:
        regions = bigwig_helper.get_regions(bed, seqlen)
        bigwig_helper.write_bigwig(data, regions, gs, tmp)
    os.replace(tmp, out_bw)
