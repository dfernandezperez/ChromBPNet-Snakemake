"""Average ChromBPNet contribution score h5 files (chrombpnet contribs_bw) across folds.

The h5 files store, for every region of `<prefix>.interpreted_regions.bed`:
  raw/seq             one-hot sequence           (N, 4, input_len)
  shap/seq            hypothetical contributions (N, 4, input_len)
  projected_shap/seq  raw * shap                 (N, 4, input_len)
Rows of the different folds must correspond to the same regions, which is checked with the
interpreted_regions bed files. raw/seq is copied and both score datasets are averaged.

Usage: python average_contribs_h5.py --h5 f0.h5 f1.h5 --beds f0.bed f1.bed --out-h5 mean.h5 [--out-bed mean.bed]
"""
import argparse
import filecmp
import shutil

import h5py
import hdf5plugin  # noqa: F401 - registers the blosc filter used by chrombpnet (deepdish) to write the h5 files
import numpy as np

SCORE_KEYS = ["shap/seq", "projected_shap/seq"]
CHUNK_ROWS = 5000


def check_same_regions(beds):
    for bed in beds[1:]:
        if not filecmp.cmp(beds[0], bed, shallow=False):
            raise ValueError(f"Interpreted regions differ between {beds[0]} and {bed}; cannot average folds row by row.")


def average_contribs(input_h5s, output_h5):
    inputs = [h5py.File(path, "r") for path in input_h5s]
    shape = inputs[0]["raw/seq"].shape
    for path, f in zip(input_h5s, inputs):
        for key in ["raw/seq"] + SCORE_KEYS:
            if f[key].shape != shape:
                raise ValueError(f"{path}:{key} has shape {f[key].shape}, expected {shape}")

    with h5py.File(output_h5, "w") as out:
        out.create_dataset("raw/seq", data=inputs[0]["raw/seq"][:], compression="lzf")
        for key in SCORE_KEYS:
            dset = out.create_dataset(key, shape=shape, dtype=inputs[0][key].dtype,
                                      chunks=(min(CHUNK_ROWS, shape[0]),) + shape[1:], compression="lzf")
            # Accumulate in float32 in chunks of rows to bound memory
            for start in range(0, shape[0], CHUNK_ROWS):
                end = min(start + CHUNK_ROWS, shape[0])
                total = np.zeros((end - start,) + shape[1:], dtype=np.float32)
                for f in inputs:
                    total += f[key][start:end].astype(np.float32)
                dset[start:end] = (total / len(inputs)).astype(dset.dtype)

    for f in inputs:
        f.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--h5", nargs="+", required=True, help="Contribution h5 files, one per fold")
    parser.add_argument("--beds", nargs="+", required=True, help="interpreted_regions.bed files, same order as --h5")
    parser.add_argument("--out-h5", required=True, help="Output averaged h5")
    parser.add_argument("--out-bed", help="Optional copy of the (shared) interpreted regions bed")
    args = parser.parse_args()

    if len(args.h5) != len(args.beds):
        parser.error("--h5 and --beds need the same number of files")
    check_same_regions(args.beds)
    average_contribs(args.h5, args.out_h5)
    if args.out_bed:
        shutil.copyfile(args.beds[0], args.out_bed)
