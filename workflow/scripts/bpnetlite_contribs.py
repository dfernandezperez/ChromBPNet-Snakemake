"""Contribution scores (DeepLIFT/SHAP) of a chrombpnet model with bpnet-lite's PyTorch API.

Drop-in replacement of parallel_contribs.py (chrombpnet's TF DeepSHAP), with the same inputs and
outputs (<out-prefix>.interpreted_regions.bed and <out-prefix>.{profile,counts}_scores.h5 with
raw/seq, shap/seq (hypothetical contributions) and projected_shap/seq). The computation is
bpnet-lite/tangermeme except for the profile target (ProfileTarget below):
  - BPNet.from_chrombpnet loads the chrombpnet TF model (*_chrombpnet_nobias.h5; bpnet-lite's
    output crops match chrombpnet's for this 8-layer model, not for the 4-layer bias model)
  - targets: bpnet-lite's CountWrapper (predicted log counts) and ProfileTarget (chrombpnet's
    profile target), not bpnet-lite's ProfileWrapper, whose scores are unstable
  - bpnetlite.attribute.deep_lift_shap with `--n-shuffles` dinucleotide-shuffled references per
    region, seeded with `--seed` (reproducible, unlike chrombpnet's)
  - tangermeme.io.extract_loci reads the sequences centered at the summits (start + 10th column)
    and skips regions whose input window falls off the chromosome, as chrombpnet does
Regions with N are skipped, as bpnet-lite's `bpnet attribute` does (tangermeme's shuffles need
one-hot sequences).

Usage: python bpnetlite_contribs.py --model nobias.h5 --regions peaks.bed --genome genome.fa \
           --out-prefix out/prefix [--heads counts profile] [--batch-size 64] [--n-shuffles 20] \
           [--seed 1234] [--dtype float32]
"""
import argparse
import os
import time
import warnings

import h5py
import hdf5plugin  # Blosc filter, the same compression as the h5 files of chrombpnet (deepdish)
import numpy as np
import pandas as pd
import torch
from bpnetlite import BPNet
from bpnetlite.attribute import deep_lift_shap
from bpnetlite.bpnet import CountWrapper
from tangermeme.io import extract_loci

NARROWPEAK_SCHEMA = ["chr", "start", "end", "1", "2", "3", "4", "5", "6", "summit"]
COMPRESSION = hdf5plugin.Blosc(cname="lz4", clevel=5, shuffle=hdf5plugin.Blosc.SHUFFLE)
BLOCK_ROWS = 2000 # Regions interpreted and written at a time
NON_ACGT = list("QWERYUIOPSDFHJKLZXVBNM") # Letters one-hot encoded as all zeros (as `bpnet attribute`)


class ProfileTarget(torch.nn.Module):
    """chrombpnet's profile target: the mean-normalized profile logits weighted by their softmax,
    sum(logits * softmax(stop_gradient(logits))), with the multipliers chrombpnet's DeepSHAP gives it.

    bpnet-lite's ProfileWrapper registers logits * softmax(logits) as one element-wise non-linearity,
    so tangermeme's rescale rule divides each position's change by the change of its own logit. The
    softmax makes a position change even when its logit barely does, so the ratio explodes at those
    positions: the scores still sum to the change in prediction (no convergence warning), but they
    agree with chrombpnet's at r ~0.7 and between seeds at ~0.46 (chrombpnet: 0.85).

    chrombpnet's DeepSHAP (shap's TF DeepExplainer) has no handler for StopGradient, so it attributes
    the product with its two-input rule: the logits get the multiplier (w_input + w_reference) / 2.
    The same is done here by averaging the weights of each input-reference pair, which relies on
    deep_lift_shap running a batch of inputs concatenated with their references.
    """

    def __init__(self, model):
        super().__init__()
        self.model = model

    def forward(self, X, X_ctl=None):
        if X.shape[0] % 2:
            raise ValueError("ProfileTarget expects inputs concatenated with their references")
        logits = self.model(X)[0].reshape(X.shape[0], -1)
        logits = logits - logits.mean(dim=-1, keepdim=True)
        w_input, w_reference = torch.softmax(logits.detach().float(), dim=-1).chunk(2)
        weights = ((w_input + w_reference) / 2).repeat(2, 1)
        return (logits * weights).sum(dim=-1, keepdim=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--model", required=True, help="chrombpnet model h5 (*_chrombpnet_nobias.h5)")
    parser.add_argument("--regions", required=True, help="10 column bed file (narrowPeak) of regions, centered at the summit")
    parser.add_argument("--genome", required=True, help="Genome fasta")
    parser.add_argument("--out-prefix", required=True)
    parser.add_argument("--heads", nargs="+", default=["counts", "profile"], choices=["counts", "profile"])
    parser.add_argument("--input-length", type=int, default=2114, help="Input length of the model (chrombpnet: 2114)")
    parser.add_argument("--n-shuffles", type=int, default=20, help="Dinucleotide-shuffled references per region (chrombpnet: 20)")
    parser.add_argument("--batch-size", type=int, default=64, help="Region-reference pairs per GPU batch")
    parser.add_argument("--seed", type=int, default=1234, help="Seed of the shuffled references")
    parser.add_argument("--dtype", default="float32", choices=["float32", "float16"],
                        help="dtype the model is cast to (as `bpnet attribute` does), which tangermeme computes in")
    args = parser.parse_args()

    device = "cuda" if torch.cuda.is_available() else "cpu"
    dtype = getattr(torch, args.dtype)
    model = BPNet.from_chrombpnet(args.model).to(device=device, dtype=dtype).eval()
    # float16 rounds the predicted log counts (~5) by up to ~0.01, so the convergence deltas of the
    # counts head reach ~0.01 with identical scores (r = 1.0000 vs float32); a missing rule gives 0.1-1
    warning_threshold = 2e-2 if dtype == torch.float16 else 1e-3
    targets = {"counts": CountWrapper(model), "profile": ProfileTarget(model)}

    X, kept = extract_loci(loci=args.regions, sequences=args.genome, in_window=args.input_length,
                           summits=True, ignore=NON_ACGT, return_mask=True)
    no_n = X.sum(dim=(1, 2)) == X.shape[-1]
    X = X[no_n]
    kept[kept.clone()] = no_n
    regions = pd.read_csv(args.regions, sep="\t", names=NARROWPEAK_SCHEMA)[kept.numpy()]
    n_regions = len(X)
    if n_regions == 0:
        raise SystemExit(f"No region of {args.regions} can be interpreted")
    print(f"Interpreting {n_regions} regions ({int((~no_n).sum())} with N skipped; {', '.join(args.heads)}) on {device}, "
          f"{args.n_shuffles} references, batch size {args.batch_size}, {args.dtype}", flush=True)

    # Written to temporary files and renamed at the end, so a killed (e.g. preempted) job leaves no
    # complete-looking output
    tmp = {head: f"{args.out_prefix}.{head}_scores.h5.tmp" for head in args.heads}
    outputs = {head: h5py.File(path, "w") for head, path in tmp.items()}
    shape = (n_regions, 4, args.input_length)
    dsets = {head: {key: f.create_dataset(key, shape=shape, dtype=h5_dtype, chunks=(min(500, n_regions), 4, args.input_length), **COMPRESSION)
                    for key, h5_dtype in [("raw/seq", np.int8), ("shap/seq", np.float16), ("projected_shap/seq", np.float16)]}
             for head, f in outputs.items()}

    start_time = time.time()
    for block_start in range(0, n_regions, BLOCK_ROWS):
        block_end = min(block_start + BLOCK_ROWS, n_regions)
        X_block = X[block_start:block_end].type(dtype)
        for head in args.heads:
            with warnings.catch_warnings():
                if head == "profile":
                    # With constant softmax weights the profile target is not a sum of DeepLIFT
                    # contributions, so its convergence deltas are not ~0 (nor are chrombpnet's)
                    warnings.filterwarnings("ignore", message="Convergence deltas too high")
                hyp = deep_lift_shap(targets[head], X_block, hypothetical=True, n_shuffles=args.n_shuffles,
                                     batch_size=args.batch_size, random_state=args.seed, device=device,
                                     warning_threshold=warning_threshold).float().numpy()
            raw = X_block.numpy().astype(np.int8)
            dsets[head]["raw/seq"][block_start:block_end] = raw
            dsets[head]["shap/seq"][block_start:block_end] = hyp.astype(np.float16)
            dsets[head]["projected_shap/seq"][block_start:block_end] = (hyp * raw).astype(np.float16)
        minutes = (time.time() - start_time) / 60
        print(f"{block_end}/{n_regions} regions ({block_end / minutes:.1f} regions/min)", flush=True)

    for f in outputs.values():
        f.close()
    for head, path in tmp.items():
        os.replace(path, f"{args.out_prefix}.{head}_scores.h5")
    # Same text as chrombpnet's interpret writes, so outputs of both implementations can be merged
    regions.to_csv(args.out_prefix + ".interpreted_regions.bed", header=False, index=False, sep="\t")
    minutes = (time.time() - start_time) / 60
    print(f"Interpretation finished in {minutes:.1f} min ({n_regions / minutes:.1f} regions/min)", flush=True)
    if device == "cuda":
        print(f"Peak GPU memory: {torch.cuda.max_memory_allocated() / 2**30:.1f} GiB", flush=True)


if __name__ == "__main__":
    main()
