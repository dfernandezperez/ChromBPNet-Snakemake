"""Run chrombpnet's DeepSHAP interpretation with several processes sharing one GPU.

chrombpnet interprets one region at a time, building the shuffled references in python, so a single
process leaves the GPU idle part of the time (~70% busy on an RTX 3090; 2 processes are ~1.2x faster,
more just compete for the GPU and each needs ~5 GB of its memory). The regions are split into
`--n-procs` contiguous parts, interpreted in parallel (TF_FORCE_GPU_ALLOW_GROWTH lets the processes
share the GPU memory) and concatenated in the original order. Outputs are the same as `chrombpnet contribs_bw` without bigwigs:
  <out-prefix>.interpreted_regions.bed and <out-prefix>.{profile,counts}_scores.h5
Must run inside the chrombpnet image.

Usage: python parallel_contribs.py --model nobias.h5 --regions peaks.bed --genome genome.fa \
           --out-prefix out/prefix [--heads profile counts] [--n-procs 4] [--threads 4]
"""
import argparse
import os
import shutil
import subprocess
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from contribs_utils import concat_beds, concat_h5, count_lines, split_lines  # noqa: E402


def tail(path, n=20):
    with open(path) as f:
        return "".join(f.readlines()[-n:])


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--model", required=True, help="Model h5 (e.g. *_chrombpnet_nobias.h5 or bias.h5)")
    parser.add_argument("--regions", required=True, help="10 column bed file (narrowPeak) of regions, centered at the summit")
    parser.add_argument("--genome", required=True, help="Genome fasta")
    parser.add_argument("--out-prefix", required=True)
    parser.add_argument("--heads", nargs="+", default=["counts", "profile"], choices=["counts", "profile"])
    parser.add_argument("--n-procs", type=int, default=1, help="Parallel DeepSHAP processes")
    parser.add_argument("--threads", type=int, default=1, help="CPUs available, shared by the processes")
    parser.add_argument("--keep-parts", action="store_true", help="Keep the per-process files")
    args = parser.parse_args()

    n_regions = count_lines(args.regions)
    if n_regions == 0:
        sys.exit(f"No regions in {args.regions}")
    n_procs = max(1, min(args.n_procs, n_regions))
    threads = str(max(1, args.threads // n_procs))

    work = args.out_prefix + ".parts"
    shutil.rmtree(work, ignore_errors=True) # Leftovers of a failed or preempted run
    os.makedirs(work)
    parts = [os.path.join(work, f"part_{i}") for i in range(n_procs)]
    split_lines(args.regions, n_procs, [p + ".bed" for p in parts])

    env = dict(os.environ, TF_FORCE_GPU_ALLOW_GROWTH="true", TF_CPP_MIN_LOG_LEVEL="1",
               OMP_NUM_THREADS=threads, OPENBLAS_NUM_THREADS=threads, MKL_NUM_THREADS=threads,
               NUMBA_NUM_THREADS=threads, TF_NUM_INTRAOP_THREADS=threads, TF_NUM_INTEROP_THREADS="1")
    print(f"Interpreting {n_regions} regions ({', '.join(args.heads)}) with {n_procs} processes "
          f"x {threads} threads", flush=True)
    start_time = time.time()
    procs = []
    for part in parts:
        cmd = [sys.executable, "-m", "chrombpnet.evaluation.interpret.interpret",
               "-g", args.genome, "-r", part + ".bed", "-m", args.model, "-o", part, "-p", *args.heads]
        with open(part + ".log", "w") as log:
            procs.append(subprocess.Popen(cmd, stdout=log, stderr=subprocess.STDOUT, env=env))

    # Fail as soon as one process fails
    while True:
        codes = [p.poll() for p in procs]
        failed = [i for i, code in enumerate(codes) if code not in (None, 0)]
        if failed:
            for p in procs:
                if p.poll() is None:
                    p.terminate()
            for i in failed:
                print(f"--- Process {i} failed (exit code {codes[i]}). End of {parts[i]}.log:\n{tail(parts[i] + '.log')}",
                      file=sys.stderr)
            sys.exit(1)
        if all(code == 0 for code in codes):
            break
        time.sleep(10)

    minutes = (time.time() - start_time) / 60
    print(f"Interpretation finished in {minutes:.1f} min ({n_regions / minutes:.1f} regions/min)", flush=True)
    for i, part in enumerate(parts):
        print(f"--- End of process {i} log:\n{tail(part + '.log', 5)}", flush=True)

    beds = [part + ".interpreted_regions.bed" for part in parts]
    for head in args.heads:
        concat_h5([f"{part}.{head}_scores.h5" for part in parts], beds, f"{args.out_prefix}.{head}_scores.h5")
    concat_beds(beds, args.out_prefix + ".interpreted_regions.bed")

    if not args.keep_parts:
        shutil.rmtree(work)
    print("Done", flush=True)


if __name__ == "__main__":
    main()
