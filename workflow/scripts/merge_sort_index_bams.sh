#!/usr/bin/env bash
# Merge (if there are several), coordinate-sort and index BAM files.
# Usage: merge_sort_index_bams.sh <threads> <output.bam> <input.bam>...
# Sorting is skipped when all the inputs are already coordinate-sorted (@HD SO:coordinate): a single
# input is symlinked and several are merged with `samtools merge`, which keeps the order.
set -euo pipefail

threads=$1
out=$2
shift 2

all_sorted=true
for bam in "$@"; do
    if ! samtools view -H "$bam" | awk '/^@HD/ && /SO:coordinate/ {found = 1} END {exit !found}'; then
        echo "$bam is not coordinate-sorted"
        all_sorted=false
    fi
done

rm -f "$out" "$out.bai"
if $all_sorted && [ $# -eq 1 ]; then
    echo "Input is coordinate-sorted: linking it"
    ln -s "$(realpath "$1")" "$out"
elif $all_sorted; then
    echo "Inputs are coordinate-sorted: merging $# BAM files"
    samtools merge -@ "$threads" -f "$out" "$@"
elif [ $# -eq 1 ]; then
    echo "Sorting the input BAM"
    samtools sort -@ "$threads" -m 7G -o "$out" "$1"
else
    echo "Merging $# BAM files and sorting"
    samtools merge -@ "$threads" -f "$out.unsorted.bam" "$@"
    samtools sort -@ "$threads" -m 7G -o "$out" "$out.unsorted.bam"
    rm "$out.unsorted.bam"
fi

echo "Indexing final BAM"
samtools index -@ "$threads" "$out"
echo "BAM processing finished"
