#!/bin/bash

# conda
source ~/miniconda3/etc/profile.d/conda.sh
conda activate deeptools

# Input
ex_list_dir="./ex_lists"

# Output directories
out_dir="/lustre/home/wallbp/GM12878/summary_101"

# Threads
threads=60

mkdir -p "${out_dir}"

echo "$(date '+%Y-%m-%d %H:%M:%S'): Getting summary: all 101"

multiBamSummary bins \
    --bamfiles "/lustre/home/juicer/Blacklist/bams/bams/star_100.bam" "/lustre/home/juicer/Blacklist/bams/star_100_sponge.sorted.bam" \
    --outFileName "${out_dir}/101.npz" \
    --numberOfProcessors max \
    --outRawCounts "${out_dir}/101.counts.tsv"

for ex_list in "$ex_list_dir"/*; do
    ex_base="$(basename "$ex_list" .bed)"

    echo "$(date '+%Y-%m-%d %H:%M:%S'): Getting summary: 101 + $ex_base"

    multiBamSummary BED-file \
        --bamfiles "/lustre/home/juicer/Blacklist/bams/bams/star_100.bam" "/lustre/home/juicer/Blacklist/bams/star_100_sponge.sorted.bam" \
        --outFileName "${out_dir}/101.${ex_base}.ex.npz" \
        --numberOfProcessors max \
        --outRawCounts "${out_dir}/101.${ex_base}.counts.ex.tsv" \
        --BED "$ex_list"

done
