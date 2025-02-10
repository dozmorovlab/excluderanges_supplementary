#!/bin/bash
#SBATCH --job-name=MACS3_call
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=48
#SBATCH --mem=64G
#SBATCH --output=MACS3_call.out
#SBATCH --error=MACS3_call.err

in_file="./Blacklist/bams/101/101_merged.bam"
out_dir="./Blacklist/bams/101/101_MACS3"
k=1000

mkdir -p "${out_dir}"

source ~/miniconda3/etc/profile.d/conda.sh
conda activate macs3

macs3 callpeak \
    --treatment "${in_file}" \
    --name "101_local" \
    --outdir "${out_dir}" \
    --gsize hs \
    --slocal 10000\
    --llocal 100000\
    --keep-dup all

macs3 callpeak \
    --treatment "${in_file}" \
    --name "101_global" \
    --outdir "${out_dir}"\
    --gsize hs \
    --nolambda \
    --keep-dup all

bedtools merge \
    -i "${out_dir}/101_local_peaks.narrowPeak" \
    -d "${k}" \
| awk -v k="${k}" \
    '($3 - $2) >= min_size' \
    > "${out_dir}/101_local_peaks.merged.filtered.bed"

bedtools merge \
    -i "${out_dir}/101_peaks.narrowPeak" \
    -d "${k}" \
| awk -v k="${k}" \
    '($3 - $2) >= k' \
    > "${out_dir}/101_peaks.merged.filtered.bed"
