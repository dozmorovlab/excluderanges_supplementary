#!/bin/bash

# conda
source ~/miniconda3/etc/profile.d/conda.sh
conda activate macs3_env

# Input
bams="./bams_sponge"

# Output directories
out_MACS3_dir="./MACS3"
out_bed_dir="./beds_sponge"
k=0

threads=64
mkdir -p "${out_MACS3_dir}"
mkdir -p "${out_bed_dir}"

# Input TSV file
input_tsv="./metadata.tsv"

# Create an associative array to store file accessions grouped by experiment target
declare -A targets

# Read the TSV file line by line (ignoring the header)
while IFS=$'\t' read -r file_accession experiment_target; do
  # Remove leading/trailing whitespace
  experiment_target=$(echo "$experiment_target" | xargs)
  file_accession=$(echo "$file_accession" | xargs)

  # Append the file path to the target's array
  targets["$experiment_target"]+="${file_accession} "
  echo "${file_accession}, ${experiment_target}"

done < <(tail -n +2 "$input_tsv" | awk -F'\t' '{ print $1 "\t" $23 }')

# Iterate over the experiment targets and call peaks
for target in "${!targets[@]}"; do
    # Correct BAM file path based on the STAR output
    bam="${bams}/${target}_sponge_Aligned.sortedByCoord.out.bam"
    echo "Processing target: ${target} with BAM: ${bam}"

    # Check if BAM file exists before proceeding
    if [[ ! -f "${bam}" ]]; then
        echo "BAM file not found for target ${target}, skipping..."
        continue
    fi

    # Create output directory for MACS3 results for each target
    macs3_output_dir="${out_MACS3_dir}_${target}"
    mkdir -p "${macs3_output_dir}"

    # Call peaks with MACS3
    macs3 callpeak \
        --treatment "${bam}" \
        --name "${target}_global" \
        --outdir "${macs3_output_dir}" \
        --gsize hs \
        --nolambda \
        --keep-dup all

    macs3 callpeak \
        --treatment "${bam}" \
        --name "${target}_local" \
        --outdir "${macs3_output_dir}" \
        --gsize hs \
        --slocal 10000\
        --llocal 100000\
        --keep-dup all

    global_narrowpeak_file="${macs3_output_dir}/${target}_global_peaks.narrowPeak"
    local_narrowpeak_file="${macs3_output_dir}/${target}_local_peaks.narrowPeak"

    # Merge peaks using bedtools
    #bedtools merge \
    #    -i "${global_narrowpeak_file}" \
    #    -d "${k}" \
    #| awk -v k="${k}" '($3 - $2) >= k' \
    #> "${out_bed_dir}/${target}_peaks.merged.filtered.bed"

    # Merge peaks using bedtools
    #bedtools merge \
    #    -i "${local_narrowpeak_file}" \
    #    -d "${k}" \
    #| awk -v k="${k}" '($3 - $2) >= k' \
    #> "${out_bed_dir}/${target}_peaks.merged.filtered.bed"

done
