#!/bin/bash

# Input TSV file (hardcoded for this example)
input_tsv="./metadata.tsv"

# STAR reference genome directory
reference_path="./star_sponge"

# Number of threads for parallel processing
threads=64

# Output directory for BAM files
bam_dir="./bams_sponge"
mkdir -p "${bam_dir}"

# Create an associative array to store file accessions grouped by experiment target
declare -A target_files

# Read the TSV file line by line (ignoring the header)
while IFS=$'\t' read -r file_accession experiment_target; do
  # Remove leading/trailing whitespace
  experiment_target=$(echo "$experiment_target" | xargs)
  file_accession=$(echo "$file_accession" | xargs)

  # Append the file path to the target's array
  target_files["$experiment_target"]+="./${file_accession}.fastq.gz "

done < <(tail -n +2 "$input_tsv" | awk -F'\t' '{ print $1 "\t" $23 }')

# Iterate over the experiment targets and perform the alignment
for target in "${!target_files[@]}"; do

    echo "Aligning ${target}"

    # Path to the merged FASTQ file
    in_file="./merged_fastqs/${target}_merged.fastq.gz"

    # Output prefix for STAR
    output_prefix="${bam_dir}/${target}_sponge_"

    # Run STAR alignment
    STAR --genomeDir "${reference_path}" \
        --readFilesIn "${in_file}" \
        --outSAMtype BAM SortedByCoordinate \
        --outFileNamePrefix "${output_prefix}" \
        --alignEndsType Local \
        --runThreadN ${threads} \
        --readFilesCommand bgzip -cd \
        --outBAMcompression 10
done
