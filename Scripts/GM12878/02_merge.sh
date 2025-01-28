#!/bin/bash

# Input TSV file (hardcoded for this example)
input_tsv="./metadata.tsv"

# Output folder
out_folder="./merged_fastqs"
mkdir -p "${out_folder}"

# threads
threads=2

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

# Iterate over the experiment targets and merge the associated files
for target in "${!target_files[@]}"; do
  # Trim any trailing whitespace from the target_files entry
  file_list=$(echo "${target_files["$target"]}" | xargs)

  output_file="${out_folder}/${target}_merged.fastq.gz"
  echo "Merging files for target: $target into $output_file"

  (for file in $file_list; do
    bgzip -cd "$file"
  done) | bgzip -c -@ "${threads}" > "$output_file"
done

echo "Merging completed."
