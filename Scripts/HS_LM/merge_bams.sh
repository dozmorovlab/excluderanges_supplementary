#!/bin/bash
#SBATCH --job-name=merge_bams_101
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=48
#SBATCH --mem=16G
#SBATCH --output=merge_bams_101.out
#SBATCH --error=merge_bams_101.err

# Define the directory containing the BAM files and the output merged file
BAM_DIR="./Blacklist/bams/BAM-files-unmerged"
OUT_FILE="./Blacklist/bams/101/101_merged.bam"
ACCESSIONS_FILE="./Blacklist/bams/101/101_accessions.txt"

# Read accessions into an array
mapfile -t ACCESSIONS < "${ACCESSIONS_FILE}"

# Find BAM files that match the accessions
# Initialize BAM_FILES array
BAM_FILES=()

# Find BAM files that match the accessions
for accession in "${ACCESSIONS[@]}"; do

  MATCHES=$(find "$BAM_DIR" -name "*${accession}*.bam")

  # Check if matches were found
  if [ -n "$MATCHES" ]; then
    # Loop through each match found
    while IFS= read -r match; do
      # Add each match to the BAM_FILES array
      BAM_FILES+=("$match")
    done <<< "$MATCHES"
  else
    echo "No BAM files found matching accession: $accession"
  fi
done

# Use samtools merge to merge BAM files
samtools merge \
    -r \
    -l 9 \
    --threads 48 \
    --write-index \
    "$OUT_FILE" \
    "${BAM_FILES[@]}"
