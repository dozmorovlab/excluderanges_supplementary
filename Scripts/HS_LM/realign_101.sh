#!/bin/bash
#SBATCH --job-name=realign_101
#SBATCH --output=./%j.out
#SBATCH --error=./%j.err
#SBATCH --mem=1G
#SBATCH --cpus-per-task=1

# BAM directory path
BAM_DIR="./Blacklist/bams/BAM-files-unmerged"

REMOVE_NAME=".bam.sorted.bam"

OUT_DIR="./realigned_101"

mkdir -p $OUT_DIR

# List of accessions to use
ACCESSIONS_FILE="./Blacklist/bams/101/101_accessions.txt"
#ACCESSIONS_FILE="./Blacklist/bams/101/test.txt"

# Job script template
JOB_SCRIPT_TEMPLATE="realign_101_template.sh"

# Reference URL
REF_URL="https://www.encodeproject.org/files/GRCh38_no_alt_analysis_set_GCA_000001405.15/@@download/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta.gz"
REF_PATH="./Blacklist/references/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz"
STAR_PATH="./Blacklist/references/star"

if [ ! -f "$REF_PATH" ]; then
  curl -o "$REF_PATH" "$REF_URL"
fi

if [ ! -d "$STAR_PATH" ]; then
  echo "Index reference for STAR before resuming"
  exit 1
fi

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

# Check if directory exists
if [ ! -d "$BAM_DIR" ]; then
  echo "Directory $BAM_DIR does not exist."
  exit 1
fi

# Loop over each BAM file in the directory
for BAM_FILE in "${BAM_FILES[@]}"; do
  # Extract the filename without the path and extension
  BASENAME=$(basename "${BAM_FILE}" "${REMOVE_NAME}")
  
  # Create a unique job script for each BAM file
  JOB_SCRIPT="job_$BASENAME.sh"
  cp "$JOB_SCRIPT_TEMPLATE" "$JOB_SCRIPT"
  
  # Replace placeholders in the job script with the actual file names
  sed -i "s|BAM_FILE|$BAM_FILE|g" "$JOB_SCRIPT"
  sed -i "s|OUT_DIR|$OUT_DIR|g" "$JOB_SCRIPT"
  sed -i "s|BASENAME|$BASENAME|g" "$JOB_SCRIPT"
  sed -i "s|REF_PATH|$REF_PATH|g" "$JOB_SCRIPT"
  sed -i "s|STAR_PATH|$STAR_PATH|g" "$JOB_SCRIPT"
  
  # Submit the job
  sbatch "$JOB_SCRIPT"

done
