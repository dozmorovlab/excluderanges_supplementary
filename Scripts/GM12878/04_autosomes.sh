#!/bin/bash

# Define autosomal chromosomes
AUTOSOMES=$(seq 1 22 | awk '{print "chr"$0}') # Adjust "chr" prefix based on your BAM file

# Directory containing BAM files
BAM_DIR="/lustre/home/wallbp/GM12878/bams"
OUTPUT_DIR="/lustre/home/wallbp/GM12878/bams_autosomes"

# Create output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Loop through all BAM files in the directory
for BAM in "$BAM_DIR"/*.bam; do
    BASENAME=$(basename "$BAM" .bam)
    TEMP_SAM="$OUTPUT_DIR/${BASENAME}.tmp.autosomes.sam"
    TEMP_BAM="$OUTPUT_DIR/${BASENAME}.tmp.autosomes.bam"
    OUTPUT="$OUTPUT_DIR/${BASENAME}.autosomes.bam"

    # Subset to autosomes and sex chromosomes (chr1 to chr22, chrX, chrY)
    samtools view \
        -b \
        -o "$TEMP_BAM" \
        -@ 64 \
        $BAM \
        chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY

    # Reheader: Only keep chromosomes chr1-22 and X, Y in the header
    samtools view -H "$TEMP_BAM" | grep -E '^@SQ.*SN:chr([1-9]|1[0-9]|2[0-2]|X|Y)\s+LN:[0-9]+$' > "$TEMP_SAM"

    # Check if "$TEMP_SAM" is empty and print debug info if it is
    if [ ! -s "$TEMP_SAM" ]; then
        echo "Debug: "$TEMP_SAM" is empty, printing header for debugging:"
        samtools view -H "$TEMP_BAM"
    fi

    samtools reheader "$TEMP_SAM" "$TEMP_BAM" > "$OUTPUT"

    rm "$TEMP_BAM"
    rm "$TEMP_SAM"
    samtools index -b "$OUTPUT"
done


echo "Done! Subset BAMs are in $OUTPUT_DIR"
