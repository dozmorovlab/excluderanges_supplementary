#!/bin/bash
#SBATCH --job-name=realign_36_BASENAME
#SBATCH --output=./%j.out
#SBATCH --error=./%j.err
#SBATCH --mem=96G
#SBATCH --cpus-per-task=24
#SBATCH --partition=cpu

# Load necessary modules
# conda create -n align
# conda activate align
# conda install bwa-mem2 star bowtie2 samtools

# Initialize Conda
source ~/miniconda3/etc/profile.d/conda.sh
conda activate align

R0="OUT_DIR/BASENAME.fastq.gz"

# Convert BAM to FASTQ
samtools fastq \
    -0  "${R0}" \
    -n \
    -O \
    -c 9 \
    -@ 24 \
    BAM_FILE

bwa-mem2 mem \
    -t 24 \
    REF_PATH \
    "${R0}" \
| samtools view -@ 24 -bS - \
| samtools sort -@ 24 -o OUT_DIR/BASENAME_36_bwamem2.bam -
samtools index -@ 24 OUT_DIR/BASENAME_36_bwamem2.bam

bwt2_x="$(dirname REF_PATH)/$(basename REF_PATH .fna.gz)"

bowtie2 \
    -x "${bwt2_x}" \
    -U \
    "${R0}" \
    -p 24 \
    --local \
| samtools view -@ 24 -bS - \
| samtools sort -@ 24 -o OUT_DIR/BASENAME_36_bowtie2.bam -
samtools index -@ 24 OUT_DIR/BASENAME_36_bowtie2.bam

STAR \
    --genomeDir STAR_PATH \
    --readFilesIn <(bgzip -cd -@ 24 "${R0}") \
    --outSAMtype BAM SortedByCoordinate \
    --outFileNamePrefix OUT_DIR/BASENAME_36_star \
    --alignEndsType Local \
    --runThreadN 24 \
    --limitBAMsortRAM 64000000000
samtools index -@ 24 OUT_DIR/BASENAME_36_starAligned.sortedByCoord.out.bam
