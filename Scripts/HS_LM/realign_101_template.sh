#!/bin/bash
#SBATCH --job-name=realign_101_BASENAME
#SBATCH --output=./%j.out
#SBATCH --error=./%j.err
#SBATCH --mem=64G
#SBATCH --cpus-per-task=16

# Load necessary modules
# conda create -n align
# conda activate align
# conda install bwa-mem2 star bowtie2 samtools

# Initialize Conda
source ~/miniconda3/etc/profile.d/conda.sh
conda activate align

R1="OUT_DIR/BASENAME_1.fastq.gz"
R2="OUT_DIR/BASENAME_2.fastq.gz"

# Convert BAM to FASTQ
samtools fastq \
    -1  "${R1}" \
    -2  "${R2}" \
    -0 /dev/null \
    -s /dev/null \
    -n \
    -O \
    -c 9 \
    -@ 16 \
    BAM_FILE

bwa-mem2 mem \
    -t 16 \
    REF_PATH \
    "${R1}" \
    "${R2}" \
| samtools view -@ 16 -bS - \
| samtools sort -@ 16 -o OUT_DIR/BASENAME_101_bwamem2.bam -
samtools index -@ 16 OUT_DIR/BASENAME_101_bwamem2.bam

bwt2_x="$(dirname REF_PATH)/$(basename REF_PATH .fna.gz)"

bowtie2 \
    -x "${bwt2_x}" \
    -U \
    "${R1}" \
    "${R2}" \
    -p 16 \
    --local \
| samtools view -@ 16 -bS - \
| samtools sort -@ 16 -o OUT_DIR/BASENAME_101_bowtie2.bam -
samtools index -@ 16 OUT_DIR/BASENAME_101_bowtie2.bam

STAR \
    --genomeDir STAR_PATH \
    --readFilesIn <(bgzip -cd -@ 16 "${R1}") \
        <(bgzip -cd -@ 16 "${R2}") \
    --outSAMtype BAM SortedByCoordinate \
    --outFileNamePrefix OUT_DIR/BASENAME_101_star \
    --alignEndsType Local \
    --runThreadN 16 \
    --limitBAMsortRAM 64000000000
samtools index -@ 16 OUT_DIR/BASENAME_101_starAligned.sortedByCoord.out.bam
