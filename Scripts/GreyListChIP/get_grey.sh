#!/usr/bin/env bash
#SBATCH --job-name=get_grey         # Job name
#SBATCH --error=./get_grey_errors.txt  # Standard error file
#SBATCH --output=./get_grey_output.txt # Standard output file
#SBATCH --array=1-6                # Array of jobs
#SBATCH --mem=16G

names=(
    ""
    "bwamem2_36"
    "bwamem2_100"
    "bowtie2_36"
    "bowtie2_100"
    "star_36"
    "star_100"
)
max_gap=1000
current_bam="${names[$SLURM_ARRAY_TASK_ID]}.bam"
current_bed="${names[$SLURM_ARRAY_TASK_ID]}_${max_gap}.bed"

cd /lustre/home/juicer/wallbp/

# Load the conda environment
source ~/miniconda3/etc/profile.d/conda.sh
conda activate R

# Run the R script
Rscript use_grey.R \
    --new_func \
    --bed $current_bed \
    --bam $current_bam \
    --kary "hg38_karyotype_no_alt.txt" \
    --n_cores 8 \
    --seed 42 \
    --max_gap $max_gap \
    --yield_size 10000000
