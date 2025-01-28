#!/bin/bash

# conda
source ~/miniconda3/etc/profile.d/conda.sh
conda activate deeptools

# Input
bams="./bams"
bams_sponge="./bams_sponge"
ex_list_dir="./ex_lists"

# Output directories
out_dir="./summary"
out_sponge_dir="./summary_sponge"
k=0

# Threads
threads=64

mkdir -p "${out_dir}"
mkdir -p "${out_sponge_dir}"

for bam in `find "$bams" -type f -name "*_Aligned.sortedByCoord.out.bam"`; do
    
    echo "$(date '+%Y-%m-%d %H:%M:%S'): Indexing: $bam"
    samtools index "$bam"

    neutral_base=$(basename "$bam" _Aligned.sortedByCoord.out.bam)

    echo "$(date '+%Y-%m-%d %H:%M:%S'): Getting summary: $bam"

    multiBamSummary bins \
        --bamfiles "$bam" \
        --outFileName "${out_dir}/${neutral_base}.no_ex.npz" \
        --numberOfProcessors max \
        --outRawCounts "${out_dir}/${neutral_base}.no_ex.counts.tsv"

    for ex_list in "$ex_list_dir"/*; do
        
        ex_base="$(basename "$ex_list" .bed)"

        echo "$(date '+%Y-%m-%d %H:%M:%S'): Getting summary: $bam + $ex_base"
        
        multiBamSummary bins \
            --bamfiles "$bam" \
            --outFileName "${out_dir}/${neutral_base}.${ex_base}.npz" \
            --numberOfProcessors max \
            --outRawCounts "${out_dir}/${neutral_base}.${ex_base}.counts.tsv" \
            --blackListFileName "$ex_list"
    done

done

for bam in `find "$bams_sponge" -type f -name "*_Aligned.sortedByCoord.out.bam"`; do
    
    echo "$(date '+%Y-%m-%d %H:%M:%S'): Indexing: $bam"
    samtools index "$bam"

    neutral_base=$(basename "$bam" _Aligned.sortedByCoord.out.bam)

    echo "$(date '+%Y-%m-%d %H:%M:%S'): Getting summary: $bam"

    multiBamSummary bins \
        --bamfiles "$bam" \
        --outFileName "${out_sponge_dir}/${neutral_base}.no_ex.npz" \
        --numberOfProcessors max \
        --outRawCounts "${out_sponge_dir}/${neutral_base}.no_ex.counts.tsv"

    for ex_list in "$ex_list_dir"/*; do
        
        ex_base="$(basename "$ex_list" .bed)"

        echo "$(date '+%Y-%m-%d %H:%M:%S'): Getting summary: $bam + $ex_base"

        multiBamSummary bins \
            --bamfiles "$bam" \
            --outFileName "${out_sponge_dir}/${neutral_base}.${ex_base}.npz" \
            --numberOfProcessors max \
            --outRawCounts "${out_sponge_dir}/${neutral_base}.${ex_base}.counts.tsv" \
            --blackListFileName "$ex_list"
    done

done

readarray -t bams_array < <(find "$bams" -type f -name "*_Aligned.sortedByCoord.out.bam")
readarray -t bams_sponge_array < <(find "$bams_sponge" -type f -name "*_Aligned.sortedByCoord.out.bam")

echo "$(date '+%Y-%m-%d %H:%M:%S'): Getting summary: all"

multiBamSummary bins \
    --bamfiles "${bams_array[@]}" \
    --outFileName "${out_dir}/all.npz" \
    --numberOfProcessors max \
    --outRawCounts "${out_dir}/all.counts.tsv"

echo "$(date '+%Y-%m-%d %H:%M:%S'): Getting summary: all (sponge)"

multiBamSummary bins \
    --bamfiles "${bams_sponge_array[@]}" \
    --outFileName "${out_sponge_dir}/all.npz" \
    --numberOfProcessors max \
    --outRawCounts "${out_sponge_dir}/all.counts.tsv"

for ex_list in "$ex_list_dir"/*; do
    
    ex_base="$(basename "$ex_list" .bed)"

    echo "$(date '+%Y-%m-%d %H:%M:%S'): Getting summary: all + $ex_base"

    multiBamSummary bins \
        --bamfiles "${bams_array[@]}" \
        --outFileName "${out_dir}/all.${ex_base}.npz" \
        --numberOfProcessors max \
        --outRawCounts "${out_dir}/all.${ex_base}.counts.tsv" \
        --blackListFileName "$ex_list"

    echo "$(date '+%Y-%m-%d %H:%M:%S'): Getting summary: all (sponge) + $ex_base"

    multiBamSummary bins \
        --bamfiles "${bams_sponge_array[@]}" \
        --outFileName "${out_sponge_dir}/all.${ex_list}.npz" \
        --numberOfProcessors max \
        --outRawCounts "${out_sponge_dir}/all.${ex_list}.counts.tsv" \
        --blackListFileName "$ex_list"
done

for ex_list in "$ex_list_dir"/*; do
    ex_base="$(basename "$ex_list" .bed)"

    echo "$(date '+%Y-%m-%d %H:%M:%S'): Getting summary: all + $ex_base"

    multiBamSummary BED-file \
        --bamfiles "${bams_array[@]}" \
        --outFileName "${out_dir}/all.${ex_base}.ex.npz" \
        --numberOfProcessors max \
        --outRawCounts "${out_dir}/all.${ex_base}.counts.ex.tsv" \
        --BED "$ex_list"

    echo "$(date '+%Y-%m-%d %H:%M:%S'): Getting summary: all (sponge) + $ex_base"

    multiBamSummary BED-file \
        --bamfiles "${bams_sponge_array[@]}" \
        --outFileName "${out_sponge_dir}/all.${ex_base}.ex.npz" \
        --numberOfProcessors max \
        --outRawCounts "${out_sponge_dir}/all.${ex_base}.counts.ex.tsv" \
        --BED "$ex_list"
done