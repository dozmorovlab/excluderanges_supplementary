#!/bin/bash

# conda
source ~/miniconda3/etc/profile.d/conda.sh
conda activate deeptools

# Input
bams="./bams_autosomes"
bams_sponge="./bams_sponge_autosomes"
ex_list_dir="./ex_lists"

# Output directories
out_dir="./summary_autosomes"

mkdir -p "${out_dir}"

readarray -t bams_array < <(find "$bams" -type f -name "*_Aligned.sortedByCoord.out.autosomes.bam")
readarray -t bams_sponge_array < <(find "$bams_sponge" -type f -name "*_Aligned.sortedByCoord.out.autosomes.bam")

#echo "$(date '+%Y-%m-%d %H:%M:%S'): Getting summary: all"

#multiBamSummary bins \
#    --bamfiles "${bams_array[@]}" \
#    --outFileName "${out_dir}/all.npz" \
#    --numberOfProcessors max \
#    --outRawCounts "${out_dir}/all.counts.tsv"

#echo "$(date '+%Y-%m-%d %H:%M:%S'): Getting summary: all (sponge)"

#multiBamSummary bins \
#    --bamfiles "${bams_sponge_array[@]}" \
#    --outFileName "${out_dir}/all.sponge.npz" \
#    --numberOfProcessors max \
#    --outRawCounts "${out_dir}/all.sponge.counts.tsv"

for ex_list in "$ex_list_dir"/*; do
    
    ex_base="$(basename "$ex_list" .bed)"

#    echo "$(date '+%Y-%m-%d %H:%M:%S'): Getting summary: all + $ex_base"

#    multiBamSummary bins \
#        --bamfiles "${bams_array[@]}" \
#        --outFileName "${out_dir}/all.${ex_base}.npz" \
#        --numberOfProcessors max \
#        --outRawCounts "${out_dir}/all.${ex_base}.counts.tsv" \
#        --blackListFileName "$ex_list"

    echo "$(date '+%Y-%m-%d %H:%M:%S'): Getting summary: all (sponge) + $ex_base"

    multiBamSummary bins \
        --bamfiles "${bams_sponge_array[@]}" \
        --outFileName "${out_dir}/all.${ex_base}.sponge.npz" \
        --numberOfProcessors max \
        --outRawCounts "${out_dir}/all.${ex_base}.sponge.counts.tsv" \
        --blackListFileName "$ex_list"
done
