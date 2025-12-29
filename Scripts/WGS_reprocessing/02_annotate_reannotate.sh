#!bin/bash
# srun -c 126 --mem=494G -J annotate --pty bash
# cd /lustre/home/juicer/Blacklist/WGS_reprocessing

cpus=126

# bgzip hub_3671779_hgUniqueHg38.bed
# tabix -p bed hub_3671779_hgUniqueHg38.bed.gz

module load bcftools

results_dir="./annotated"
results_dir_reannotated="./reannotated"

mkdir -p "$results_dir_reannotated"

selected_chrs="-r chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22"

# BED file with regions to exclude
exclude_bed="hub_3671779_hgUniqueHg38.bed.gz"

# Extract ids ##############################################################

for type in snps indels; do

    echo "$(date) Extracting ${type} ids (excluding hgUnique regions)..."

    # T2T_CHM13
    in_vcf="${results_dir}/T2T_CHM13.annotated.vcf.gz"
    out_ids="${results_dir_reannotated}/ids/T2T_CHM13/autosomes.${type}.ids.txt"
    mkdir -p "$(dirname "$out_ids")"

    bcftools view \
        -v "$type" \
        $selected_chrs \
        -T "^${exclude_bed}" \
        --threads "$cpus" \
        --no-header \
        "$in_vcf" | \
    cut -f3 > \
    "$out_ids"

done
