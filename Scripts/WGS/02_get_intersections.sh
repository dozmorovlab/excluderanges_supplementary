#!/usr/bin/env bash

#Using the following file names:
#0000.vcf	for records private to	A
#0001.vcf	for records private to	B
#0002.vcf	for records from A shared by both A & B
#0003.vcf	for records from B shared by both A & B

benchmark=./data/benchmark.norm.vcf.gz
hg38=./data/hg38.norm.vcf.gz
sponge=./data/sponge.norm.vcf.gz

bcftools isec "$hg38" "$sponge" -p ./data/isec

# Rename Files
mv ./data/isec/0000.vcf ./data/isec/lost.vcf
mv ./data/isec/0001.vcf ./data/isec/gained.vcf
mv ./data/isec/0002.vcf ./data/isec/hg38_shared.vcf
mv ./data/isec/0003.vcf ./data/isec/sponge_shared.vcf

cd ./data/isec
OUTPUT_FILE="./../isec_counts.csv"
echo "File,TOTAL,SNPS,INDELS" > $OUTPUT_FILE

# Count variants
for vcf in lost.vcf gained.vcf hg38_shared.vcf sponge_shared.vcf; do
    if [[ -f "$vcf" ]]; then
        TOTAL=$(bcftools view "$vcf" -H | wc -l)
        SNPS=$(bcftools view -v snps "$vcf" -H | wc -l)
        INDELS=$(bcftools view -v indels "$vcf" -H | wc -l)

        echo "$vcf,$TOTAL,$SNPS,$INDELS" >> $OUTPUT_FILE
    else
        echo "Warning: $vcf not found, skipping."
    fi
done

cd ../..

bgzip ./data/isec/lost.vcf
bgzip ./data/isec/gained.vcf

bcftools index ./data/isec/lost.vcf.gz
bcftools index ./data/isec/gained.vcf.gz

bcftools isec ./data/isec/lost.vcf.gz "$benchmark" -p ./data/lost
bcftools isec ./data/isec/gained.vcf.gz "$benchmark" -p ./data/gained

# Rename Files
mv ./data/lost/0000.vcf ./data/lost/lost_unique.vcf
mv ./data/lost/0001.vcf ./data/lost/benchmark_unique.vcf
mv ./data/lost/0002.vcf ./data/lost/lost_shared.vcf
mv ./data/lost/0003.vcf ./data/lost/benchmark_shared.vcf
mv ./data/gained/0000.vcf ./data/gained/gained_unique.vcf
mv ./data/gained/0001.vcf ./data/gained/benchmark_unique.vcf
mv ./data/gained/0002.vcf ./data/gained/gained_shared.vcf
mv ./data/gained/0003.vcf ./data/gained/benchmark_shared.vcf

cd ./data/lost
OUTPUT_FILE="./../lost_counts.csv"
echo "File,TOTAL,SNPS,INDELS" > $OUTPUT_FILE

for vcf in lost_unique.vcf lost_shared.vcf benchmark_shared.vcf; do
    if [[ -f "$vcf" ]]; then
        TOTAL=$(bcftools view "$vcf" -H | wc -l)
        SNPS=$(bcftools view -v snps "$vcf" -H | wc -l)
        INDELS=$(bcftools view -v indels "$vcf" -H | wc -l)

        echo "$vcf,$TOTAL,$SNPS,$INDELS" >> $OUTPUT_FILE
    else
        echo "Warning: $vcf not found, skipping."
    fi
done

cd ../..

cd ./data/gained
OUTPUT_FILE="./../gained_counts.csv"
echo "File,TOTAL,SNPS,INDELS" > $OUTPUT_FILE

for vcf in gained_unique.vcf gained_shared.vcf benchmark_shared.vcf; do
    if [[ -f "$vcf" ]]; then
        TOTAL=$(bcftools view "$vcf" -H | wc -l)
        SNPS=$(bcftools view -v snps "$vcf" -H | wc -l)
        INDELS=$(bcftools view -v indels "$vcf" -H | wc -l)

        echo "$vcf,$TOTAL,$SNPS,$INDELS" >> $OUTPUT_FILE
    else
        echo "Warning: $vcf not found, skipping."
    fi
done

cd ../..

lost.vcf,112623,97584,15039
gained.vcf,23793,19576,4217
hg38_shared.vcf,4574277,3774620,799657
sponge_shared.vcf,4574277,3774620,799657
