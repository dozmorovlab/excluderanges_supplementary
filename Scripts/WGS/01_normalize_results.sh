#!/usr/bin/env bash

threads=24

# benchmark
bcftools norm ../HG001_GRCh38_1_22_v4.2.1_benchmark.vcf.gz \
    --atomize \
    --rm-dup all \
    --regions chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22 \
    --threads $threads \
    --write-index tbi \
    --output ./data/benchmark.norm.vcf.gz

# hg38
#bcftools norm nf_hg38/variant_calling/haplotypecaller/S1/S1.haplotypecaller.filtered.vcf.gz \
bcftools norm ../nf_hg38/variant_calling/haplotypecaller/S1/S1.haplotypecaller.filtered.vcf.gz \
    --atomize \
    --rm-dup all \
    --regions chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22 \
    --threads $threads \
    --write-index \
    --output ./data/hg38.norm.vcf.gz

# hg38 + sponge
bcftools norm ../nf_hg38_sponge/variant_calling/haplotypecaller/S1/S1.haplotypecaller.filtered.vcf.gz \
    --atomize \
    --rm-dup all \
    --regions chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22 \
    --threads $threads \
    --write-index \
    --output ./data/sponge.norm.vcf.gz
