#!/usr/bin/env bash

NXF_SINGULARITY_CACHEDIR="./work/singularity"

nextflow run nf-core/sarek -r 3.5.0 \
    -c nf_sarek.config \
    -profile singularity \
    --input "./samplesheet.csv" \
    --outdir "./nf_hg38_sponge" \
    --tools haplotypecaller \
    --skip_tools fastqc \
    --aligner bwa-mem2 \
    --igenomes_ignore \
    --no_intervals \
    --download_cache \
    --fasta ./hg38_with_sponge.fasta \
    --known_snps ./Homo_sapiens_assembly38.dbsnp138.vcf.gz \
    --known_indels ./Homo_sapiens_assembly38.known_indels.vcf.gz \
    -resume

# Known Sites Sources
# https://storage.cloud.google.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf
# https://storage.cloud.google.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.known_indels.vcf.gz
# https://gatk.broadinstitute.org/hc/en-us/community/posts/360075305092-Known-Sites-for-BQSR