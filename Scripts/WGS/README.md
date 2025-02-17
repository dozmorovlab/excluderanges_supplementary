# WGS

Scripts used to evaluate how our "sponge" affects WGS analysis.

[00_WGS.sh](00_WGS.sh): Uses nextflow core's [Sarek](https://nf-co.re/sarek/3.5.0) to call variants from WGS data.

[01_normalize_results.sh](01_normalize_results.sh): Normalizes the VCF files with and without our "Sponge" so that SNPs can be counted on autosomes.

[03_plot.R](03_plot.R): An `R` script for plotting Lost or Gained SNPs and INDELs.
