# Beyond Blacklists: A Critical Assessment of Exclusion Set Generation Strategies and Alternative Approaches

- [Preprint](#preprint)
- [Abstract](#abstract)
- [Figures](#figures)
  - [Figures](#figures-1)
  - [Tables](#tables)
  - [The Blacklist algorithm](#the-blacklist-algorithm)
  - [hg38 sponge sequences](#hg38-sponge-sequences)
  - [Scripts](#scripts)

<!-- README.md is generated from README.Rmd. Please edit this file -->
<!-- rmarkdown::render("README.Rmd", output_format = "github_document") -->
<!-- Manually move the first header above the TOC -->

![Lifecycle:
stable](https://img.shields.io/badge/lifecycle-stable-brightgreen.svg)
![License:
GPL-3.0](https://img.shields.io/badge/license-GPL%203-blue.svg)

## Preprint

## Abstract

Short-read sequencing data can be affected by alignment artifacts in
certain genomic regions. Removing reads overlapping these exclusion
regions, previously known as Blacklists, help to potentially improve
biological signal. Tools like the widely used Blacklist software
facilitate this process, but their algorithmic details and parameter
choices are not always clearly documented, affecting reproducibility and
biological relevance. We examined the Blacklist software and found that
pre-generated exclusion sets were difficult to reproduce due to
variability in input data, aligner choice, and read length. We also
identified and addressed a coding issue that led to over-annotation of
high-signal regions. We further explored the use of “sponge”
sequences—unassembled genomic regions such as satellite DNA, ribosomal
DNA, and mitochondrial DNA—as an alternative approach. Aligning reads to
a genome that includes sponge sequences reduced signal correlation in
ChIP-seq data comparably to Blacklist-derived exclusion sets while
preserving biological signal. Sponge-based alignment also had minimal
impact on RNA-seq gene counts, suggesting broader applicability beyond
chromatin profiling. These results highlight the limitations of fixed
exclusion sets and suggest that sponge sequences offer a flexible,
alignment-guided strategy for reducing artifacts and improving
functional genomics analyses.

## Figures

![Figure 1](Figures/figures/Figure_1.svg) **Figure 1. Differences
between the GitHub version of the hg38 exclusion set, our hg38 exclusion
set generated with the Blacklist software, and the reference Kundaje
Unified set.** (A) Differences in count, coverage, and width
distribution. (B) Differences in gap coverage. (C) Pairwise overlaps of
region counts and genome coverage between sets.

![Figure 2](Figures/figures/Figure_2.svg) **Figure 2. Differences
between the GitHub version of the hg38 exclusion set and manually
defined gold-standard exclusion sets.** (A) Differences in count. (B)
Pairwise overlaps of region counts between the GitHub exclusion set and
manually defined gold-standard exclusion sets. (C) Differences in width
distribution. (D) Multidimensional scaling plot of Jaccard count
similarity among the GitHub exclusion set and manually defined gold
standards.

![Figure 3](Figures/figures/Figure_3.svg) **Figure 3. Differences
between the hg38 exclusion sets, the reference Kundaje Unified set, and
sets generated from 36bp single-end reads realigned with different
aligners.** A) Count, coverage, and width distribution differences; B)
Gap coverage differences; C) Pairwise overlaps of counts and coverage
between sets.

![Figure 4](Figures/figures/Figure_4.svg) **Figure 4. Biological
characterization of genes affected by exclusion regions.** A) Count and
coverage of protein-coding, long noncoding, and other transcripts
affected by exclusion sets; B) Representative comparison of exclusion
set coverage over a cluster of ubiquitin-specific peptidase 17-like
family member genes; C) KEGG pathways enriched in genes overlapped by
exclusion sets.

![Figure 5](Figures/figures/Figure_5.svg) **Figure 5. Changes in
ChIP-seq signal correlation with and without reads overlapping exclusion
sets or aligned to the sponge.** Data for the Gm12878 cell line is
shown. A) Heatmaps of correlation differences for each exclusion set,
sorted by mean (red/blue gradient corresponds to decreases/increases in
correlations, respectively); B) Correlation difference distributions for
the top three most affected transcription factors.

### [Figures](Figures)

Supplementary figures and scripts needed to make the figures, and the
figures themselves.

### [Tables](Tables)

**Table 1. Characteristics of the hg38 exclusion sets.** Total number of
regions, average region width, number of regions overlapping gaps, and
proportion of gap coverage.
<img src="README_files/figure-gfm/table function-1.png" width="1219" /><img src="README_files/figure-gfm/table function-2.png" width="1219" />

### The [Blacklist](https://github.com/Boyle-Lab/Blacklist) algorithm

- [Supplementary Note](Figures/Supplementary_Note.pdf)
- Blacklist [original and corrected C++
  code](Scripts/Blacklist/Blacklist_corrected)
- Blacklist [Python
  reimplementation](Scripts/Blacklist/Blacklist_python)

![Figure S5](Figures/figures/Supplementary_Figure_S5.svg)
**Supplementary Figure S5. The Blacklist Algorithm.** A) **Binning
parameters:** Two internal parameters, “binSize” and “binOverlap,”
define how bins are created. B) **Read sorting:** The start position of
each read is sorted according to the read’s length and compared using
the given pre-calculated mappability vector. C) **The “binsMap”
vector:** This vector is generated using the “uniqueLength” threshold.
D) **Normalization:** The signal (“Unique”) and mappability (“Multi”)
vectors are divided, quantile normalized, and the median is taken per
bin. E) **Region calling:** This step uses several predefined internal
parameters/thresholds to call and annotate excludable regions.

### [hg38 sponge sequences](Sponge)

[GRCh38_sponge.fa.gz](Sponge/GRCh38_sponge.fa.gz), 25Mb

**Source:** Miga, Karen H., Christopher Eisenhart, and W. James Kent.
“Utilizing mapping targets of sequences underrepresented in the
reference assembly to reduce false positive alignments.” Nucleic acids
research (2015), <https://doi.org/10.1093/nar/gkv671>, [Supplementary
Data](https://academic.oup.com/nar/article/43/20/e133/1398360#81910257)

### [Scripts](Scripts)

A directory containing various scripts used to generate the data.
