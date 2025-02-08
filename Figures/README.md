The data and scripts needed to make the figures, and the figures themselves used in our publication.

# Supplementary Figures

![Figure S1](Figures/figures/Supplementary_Figure_S1.svg)
**Supplementary Figure S1. Differences between the GitHub version of the mm10 exclusion set and our mm10 exclusion set generated with the Blacklist software.** (A) Differences in count, coverage, and width distribution. (B) Differences in gap coverage. (C) Pairwise overlaps of region counts and genome coverage between sets.

![Figure S2](Figures/figures/Supplementary_Figure_S2.svg)
**Supplementary Figure S2. The diversity of BAM files used to generate the hg38 exclusion set.** (A) Read length proportions among 1,225 BAM files. (B) Number of mapped reads for single- and paired-end files. (C) Number of single- and paired-end files per donor. (D) Total number of mapped and unmapped reads per donor.

![Figure S3](Figures/figures/Supplementary_Figure_S3.svg)
**Supplementary Figure S3. Exclusion set characteristics generated from the 274 unmerged 101bp BAM files and the corresponding merged BAM file.** (A) Differences in count, coverage, and width distribution. (B) Differences in gap coverage. (C) Pairwise overlaps of region counts and genome coverage between sets.

![Figure S4](Figures/figures/Supplementary_Figure_S4.svg)
**Supplementary Figure S4. Differences between the GitHub version of the hg38 exclusion set and manually defined gold-standard exclusion sets.** (A) Differences in width. (B) Pairwise width overlaps between the GitHub exclusion set and manually defined gold-standard exclusion sets. (C) Differences in gap coverage. (D) Multidimensional scaling plot of Forbes width similarity among the GitHub exclusion set and manually defined gold standards.

![Figure S5](Figures/figures/Supplementary_Figure_S5.svg)
**Supplementary Figure S5. The Blacklist Algorithm.**  A) **Binning parameters:** Two internal parameters, "binSize" and "binOverlap," define how bins are created. B) **Read sorting:** The start position of each read is sorted according to the read’s length and compared using the given pre-calculated mappability vector. C) **The “binsMap” vector:** This vector is generated using the "uniqueLength" threshold. D) **Normalization:** The signal ("Unique") and mappability ("Multi") vectors are divided, quantile normalized, and the median is taken per bin. E) **Region calling:** This step uses several predefined internal parameters/thresholds to call and annotate excludable regions.  

![Figure S6](Figures/figures/Supplementary_Figure_S6.svg)
**Supplementary Figure S6. Relationships among realigned 36bp exclusion sets.** A, B) Jaccard count overlap, multidimensional scaling, and hierarchical clustering; C, D) Forbes width overlap, MDS, and clustering.

![Figure S7](Figures/figures/Supplementary_Figure_S7.svg)
**Supplementary Figure S7. Relationships among realigned 36bp and 101bp exclusion sets.** A) Number and width distribution; B) Jaccard count overlap and Forbes width overlap MDS plots.

![Figure S8](Figures/figures/Supplementary_Figure_S8.svg)
**Supplementary Figure S8. The effect of varying the number of BAM files, the "bridge" parameter, and the "k-mer" parameter on the output of the Blacklist software.** The effect on the number of exclusion regions, sum of region widths, and mean and median region widths.

![Figure S9](Figures/figures/Supplementary_Figure_S9.svg)
**Supplementary Figure S9. Jaccard count overlap (A-D) and Forbes similarity (E-H) MDS plots colored by parameters.** A), E) Number of BAM files; B), F) Bridge parameter; C), G) k-mer parameter; D), H) Parameter sweep results (gray) vs. reference exclusion sets.

![Figure S10](Figures/figures/Supplementary_Figure_S10.svg)
**Supplementary Figure S10. Relationships among Blacklist, manually created, and Nordin CUT&RUN and GreyListChIP exclusion sets.** A, B) Jaccard count overlap, multi-dimensional scaling, and hierarchical clustering; C, D) Forbes width overlap, MDS, and clustering.

![Figure S11](Figures/figures/Supplementary_Figure_S11.svg)
**Supplementary Figure S11. Transcription factor-specific changes in ChIP-seq signal correlation with and without reads overlapping exclusion sets or aligned to the sponge.** Data for the Gm12878 cell line is shown. Correlation difference distributions, red transcription factor names highlights cases where exclusion lists other than "Blacklist" and "Sponge" have the strongest decrease in correlations.

![Figure S12](Figures/figures/Supplementary_Figure_S12.svg)
**Supplementary Figure S12. The effect of "sponge" sequences on alignment.** A) The number and proportion of 36 bp and 101 bp reads aligned to "sponge" sequences from the whole genome and from various excludable sets; B) The total number of reads overlapping annotated gene regions when aligned with and without the "sponge" sequences.
