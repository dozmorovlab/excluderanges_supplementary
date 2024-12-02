- [hg38_BL-reproducibility](hg38_BL-reproducibility): A set of Blacklists made using [01_generate-excludable-set-BL.sh](../../Blacklist/hg38_BL-reproducibility). Each set was made from 100 randomly sampled BAM files from this [table](../../data/ENCODE-IDS/hg38.input_bams.txt), kmer parameter = 36.

- [hg38_BL-various-sizes](hg38_BL-various-sizes): A set of Blacklists made using [01_generate-excludable-set-BL.sh](../../Blacklist/hg38_BL-reproducibility). Kundaje unified (KU), and Boyle's Blacklist (BL) are duplicated. All other sets have names <version>.<nBams>.bed. Each set is made using `nBams` accesion merged BAM files.
    - currenly making bed files from merged BAMs.
        - [x] 10 nBams
        - [x] 50 nBams
        - [x] 100 nBams
        - [x] 200 nBams
        - [x] 250 nBams (all identical, checked md5sums. This is because there 250 merged bam files, so each of these is made using the same bams.)

- [hg38_gaps](hg38_gaps): UCSC gaps for hg38, centromeres, telomeres, and shortarms. See `download_data()` in [utils.R](../../utils.R) for download details.

- [annotations](annotations): A set of Blacklists to compare original annotations with corrected annotations
    - [250.bed](annotations/250.bed): Generated using all 250 `nBams` and original Blacklist program from GitHub
    - [250e.bed](annotations/250e.bed): Generated using all 250 `nBams` and our ExpandedBlacklist program

- [graphs_original_v_expanded](graphs_original_v_expanded): A set of Blacklists to compare original and corrected distributions
    - [250.bed](graphs_original_v_expanded/250.bed): Generated using all 250 `nBams` and original Blacklist program from GitHub
    - [250e.bed](graphs_original_v_expanded/250e.bed): Generated using all 250 `nBams` and our ExpandedBlacklist program

- [graphs_ours_v_github](graphs_ours_v_github): A set of Blacklists to compare original and corrected distributions
    - [250.bed](graphs_ours_v_github/250.bed): Generated using all 250 `nBams` and original Blacklist program from GitHub
    - [hg38-blacklistv2.bed](graphs_ours_v_github/hg38-blacklistv2.bed): Original Blacklist from GitHub

- [tables](tables): A set of Blacklists to compare original GitHub blacklists with our generated blacklists
    - [250.bed](tables/250.bed): Generated using all 250 `nBams` and original Blacklist program from GitHub
    - [hg38-blacklistv2.bed](tables/hg38-blacklistv2.bed): Original Blacklist from GitHub

- [hg38](hg38): A combination of hg38 original, expanded, and github
    - [250.bed](hg38/250.bed): Generated using all 250 `nBams` and original Blacklist program from GitHub
    - [250e.bed](hg38/250e.bed): Generated using all 250 `nBams` and our ExpandedBlacklist program
    - [hg38-blacklistv2.bed](hg38/hg38-blacklistv2.bed): Original Blacklist from GitHub

- [mm10](mm10): A combination of mm10 original, expanded, and github
    - [266_mm10og.bed](mm10/266_mm10og.bed): Generated using all 266 `nBams` and original Blacklist program from GitHub
    - [266_mm10e.bed](mm10/266_mm10e.bed): Generated using all 266 `nBams` and our ExpandedBlacklist program
    - [mm10-blacklistv2.bed](mm10/mm10-blacklistv2.bed): Original Blacklist from GitHub