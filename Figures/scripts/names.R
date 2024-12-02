names_dict <- c(

  # Blacklists
  "hg38-blacklistv2" = "GitHub Blacklist",
  "250" = "Generated Blacklist",
  "250e" = "hg38 Our Expanded List",
  "KU" = "Kundaje Unified",
  "hg38.Nordin.CandRblacklist_hg38" = "Nordin CUT&RUN",

  # HS, LM, Both
  "101_MACS3_fc0.99_local_cent.merged.filtered" = "High Signal",
  "k100_b1000_f1000_0.01.multi.LM" = "Low Mappability",
  "HS_LM_CM" = "HS + LM",

  # GreyListChIP
  "GreyListChIP STAR 1k" = "GreyListChIP",

  # mm10
  "mm10-blacklistv2" = "mm10 GitHub Blacklist",
  "266_mm10og" = "mm10 Generated Blacklist",

  # Aligners
  "s1000_v100_b20000_k36_wt0.99_st0.999_36_bowtie2" = "Bowtie2 36bp",
  "s1000_v100_b20000_k36_wt0.99_st0.999_36_bwamem2" = "Bwa-mem2 36bp",
  "s1000_v100_b20000_k36_wt0.99_st0.999_36_starAligned.sortedByCoord.out" = "STAR 36bp",
  "s1000_v100_b20000_k36_wt0.99_st0.999_36.bam.sorted" = "Bwa-samse 36bp",

  "s1000_v100_b20000_k36_wt0.99_st0.999_101_bowtie2" = "Bowtie2 101bp",
  "s1000_v100_b20000_k36_wt0.99_st0.999_101_bwamem2" = "Bwa-mem2 101bp",
  "s1000_v100_b20000_k36_wt0.99_st0.999_101_starAligned.sortedByCoord.out" = "STAR 101bp",
  "s1000_v100_b20000_k36_wt0.99_st0.999_101.bam.sorted" = "Bwa-sampe 101bp",
  "s1000_v100_b20000_k36_wt0.99_st0.999_101_merged" = "Bwa-sampe Merged 101bp",
  #"s1000_v100_b20000_k101_wt0.99_st0.999_101_bowtie2" = "Bowtie2 101-BAMs k101 Blacklist-defaults",
  #"s1000_v100_b20000_k101_wt0.99_st0.999_101_bwamem2" = "Bwa-mem2 101-BAMs k101 Blacklist-defaults",
  #"s1000_v100_b20000_k101_wt0.99_st0.999_101_starAligned.sortedByCoord.out" = "STAR 101-BAMs k101 Blacklist-defaults",
  #"s1000_v100_b20000_k101_wt0.99_st0.999_101.bam.sorted" = "Bwa-sampe 101-BAMs k101 Blacklist-defaults",

  # Parameters
  "s1000_v100_b1000_k36_wt0.99_st0.999_n10_36_starAligned.sortedByCoord.out" = "STAR 36-BAMs b1000 k36 n10",
  "s1000_v100_b1000_k36_wt0.99_st0.999_n50_36_starAligned.sortedByCoord.out" = "STAR 36-BAMs b1000 k36 n50",
  "s1000_v100_b1000_k36_wt0.99_st0.999_n100_36_starAligned.sortedByCoord.out" = "STAR 36-BAMs b1000 k36 n100",
  "s1000_v100_b1000_k36_wt0.99_st0.999_n200_36_starAligned.sortedByCoord.out" = "STAR 36-BAMs b1000 k36 n200",
  "s1000_v100_b1000_k36_wt0.99_st0.999_n300_36_starAligned.sortedByCoord.out" = "STAR 36-BAMs b1000 k36 n300",
  "s1000_v100_b1000_k50_wt0.99_st0.999_n10_36_starAligned.sortedByCoord.out" = "STAR 36-BAMs b1000 k50 n10",
  "s1000_v100_b1000_k50_wt0.99_st0.999_n50_36_starAligned.sortedByCoord.out" = "STAR 36-BAMs b1000 k50 n50",
  "s1000_v100_b1000_k50_wt0.99_st0.999_n100_36_starAligned.sortedByCoord.out" = "STAR 36-BAMs b1000 k50 n100",
  "s1000_v100_b1000_k50_wt0.99_st0.999_n200_36_starAligned.sortedByCoord.out" = "STAR 36-BAMs b1000 k50 n200",
  "s1000_v100_b1000_k50_wt0.99_st0.999_n300_36_starAligned.sortedByCoord.out" = "STAR 36-BAMs b1000 k50 n300",
  "s1000_v100_b1000_k100_wt0.99_st0.999_n10_36_starAligned.sortedByCoord.out" = "STAR 36-BAMs b1000 k100 n10",
  "s1000_v100_b1000_k100_wt0.99_st0.999_n50_36_starAligned.sortedByCoord.out" = "STAR 36-BAMs b1000 k100 n50",
  "s1000_v100_b1000_k100_wt0.99_st0.999_n100_36_starAligned.sortedByCoord.out" = "STAR 36-BAMs b1000 k100 n100",
  "s1000_v100_b1000_k100_wt0.99_st0.999_n200_36_starAligned.sortedByCoord.out" = "STAR 36-BAMs b1000 k100 n200",
  "s1000_v100_b1000_k100_wt0.99_st0.999_n300_36_starAligned.sortedByCoord.out" = "STAR 36-BAMs b1000 k100 n300",
  "s1000_v100_b10000_k36_wt0.99_st0.999_n10_36_starAligned.sortedByCoord.out" = "STAR 36-BAMs b10000 k36 n10",
  "s1000_v100_b10000_k36_wt0.99_st0.999_n50_36_starAligned.sortedByCoord.out" = "STAR 36-BAMs b10000 k36 n50",
  "s1000_v100_b10000_k36_wt0.99_st0.999_n100_36_starAligned.sortedByCoord.out" = "STAR 36-BAMs b10000 k36 n100",
  "s1000_v100_b10000_k36_wt0.99_st0.999_n200_36_starAligned.sortedByCoord.out" = "STAR 36-BAMs b10000 k36 n200",
  "s1000_v100_b10000_k36_wt0.99_st0.999_n300_36_starAligned.sortedByCoord.out" = "STAR 36-BAMs b10000 k36 n300",
  "s1000_v100_b10000_k50_wt0.99_st0.999_n10_36_starAligned.sortedByCoord.out" = "STAR 36-BAMs b10000 k50 n10",
  "s1000_v100_b10000_k50_wt0.99_st0.999_n50_36_starAligned.sortedByCoord.out" = "STAR 36-BAMs b10000 k50 n50",
  "s1000_v100_b10000_k50_wt0.99_st0.999_n100_36_starAligned.sortedByCoord.out" = "STAR 36-BAMs b10000 k50 n100",
  "s1000_v100_b10000_k50_wt0.99_st0.999_n200_36_starAligned.sortedByCoord.out" = "STAR 36-BAMs b10000 k50 n200",
  "s1000_v100_b10000_k50_wt0.99_st0.999_n300_36_starAligned.sortedByCoord.out" = "STAR 36-BAMs b10000 k50 n300",
  "s1000_v100_b10000_k100_wt0.99_st0.999_n10_36_starAligned.sortedByCoord.out" = "STAR 36-BAMs b10000 k100 n10",
  "s1000_v100_b10000_k100_wt0.99_st0.999_n50_36_starAligned.sortedByCoord.out" = "STAR 36-BAMs b10000 k100 n50",
  "s1000_v100_b10000_k100_wt0.99_st0.999_n100_36_starAligned.sortedByCoord.out" = "STAR 36-BAMs b10000 k100 n100",
  "s1000_v100_b10000_k100_wt0.99_st0.999_n200_36_starAligned.sortedByCoord.out" = "STAR 36-BAMs b10000 k100 n200",
  "s1000_v100_b10000_k100_wt0.99_st0.999_n300_36_starAligned.sortedByCoord.out" = "STAR 36-BAMs b10000 k100 n300",
  "s1000_v100_b20000_k36_wt0.99_st0.999_n10_36_starAligned.sortedByCoord.out" = "STAR 36-BAMs b20000 k36 n10",
  "s1000_v100_b20000_k36_wt0.99_st0.999_n50_36_starAligned.sortedByCoord.out" = "STAR 36-BAMs b20000 k36 n50",
  "s1000_v100_b20000_k36_wt0.99_st0.999_n100_36_starAligned.sortedByCoord.out" = "STAR 36-BAMs b20000 k36 n100",
  "s1000_v100_b20000_k36_wt0.99_st0.999_n200_36_starAligned.sortedByCoord.out" = "STAR 36-BAMs b20000 k36 n200",
  "s1000_v100_b20000_k36_wt0.99_st0.999_n300_36_starAligned.sortedByCoord.out" = "STAR 36-BAMs b20000 k36 n300",
  "s1000_v100_b20000_k50_wt0.99_st0.999_n10_36_starAligned.sortedByCoord.out" = "STAR 36-BAMs b20000 k50 n10",
  "s1000_v100_b20000_k50_wt0.99_st0.999_n50_36_starAligned.sortedByCoord.out" = "STAR 36-BAMs b20000 k50 n50",
  "s1000_v100_b20000_k50_wt0.99_st0.999_n100_36_starAligned.sortedByCoord.out" = "STAR 36-BAMs b20000 k50 n100",
  "s1000_v100_b20000_k50_wt0.99_st0.999_n200_36_starAligned.sortedByCoord.out" = "STAR 36-BAMs b20000 k50 n200",
  "s1000_v100_b20000_k50_wt0.99_st0.999_n300_36_starAligned.sortedByCoord.out" = "STAR 36-BAMs b20000 k50 n300",
  "s1000_v100_b20000_k100_wt0.99_st0.999_n10_36_starAligned.sortedByCoord.out" = "STAR 36-BAMs b20000 k100 n10",
  "s1000_v100_b20000_k100_wt0.99_st0.999_n50_36_starAligned.sortedByCoord.out" = "STAR 36-BAMs b20000 k100 n50",
  "s1000_v100_b20000_k100_wt0.99_st0.999_n100_36_starAligned.sortedByCoord.out" = "STAR 36-BAMs b20000 k100 n100",
  "s1000_v100_b20000_k100_wt0.99_st0.999_n200_36_starAligned.sortedByCoord.out" = "STAR 36-BAMs b20000 k100 n200",
  "s1000_v100_b20000_k100_wt0.99_st0.999_n300_36_starAligned.sortedByCoord.out" = "STAR 36-BAMs b20000 k100 n300"

)