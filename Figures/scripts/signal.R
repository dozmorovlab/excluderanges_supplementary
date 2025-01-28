library(readr)
library(ggplot2)
library(tidyr)
library(dplyr)
library(patchwork)

setwd(file.path(here::here(), "Figures"))

results_dir <- file.path("results")

bp <- 101

get_sig_plot <- function(bp) {

  data_folder <- file.path("data", paste0("summary_", bp))
  og_bam <- paste0("'star_", bp, ".bam'")
  sponge_bam <- paste0("'star_", bp, "_sponge.bam'")
  #sponge_bam <- paste0("'star_", bp, "_sponge.sorted.bam'")

  names_dict <- setNames(
    c(
      "Whole Genome",
      "GitHub Blacklist",
      "Generated Blacklist",
      "Kundaje Unified",
      "High Signal",
      "Low Mappability",
      "HS + LM",
      "Nordin CUT&RUN",
      "GreyListChIP"
    ),
    c(
      paste0("data/summary_", bp, "/", bp, ".counts.tsv"),
      paste0("data/summary_", bp, "/", bp, ".hg38-blacklistv2.counts.ex.tsv"),
      paste0("data/summary_", bp, "/", bp, ".250.counts.ex.tsv"),
      paste0("data/summary_", bp, "/", bp, ".KU.counts.ex.tsv"),
      paste0("data/summary_", bp, "/", bp, ".101_MACS3_fc0.99_local_cent.merged.filtered.counts.ex.tsv"),
      paste0("data/summary_", bp, "/", bp, ".k100_b1000_f1000_0.01.multi.LM.counts.ex.tsv"),
      paste0("data/summary_", bp, "/", bp, ".HS_LM_CM.counts.ex.tsv"),
      paste0("data/summary_", bp, "/", bp, ".hg38.Nordin.CandRblacklist_hg38.counts.ex.tsv"),
      paste0("data/summary_", bp, "/", bp, ".GreyListChIP STAR 1k.counts.ex.tsv")
    )
  )

  make_datatable <- function(df) {
    return(
      DT::datatable(
        df,
        options = list(pageLength = nrow(df)),
        rownames = FALSE
      )
    )
  }

  categorize_chr <- function(chr) {
    # Check categories in order of specificity
    if (grepl("_random$", chr)) {
      return("Random")
    } else if (grepl("^chrUn", chr)) {
      return("Unplaced")
    } else if (grepl("^chr[0-9]+$", chr)) {
      return("Autosome")
    } else if (grepl("^chr[XY]$", chr)) {
      return("Sex Chromosome")
    } else if (grepl("^chrM$", chr)) {
      return("Mitochondrial")
    } else if (grepl("^chrEBV$", chr)) {
      return("Viral (EBV)")
    } else {
      return("Unknown")
    }
  }

  # Get a list of all .tsv files in the directory
  tsv_files <- list.files(
    path = data_folder,
    pattern = "\\.tsv$",
    full.names = TRUE
  )

  # Read all .tsv files into a list of data frames
  data_list <- sapply(
    tsv_files,
    read_tsv,
    simplify = FALSE
  )

  updated_data_list <- lapply(names(data_list), function(file_name) {
    df <- data_list[[file_name]]
    df$chr_type <- sapply(df$`#'chr'`, categorize_chr)
    return(df)
  })
  names(updated_data_list) <- names(data_list)

  sponge_category_counts <- as.data.frame(sapply(names(data_list), function(file_name) {
    df <- updated_data_list[[file_name]]
    sapply(
      c(
        "Autosome",
        "Sex Chromosome",
        "Mitochondrial",
        "Viral (EBV)",
        "Unplaced",
        "Random",
        "Unknown"
      ),
      function(type) {
        sum(df[df$chr_type == type, ][sponge_bam])
      }
    )
  }))

  sponge_count <- sapply(names(data_list), function(file_name) {
    df <- data_list[[file_name]]
    out_count <- sum(df[og_bam]) - sum(df[sponge_bam])
    return(out_count)
  })
  hg38_count <- sapply(names(data_list), function(file_name) {
    df <- data_list[[file_name]]
    out_count <- sum(df[sponge_bam])
    return(out_count)
  })
  total_count <- sapply(names(data_list), function(file_name) {
    df <- data_list[[file_name]]
    out_count <- sum(df[og_bam])
    return(out_count)
  })
  summary_df <- data.frame(
    file = names_dict[names(sponge_count)],
    `Percent Sponge` = sponge_count / total_count,
    `Percent hg38` = 1 - (sponge_count / total_count),
    `hg38 Count` = hg38_count,
    `Sponge Count` = sponge_count,
    `Autosome` = t(sponge_category_counts["Autosome", ]),
    `Percent Autosome` = t(sponge_category_counts["Autosome", ]) / total_count,
    `Sex Chromosome` = t(sponge_category_counts["Sex Chromosome", ]),
    `Percent Sex Chromosome` = t(sponge_category_counts["Sex Chromosome", ]) / total_count,
    `Mitochondrial` = t(sponge_category_counts["Mitochondrial", ]),
    `Percent Mitochondrial` = t(sponge_category_counts["Mitochondrial", ]) / total_count,
    `Viral (EBV)` = t(sponge_category_counts["Viral (EBV)", ]),
    `Percent Viral (EBV)` = t(sponge_category_counts["Viral (EBV)", ]) / total_count,
    `Unplaced` = t(sponge_category_counts["Unplaced", ]),
    `Percent Unplaced` = t(sponge_category_counts["Unplaced", ]) / total_count,
    `Random` = t(sponge_category_counts["Random", ]),
    `Percent Random` = t(sponge_category_counts["Random", ]) / total_count,
    `Unknown` = t(sponge_category_counts["Unknown", ]),
    `Percent Unknown` = t(sponge_category_counts["Unknown", ]) / total_count,
    check.names = FALSE,
    fix.empty.names = FALSE
  )
  colnames(summary_df) <- c(
    "file", "Percent Sponge", "Percent hg38", "hg38 Count", "Sponge Count",
    "Autosome", "Percent Autosome", "Sex Chromosome", "Percent Sex Chromosome",
    "Mitochondrial", "Percent Mitochondrial", "Viral (EBV)", "Percent Viral (EBV)",
    "Unplaced", "Percent Unplaced", "Random", "Percent Random",
    "Unknown", "Percent Unknown"
  )

  summary_df$file <- factor(summary_df$file, levels = rev(names_dict))
  # ----------------------------
  read_counts_long_df <- summary_df %>%
    pivot_longer(
      cols = c("hg38 Count", "Sponge Count"),
      names_to = "type",
      values_to = "count"
    )
  read_counts_long_df$type <- factor(
    read_counts_long_df$type,
    levels = c("hg38 Count", "Sponge Count")
  )
  counts_plot <- ggplot(
    read_counts_long_df,
    aes(x = file, y = count, fill = type)
  ) +
    geom_bar(stat = "identity") +
    labs(
      x = NULL,
      y = "Read count",
      title = "Read counts",
      subtitle = paste0(bp, "bp BAM")
    ) +
    theme_minimal() +
    coord_flip() +
    scale_fill_manual(
      values = c("hg38 Count" = "grey", "Sponge Count" = "black")
    )

  read_per_long_df <- summary_df %>%
    pivot_longer(
      cols = c("Percent Sponge", "Percent hg38"),
      names_to = "type",
      values_to = "count"
    )
  read_per_long_df$type <- factor(
    read_per_long_df$type,
    levels = c("Percent hg38", "Percent Sponge")
  )
  percent_plot <- ggplot(
    read_per_long_df,
    aes(x = file, y = count, fill = type)
  ) +
    geom_bar(stat = "identity") +
    labs(
      x = NULL,
      y = NULL,
      title = "Read percents",
      subtitle = paste0(bp, "bp BAM")
    ) +
    theme_minimal() +
    coord_flip() +
    scale_fill_manual(
      values = c("Percent hg38" = "grey", "Percent Sponge" = "black")
    ) +
    scale_y_continuous(labels = function(x) paste0(x * 100, "%"))

  write.csv(summary_df, file.path(results_dir, paste0(bp, "_summary.csv")))

  return(list("counts" = counts_plot, "percents" = percent_plot))
}

plot_36_list <- get_sig_plot(36)
plot_101_list <- get_sig_plot(101)

c_plot <- plot_36_list$counts / plot_101_list$counts +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom", legend.title = element_blank()) &
  labs(title = NULL) &
  guides(fill = guide_legend(reverse = TRUE))

p_plot <- plot_36_list$percents / plot_101_list$percents +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom", legend.title = element_blank()) &
  labs(title = NULL) &
  guides(fill = guide_legend(reverse = TRUE))

signal_plot <- c_plot | p_plot
