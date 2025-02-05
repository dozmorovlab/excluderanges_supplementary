# Define the list of required packages
required_packages <- c(
  "limma",
  "data.table",
  "tidyverse",
  "edgeR",
  "TxDb.Hsapiens.UCSC.hg38.knownGene",
  "TxDb.Mmusculus.UCSC.mm10.knownGene",
  "GenomicRanges",
  "ggplot2",
  "ggpattern",
  "openxlsx",
  "IRanges",
  "dplyr",
  "ggvenn",
  "patchwork",
  "ggrepel",
  "ggdendro",
  "cowplot",
  "magick",
  "here",
  "rtracklayer",
  "dplyr",
  "ComplexHeatmap",
  "gridExtra",
  "gdata",
  "AnnotationHub",
  "readr",
  "ggsci",
  "tidyr",
  "biomaRt",
  "jsonlite",
  "reshape2",
  "viridis",
  "uwot",
  "ggrepel",
  "knitr",
  "pander",
  "readxl",
  "stringr",
  "ggridges",
  "svglite"
)

# Function to check if a package is installed and install it if not
install_if_missing <- function(package) {
  if (!requireNamespace(package, quietly = TRUE)) {
    BiocManager::install(package)
  }
}

# Ensure BiocManager is installed
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# Install missing packages
invisible(lapply(required_packages, install_if_missing))

# Load all packages
invisible(lapply(required_packages, library, character.only = TRUE))

setwd(file.path(here::here(), "Figures"))

# Get functions
source(file.path("scripts", "utils.R"))
source(file.path("scripts", "names.R"))
text_size <- 3
point_size <- 2

figure_1_dir <- file.path("data", "figure_1")
figure_2_dir <- file.path("data", "figure_2")
figure_2_exp_dir <- file.path("data", "figure_2_exp")
figure_mm10_dir <- file.path("data", "figure_mm10")
figure_aligners_dir <- file.path("data", "figure_aligners")
figure_aligners_2_dir <- file.path("data", "figure_aligners_2")
figure_parameters_dir <- file.path("data", "figure_parameters")
figure_merged_dir <- file.path("data", "figure_merged")

# load BAM and donor data
bam_data_path <- file.path("data", "bam_accession_table.csv")
donor_data_path <- file.path("data", "donor_accession_table.csv")

bam_data_df <- read.csv(bam_data_path, check.names = FALSE)
donor_data_df <- read.csv(donor_data_path, check.names = FALSE)

# Use this genome
genome <- "hg38"
#genome <- "mm10"

# Chromosomes to keep
include_chr <- c(
  "chr1", "chr2", "chr3", "chr4", "chr5",
  "chr6", "chr7", "chr8", "chr9", "chr10",
  "chr11", "chr12", "chr13", "chr14", "chr15",
  "chr16", "chr17", "chr18", "chr19", "chr20",
  "chr21", "chr22", "chrX", "chrY"
)

mm10_include_chr <- c(
  "chr1", "chr2", "chr3", "chr4", "chr5",
  "chr6", "chr7", "chr8", "chr9", "chr10",
  "chr11", "chr12", "chr13", "chr14", "chr15",
  "chr16", "chr17", "chr18", "chr19", "chrX", "chrY"
)

genome_size <- get_genome_size(
  genome_id = "hg38",
  only_autosomal = FALSE,
  include_chr = include_chr
)
mm10_genome_size <- get_genome_size(
  genome_id = "mm10",
  only_autosomal = FALSE,
  include_chr = include_chr
)

#gap_dir <- file.path("data", "hg38_gaps")
gap_dir <- file.path("data", paste0(genome, "_gaps"))
mm10_gap_dir <- file.path("data", paste0("mm10", "_gaps"))

exons_file <- file.path(".", "data", "exons", paste0(genome, "_exons.bed"))
genes_file <- file.path(".", "data", "genes", paste0(genome, "_genes.bed"))
cds_file <- file.path(".", "data", "cds", paste0(genome, "_cds.bed"))

mm10_exons_file <- file.path(".", "data", "exons", paste0("mm10", "_exons.bed"))
mm10_genes_file <- file.path(".", "data", "genes", paste0("mm10", "_genes.bed"))
mm10_cds_file <- file.path(".", "data", "cds", paste0("mm10", "_cds.bed"))

# Paths to centromeres, telomeres, and shortarms
path_to_centromeres <- file.path(
  gap_dir, paste0(genome, ".UCSC.centromere_merged.bed")
)
path_to_telomeres <- file.path(
  gap_dir, paste0(genome, ".UCSC.telomere.bed")
)
path_to_shortarms <- file.path(
  gap_dir, paste0(genome, ".UCSC.short_arms.bed")
)

mm10_path_to_centromeres <- file.path(
  mm10_gap_dir, paste0("mm10", ".UCSC.centromere_merged.bed")
)
mm10_path_to_telomeres <- file.path(
  mm10_gap_dir, paste0("mm10", ".UCSC.telomere.bed")
)
mm10_path_to_shortarms <- file.path(
  mm10_gap_dir, paste0("mm10", ".UCSC.short_arms.bed")
)

source(file.path("scripts", "main.R"))
try(source(file.path("scripts", "bio.R")), silent = TRUE)
#source(file.path("scripts", "embeddings_visuals.R"))
source(file.path("scripts", "names.R"))
#source(file.path("scripts", "ChIP_heatmaps.R"))
source(file.path("scripts", "GSEA.R"))
source(file.path("scripts", "BAMs.R"))
source(file.path("scripts", "signal.R"))
source(file.path("scripts", "summary_heatmaps.R"))
source(file.path("scripts", "RNA-seq_sponge.R"))
source(file.path("scripts", "WGS.R"))

## Figure 1
### A

exclusion_set_fig_1 <- names_dict[unique_names_fig_1]

# All
num_total_fig_1 <- get_num_regions(n_bams_list_fig_1$all)

# Gaps
num_overlap_gaps_fig_1 <- get_num_regions(n_bams_list_fig_1$overlap_gaps)

# All
t_w_total_fig_1 <- get_c_w(n_bams_list_fig_1$all)

# Gaps
t_w_overlap_gaps <- get_t_w(
  n_bams_list_fig_1$all,
  gaps_gr
)

gaps_counts_fig_1 <- data.frame(
  "e_s" = exclusion_set_fig_1,
  "total" = num_total_fig_1,
  "all_gaps" = num_overlap_gaps_fig_1,
  "not_gaps" = num_total_fig_1 - num_overlap_gaps_fig_1,

  # Keep column names unaltered
  check.names = FALSE
)
gaps_widths_fig_1 <- data.frame(
  "e_s" = exclusion_set_fig_1,
  "total" = t_w_total_fig_1,
  "all_gaps" = t_w_overlap_gaps,
  "not_gaps" = t_w_total_fig_1 - t_w_overlap_gaps,

  # Keep column names unaltered
  check.names = FALSE
)

plot_1a_1 <- get_plot_1a(gaps_counts_fig_1) +
  labs(
    y = "Number of Regions"
  )
plot_1a_2 <- get_plot_1a(gaps_widths_fig_1) +
  labs(
    y = "Coverage"
  )

bed_list_1b <- c("250.bed", "hg38-blacklistv2.bed", "KU.bed")

widths_1b <- sapply(bed_list_1b, function(l) {
  width(n_bams_list_fig_1$all[l][[1]])
})

df_1b <- data.frame(
  list = c(),
  widths = c()
)
for (i in seq_along(widths_1b)) {
  current_df <- data.frame(
    list = unname(
      names_dict[sub(paste0(".bed", "$"), "", names(widths_1b)[i])]
    ),
    widths = widths_1b[i]
  )
  names(current_df) <- c("list", "widths")
  df_1b <- rbind(df_1b, current_df)
}

df_1b$list <- factor(
  df_1b$list,
  levels = rev(unname(names_dict))
)

plot_1b <- ggplot(
  df_1b,
  aes(
    x = widths,
    y = list,
    fill = list
  )
) +
  ggridges::geom_density_ridges(
    color = NA
  ) +
  theme_minimal() +
  labs(
    title = NULL,
    y = NULL,
    x = expression(paste("Width (log", italic()[10], ")"))
  ) +
  scale_x_log10() +
  scale_fill_manual(
    name = NULL,
    values = list_colors_fig_1,
    breaks = NULL
  ) +
  theme(
    strip.text = element_blank(),
    strip.background = element_blank(),
    legend.position = "none"
  )

### C
bed_list_1c <- c("250.bed", "hg38-blacklistv2.bed", "KU.bed")
gr_list_1c <- n_bams_list_fig_1$all[bed_list_1c]

overlaps_1c <- lapply(seq_along(bed_list_1c), function(i) {
  gr_overlaps_apply(
    gr_list_1c[bed_list_1c[bed_list_1c != bed_list_1c[[i]]]],
    n_bams_list_fig_1$all[[bed_list_1c[[i]]]],
    io = FALSE
  )
})
names(overlaps_1c) <- bed_list_1c

counts_1c <- lapply(
  overlaps_1c, function(l_1) {
    lapply(l_1, function(l_2) {
      get_num_regions(l_2, unique_names_ = NULL)
    })
  }
)
widths_1c <- lapply(
  overlaps_1c, function(l_1) {
    lapply(l_1, function(l_2) {
      get_widths(l_2, unique_names_ = NULL)
    })
  }
)

plot_1c_1 <- make_bars_1c(
  counts_1c,
  mode = "counts"
) +
  labs(
    y = "Number of Regions"
  )
plot_1c_2 <- make_bars_1c(
  widths_1c,
  mode = "widths"
) +
  labs(
    y = "Width of Regions"
  )

cent_df <- num_overlap_cent_fig_1
telo_df <- num_overlap_telo_fig_1
s_arms_df <- num_overlap_s_arms_fig_1

n_bams_names_fig_1 <- unname(names_dict[
  sub(paste0(".bed", "$"), "", names(n_bams_all_fig_1))
])

plot_1d_1 <- make_bars_1d(
  num_overlap_cent_fig_1,
  num_overlap_telo_fig_1,
  num_overlap_s_arms_fig_1,
  n_bams_names_fig_1
) +
  labs(y = "Number of Regions")

plot_1d_2 <- make_bars_1d(
  t_w_overlap_cent_fig_1,
  t_w_overlap_telo_fig_1,
  t_w_overlap_s_arms_fig_1,
  n_bams_names_fig_1
) +
  labs(y = "Coverage")

fig_1a <- wrap_plots((
  plot_1a_1 +
    theme(legend.position = "bottom")
), (
  plot_1a_2 +
    theme(legend.position = "none")
), nrow = 2)
fig_1a[[1]] <- fig_1a[[1]] + labs(tag = "A")

fig_1b <- plot_1b + labs(tag = "D")

fig_1c <- wrap_plots((
  plot_1c_1 +
    theme(legend.position = "bottom")
), (
  plot_1c_2 +
    theme(legend.position = "none")
), nrow = 2)
fig_1c[[1]] <- fig_1c[[1]] + labs(tag = "C")

fig_1d <- wrap_plots((
  plot_1d_1 &
    theme(
      legend.position = "none",
      legend.direction = "vertical"
    )
), (
  plot_1d_2 &
    theme(legend.position = "none")
), nrow = 2)
fig_1d[[1]][[1]] <- fig_1d[[1]][[1]] + labs(tag = "B")

fig_1 <- (fig_1a | fig_1d | fig_1c) +
  plot_annotation(title = NULL)

export_fig(
  fig_1,
  "Figure_1"
)

## Figure 2
### A

exclusion_set_fig_2 <- names_dict[unique_names_fig_2]

# All
num_total_fig_2 <- get_num_regions(n_bams_list_fig_2$all)

# Gaps
num_overlap_gaps_fig_2 <- get_num_regions(n_bams_list_fig_2$overlap_gaps)

# All
t_w_total_fig_2 <- get_c_w(n_bams_list_fig_2$all)

# Gaps
t_w_overlap_gaps <- get_t_w(
  n_bams_list_fig_2$all,
  gaps_gr
)

gaps_counts_fig_2 <- data.frame(
  "e_s" = exclusion_set_fig_2,
  "total" = num_total_fig_2,
  "all_gaps" = num_overlap_gaps_fig_2,
  "not_gaps" = num_total_fig_2 - num_overlap_gaps_fig_2,

  # Keep column names unaltered
  check.names = FALSE
)
gaps_widths_fig_2 <- data.frame(
  "e_s" = exclusion_set_fig_2,
  "total" = t_w_total_fig_2,
  "all_gaps" = t_w_overlap_gaps,
  "not_gaps" = t_w_total_fig_2 - t_w_overlap_gaps,

  # Keep column names unaltered
  check.names = FALSE
)

gaps_df <- gaps_widths_fig_2
get_plot_2a <- function(gaps_df) {
  gaps_long <- gaps_df[
    c("e_s", "all_gaps", "not_gaps")
  ] %>%
    tidyr::pivot_longer(
      cols = c("all_gaps", "not_gaps"),
      names_to = "variable",
      values_to = "value"
    )
  gaps_long$e_s <- factor(gaps_long$e_s, levels = rev(unname(names_dict)))
  gaps_long$group <- ifelse(
    gaps_long$variable == "all_gaps",
    "Gaps",
    as.character(gaps_long$e_s)
  )
  gaps_long$group <- factor(
    gaps_long$group,
    levels = c("Gaps", setdiff(names(list_colors_fig_2), "Gaps"))
  )

  out_plot <- ggplot(
    gaps_long,
    aes(
      x = e_s,
      y = value,
      fill = group
    )
  ) +
    theme_minimal() +
    geom_bar(
      stat = "identity"
    ) +
    scale_fill_manual(
      name = NULL,
      values = list_colors_fig_2,
      breaks = c("Gaps")
    ) +
    theme(
      #axis.text.x = element_blank()
      strip.text = element_blank(),
      strip.background = element_blank(),
      legend.position = "bottom"
      #legend.position = "none"
    ) +
    labs(
      x = NULL,
      y = NULL
    ) +
    coord_flip()

  return(out_plot)
}

plot_2a_1 <- get_plot_2a(gaps_counts_fig_2) +
  labs(
    y = "Number of Regions"
  )
plot_2a_2 <- get_plot_2a(gaps_widths_fig_2) +
  labs(
    y = "Coverage"
  )

bed_list_2b <- names(n_bams_list_fig_2$all)

widths_2b <- sapply(bed_list_2b, function(l) {
  width(n_bams_list_fig_2$all[l][[1]])
})

df_2b <- data.frame(
  list = c(),
  widths = c()
)
for (i in seq_along(widths_2b)) {
  current_df <- data.frame(
    list = unname(
      names_dict[sub(paste0(".bed", "$"), "", names(widths_2b)[i])]
    ),
    widths = widths_2b[i]
  )
  names(current_df) <- c("list", "widths")
  df_2b <- rbind(df_2b, current_df)
}

df_2b$list <- factor(
  df_2b$list,
  levels = rev(unname(names_dict))
)

plot_2b <- ggplot(
  df_2b,
  aes(
    x = widths,
    y = list,
    fill = list
  )
) +
  ggridges::geom_density_ridges(
    color = NA
  ) +
  theme_minimal() +
  labs(
    title = NULL,
    y = NULL,
    x = expression(paste("Width (log", italic()[10], ")"))
  ) +
  scale_x_log10() +
  scale_fill_manual(
    name = NULL,
    values = list_colors_fig_2,
    breaks = NULL
  ) +
  theme(
    strip.text = element_blank(),
    strip.background = element_blank(),
    legend.position = "none"
  )

drop_lists <- c("250.bed")
bed_list_2c <- bed_list_fig_2[! bed_list_fig_2 %in% drop_lists]

overlaps_2c <- sapply(bed_list_2c, function(l) {
  gr_overlaps_apply_2(
    gr_list_fig_2[bed_list_fig_2[bed_list_fig_2 != l]],
    n_bams_list_fig_2$all[[l]],
    io = FALSE
  )
}, USE.NAMES = TRUE, simplify = FALSE)

counts_2c <- sapply(
  overlaps_2c,
  function(l_1) {
    sapply(l_1, function(l_2) {
      get_num_regions(l_2, unique_names_ = NULL)
    }, simplify = FALSE)
  },
  simplify = FALSE
)
widths_2c <- sapply(
  overlaps_2c,
  function(l_1) {
    sapply(l_1, function(l_2) {
      get_widths(l_2, unique_names_ = NULL)
    }, simplify = FALSE)
  },
  simplify = FALSE
)

in_list <- counts_2c
make_bars_2c <- function(
  in_list,
  mode = "counts"
) {

  in_combos <- combn(names(in_list), 2)
  combined_df <- data.frame(
    list = c(),
    Ex = c(),
    Shared = c(),
    pair = c()
  )

  # Loop through pairs
  i <- 1
  for (combo in split(in_combos, col(in_combos))) {
    bed_a_old_name <- combo[1]
    bed_a_new_name <- unname(
      names_dict[sub(paste0(".bed", "$"), "", bed_a_old_name)]
    )

    bed_b_old_name <- combo[2]
    bed_b_new_name <- unname(
      names_dict[sub(paste0(".bed", "$"), "", bed_b_old_name)]
    )

    # Extract values from in_list
    ex_a <- in_list[[bed_b_old_name]][["inside"]][[bed_a_old_name]]
    shared_a <- in_list[[bed_b_old_name]][["trans_a"]][[bed_a_old_name]]
    ex_b <- in_list[[bed_b_old_name]][["outside"]][[bed_a_old_name]]
    shared_b <- in_list[[bed_b_old_name]][["trans_b"]][[bed_a_old_name]]

    if (mode != "counts") {
      # True widths
      ex_a <- in_list[[bed_b_old_name]][["tw_a"]][[bed_a_old_name]]
      shared_a <- in_list[[bed_b_old_name]][["tw_c"]][[bed_a_old_name]]
      ex_b <- in_list[[bed_b_old_name]][["tw_b"]][[bed_a_old_name]]
      shared_b <- in_list[[bed_b_old_name]][["tw_c"]][[bed_a_old_name]]
    }

    current_df <- data.frame(
      list = c(bed_a_new_name, bed_b_new_name),
      Ex = c(ex_a, ex_b),
      Shared = c(shared_a, shared_b),
      pair = c(i, i)
    )
    combined_df <- rbind(combined_df, current_df)
    i <- i + 1
  }

  # Assign pair factor for grouping
  #df_a$pair <- factor(rownames(df_a), levels = rownames(df_a))

  combined_df_long <- combined_df %>%
    tidyr::pivot_longer(
      cols = c("Ex", "Shared"),
      names_to = "type",
      values_to = "value"
    )

  combined_df_long$list <- factor(
    combined_df_long$list,
    levels = rev(unname(names_dict))
  )

  combined_df_long$group <- ifelse(
    combined_df_long$type == "Shared",
    "Shared",
    as.character(combined_df_long$list)
  )
  ranking <- rank(
    as.integer(combined_df_long$list)[
      c(TRUE, FALSE, FALSE, FALSE)
    ] +
      as.integer(combined_df_long$list)[
        c(FALSE, FALSE, TRUE, FALSE)
      ],
    ties.method = "first"
  )
  combined_df_long$pair <- factor(
    combined_df_long$pair,
    levels = rev(order(ranking))
  )

  plot <- ggplot(
    combined_df_long,
    aes(
      x = list,
      y = value,
      fill = group
    )
  ) +
    geom_bar(
      stat = "identity",
      position = position_stack()
    ) +
    facet_grid(
      pair ~ .,
      scales = "free_y",
      space = "free_y"
    ) +
    scale_fill_manual(
      values = list_colors_fig_2,
      name = NULL,
      breaks = c("Shared")
    ) +
    labs(
      x = NULL,
      y = NULL
    ) +
    theme_minimal() +
    theme(
      #axis.text.x = element_blank(),
      strip.text = element_blank(),
      strip.background = element_blank(),
      legend.position = "bottom"
      #legend.position = "none"
    ) +
    coord_flip()

  return(plot)
}

plot_2c_1 <- make_bars_2c(
  counts_2c,
  mode = "counts"
) +
  labs(
    y = "Number of Regions"
  )
plot_2c_2 <- make_bars_2c(
  widths_2c,
  mode = "widths"
) +
  labs(
    y = "Width of Regions"
  )

### D
text_size <- 3
point_size <- 2

plot_2d_1 <- get_p_c_plots(
  j_c_distance_matrix_fig_2,
  "Jaccard Count Overlap",
  text_size,
  point_size,
  parameter_test = FALSE
)
plot_2d_2 <- get_p_c_plots(
  j_w_distance_matrix_fig_2,
  "Jaccard Width Overlap",
  text_size,
  point_size,
  parameter_test = FALSE
)
plot_2d_3 <- get_p_c_plots(
  f_w_distance_matrix_fig_2,
  "Forbes Width Overlap",
  text_size,
  point_size,
  parameter_test = FALSE
)
#plot_2d_embeddings <- embeddings_plot_list[["figure_S4D.csv"]] +
#  labs(
#    title = NULL,
#    subtitle = NULL
#  ) +
#  scale_color_manual(
#    values = list_colors_fig_2[
#      ggplot_build(
#        embeddings_plot_list[["figure_S4D.csv"]]
#      )$data[[2]]$label
#    ]
#  ) +
#  theme(
#    legend.position = "none"
#  ) + coord_fixed(ratio = 1)

fig_2a <- plot_2a_1 +
  theme(legend.position = "bottom") +
  labs(tag = "A")

fig_2b <- plot_2b + labs(tag = "C")

fig_2c <- plot_2c_1 +
  theme(legend.position = "bottom") +
  labs(tag = "B")

fig_2d <- plot_2d_1 +
  labs(tag = "D")

fig_2 <- (fig_2a / fig_2b | fig_2c / fig_2d) +
  plot_annotation(title = NULL)

export_fig(
  fig_2,
  "Figure_2"
)

## Figure 2 expanded
### A

exclusion_set_fig_2_exp <- names_dict[unique_names_fig_2_exp]

# All
num_total_fig_2_exp <- get_num_regions(n_bams_list_fig_2_exp$all)

# Gaps
num_overlap_gaps_fig_2_exp <- get_num_regions(n_bams_list_fig_2_exp$overlap_gaps)

# All
t_w_total_fig_2_exp <- get_c_w(n_bams_list_fig_2_exp$all)

# Gaps
t_w_overlap_gaps <- get_t_w(
  n_bams_list_fig_2_exp$all,
  gaps_gr
)

gaps_counts_fig_2_exp <- data.frame(
  "e_s" = exclusion_set_fig_2_exp,
  "total" = num_total_fig_2_exp,
  "all_gaps" = num_overlap_gaps_fig_2_exp,
  "not_gaps" = num_total_fig_2_exp - num_overlap_gaps_fig_2_exp,

  # Keep column names unaltered
  check.names = FALSE
)
gaps_widths_fig_2_exp <- data.frame(
  "e_s" = exclusion_set_fig_2_exp,
  "total" = t_w_total_fig_2_exp,
  "all_gaps" = t_w_overlap_gaps,
  "not_gaps" = t_w_total_fig_2_exp - t_w_overlap_gaps,

  # Keep column names unaltered
  check.names = FALSE
)

gaps_df <- gaps_widths_fig_2_exp
get_plot_2a <- function(gaps_df) {
  gaps_long <- gaps_df[
    c("e_s", "all_gaps", "not_gaps")
  ] %>%
    tidyr::pivot_longer(
      cols = c("all_gaps", "not_gaps"),
      names_to = "variable",
      values_to = "value"
    )
  gaps_long$e_s <- factor(gaps_long$e_s, levels = rev(unname(names_dict)))
  gaps_long$group <- ifelse(
    gaps_long$variable == "all_gaps",
    "Gaps",
    as.character(gaps_long$e_s)
  )
  gaps_long$group <- factor(
    gaps_long$group,
    levels = c("Gaps", setdiff(names(list_colors_fig_2_exp), "Gaps"))
  )

  out_plot <- ggplot(
    gaps_long,
    aes(
      x = e_s,
      y = value,
      fill = group
    )
  ) +
    theme_minimal() +
    geom_bar(
      stat = "identity"
    ) +
    scale_fill_manual(
      name = NULL,
      values = list_colors_fig_2_exp,
      breaks = c("Gaps")
    ) +
    theme(
      #axis.text.x = element_blank()
      strip.text = element_blank(),
      strip.background = element_blank(),
      legend.position = "bottom"
      #legend.position = "none"
    ) +
    labs(
      x = NULL,
      y = NULL
    ) +
    coord_flip()

  return(out_plot)
}

plot_2_exp_a_1 <- get_plot_2a(gaps_counts_fig_2_exp) +
  labs(
    y = "Number of Regions"
  )
plot_2_exp_a_2 <- get_plot_2a(gaps_widths_fig_2_exp) +
  labs(
    y = "Coverage"
  )

### B

bed_list_2b <- names(n_bams_list_fig_2_exp$all)

widths_2b <- sapply(bed_list_2b, function(l) {
  width(n_bams_list_fig_2_exp$all[l][[1]])
})

df_2b <- data.frame(
  list = c(),
  widths = c()
)
for (i in seq_along(widths_2b)) {
  current_df <- data.frame(
    list = unname(
      names_dict[sub(paste0(".bed", "$"), "", names(widths_2b)[i])]
    ),
    widths = widths_2b[i]
  )
  names(current_df) <- c("list", "widths")
  df_2b <- rbind(df_2b, current_df)
}

df_2b$list <- factor(
  df_2b$list,
  levels = rev(unname(names_dict))
)

plot_2_exp_b <- ggplot(
  df_2b,
  aes(
    x = widths,
    y = list,
    fill = list
  )
) +
  ggridges::geom_density_ridges(
    color = NA
  ) +
  theme_minimal() +
  labs(
    title = NULL,
    y = NULL,
    x = expression(paste("Width (log", italic()[10], ")"))
  ) +
  scale_x_log10() +
  scale_fill_manual(
    name = NULL,
    values = list_colors_fig_2_exp,
    breaks = NULL
  ) +
  theme(
    strip.text = element_blank(),
    strip.background = element_blank(),
    legend.position = "none"
  )

### C
drop_lists <- c("250.bed")
bed_list_2c <- bed_list_fig_2_exp[! bed_list_fig_2_exp %in% drop_lists]

overlaps_2c <- sapply(bed_list_2c, function(l) {
  gr_overlaps_apply_2(
    gr_list_fig_2_exp[bed_list_fig_2_exp[bed_list_fig_2_exp != l]],
    n_bams_list_fig_2_exp$all[[l]],
    io = FALSE
  )
}, USE.NAMES = TRUE, simplify = FALSE)

counts_2c <- sapply(
  overlaps_2c,
  function(l_1) {
    sapply(l_1, function(l_2) {
      get_num_regions(l_2, unique_names_ = NULL)
    }, simplify = FALSE)
  },
  simplify = FALSE
)
widths_2c <- sapply(
  overlaps_2c,
  function(l_1) {
    sapply(l_1, function(l_2) {
      get_widths(l_2, unique_names_ = NULL)
    }, simplify = FALSE)
  },
  simplify = FALSE
)

in_list <- counts_2c
make_bars_2c <- function(
  in_list,
  mode = "counts"
) {

  in_combos <- combn(names(in_list), 2)
  combined_df <- data.frame(
    list = c(),
    Ex = c(),
    Shared = c(),
    pair = c()
  )

  # Loop through pairs
  i <- 1
  for (combo in split(in_combos, col(in_combos))) {
    bed_a_old_name <- combo[1]
    bed_a_new_name <- unname(
      names_dict[sub(paste0(".bed", "$"), "", bed_a_old_name)]
    )

    bed_b_old_name <- combo[2]
    bed_b_new_name <- unname(
      names_dict[sub(paste0(".bed", "$"), "", bed_b_old_name)]
    )

    # Extract values from in_list
    ex_a <- in_list[[bed_b_old_name]][["inside"]][[bed_a_old_name]]
    shared_a <- in_list[[bed_b_old_name]][["trans_a"]][[bed_a_old_name]]
    ex_b <- in_list[[bed_b_old_name]][["outside"]][[bed_a_old_name]]
    shared_b <- in_list[[bed_b_old_name]][["trans_b"]][[bed_a_old_name]]

    if (mode != "counts") {
      # True widths
      ex_a <- in_list[[bed_b_old_name]][["tw_a"]][[bed_a_old_name]]
      shared_a <- in_list[[bed_b_old_name]][["tw_c"]][[bed_a_old_name]]
      ex_b <- in_list[[bed_b_old_name]][["tw_b"]][[bed_a_old_name]]
      shared_b <- in_list[[bed_b_old_name]][["tw_c"]][[bed_a_old_name]]
    }

    current_df <- data.frame(
      list = c(bed_a_new_name, bed_b_new_name),
      Ex = c(ex_a, ex_b),
      Shared = c(shared_a, shared_b),
      pair = c(i, i)
    )
    combined_df <- rbind(combined_df, current_df)
    i <- i + 1
  }

  # Assign pair factor for grouping
  #df_a$pair <- factor(rownames(df_a), levels = rownames(df_a))

  combined_df_long <- combined_df %>%
    tidyr::pivot_longer(
      cols = c("Ex", "Shared"),
      names_to = "type",
      values_to = "value"
    )

  combined_df_long$list <- factor(
    combined_df_long$list,
    levels = rev(unname(names_dict))
  )

  combined_df_long$group <- ifelse(
    combined_df_long$type == "Shared",
    "Shared",
    as.character(combined_df_long$list)
  )
  ranking <- rank(
    as.integer(combined_df_long$list)[
      c(TRUE, FALSE, FALSE, FALSE)
    ] +
      as.integer(combined_df_long$list)[
        c(FALSE, FALSE, TRUE, FALSE)
      ],
    ties.method = "first"
  )
  combined_df_long$pair <- factor(
    combined_df_long$pair,
    levels = rev(order(ranking))
  )

  plot <- ggplot(
    combined_df_long,
    aes(
      x = list,
      y = value,
      fill = group
    )
  ) +
    geom_bar(
      stat = "identity",
      position = position_stack()
    ) +
    facet_grid(
      pair ~ .,
      scales = "free_y",
      space = "free_y"
    ) +
    scale_fill_manual(
      values = list_colors_fig_2_exp,
      name = NULL,
      breaks = c("Shared")
    ) +
    labs(
      x = NULL,
      y = NULL
    ) +
    theme_minimal() +
    theme(
      #axis.text.x = element_blank(),
      strip.text = element_blank(),
      strip.background = element_blank(),
      legend.position = "bottom"
      #legend.position = "none"
    ) +
    coord_flip()

  return(plot)
}

plot_2_exp_c_1 <- make_bars_2c(
  counts_2c,
  mode = "counts"
) +
  labs(
    y = "Number of Regions"
  )
plot_2_exp_c_2 <- make_bars_2c(
  widths_2c,
  mode = "widths"
) +
  labs(
    y = "Width of Regions"
  )

### D
plot_2_exp_d_1 <- get_p_c_plots(
  j_c_distance_matrix_fig_2_exp,
  "Jaccard Count Overlap",
  text_size,
  point_size,
  parameter_test = FALSE
)
plot_2_exp_d_2 <- get_p_c_plots(
  j_w_distance_matrix_fig_2_exp,
  "Jaccard Width Overlap",
  text_size,
  point_size,
  parameter_test = FALSE
)
plot_2_exp_d_3 <- get_p_c_plots(
  f_w_distance_matrix_fig_2_exp,
  "Forbes Width Overlap",
  text_size,
  point_size,
  parameter_test = FALSE
)
#plot_2_exp_d_embeddings <- fig_2_expanded_plot +
#  labs(
#    title = NULL,
#    subtitle = NULL
#  ) +
#  scale_color_manual(
#    values = list_colors_fig_2_exp[
#      ggplot_build(
#        fig_2_expanded_plot
#      )$data[[2]]$label
#    ]
#  ) +
#  theme(
#    legend.position = "none"
#  )

fig_2_exp_a <- plot_2_exp_a_1 +
  theme(legend.position = "bottom") +
  labs(tag = "A")

fig_2_exp_b <- plot_2_exp_b + labs(tag = "C")

fig_2_exp_c <- plot_2_exp_c_1 +
  theme(legend.position = "bottom") +
  labs(tag = "B")

fig_2_exp_d <- plot_2_exp_d_1 +
  labs(tag = "D")

# Figure SX
fig_2_exp <- (fig_2_exp_a / fig_2_exp_b | fig_2_exp_c / fig_2_exp_d) +
  plot_annotation(title = NULL)

## Figure 2 - S2
n_bams_names_fig_2 <- unname(names_dict[
  sub(paste0(".bed", "$"), "", names(n_bams_all_fig_2))
])

plot_2_S2b <- make_bars_1d(
  t_w_overlap_cent_fig_2,
  t_w_overlap_telo_fig_2,
  t_w_overlap_s_arms_fig_2,
  n_bams_names_fig_2,
  colors_ = list_colors_fig_2
) +
  labs(y = "Coverage")

fig_2_S2a <- plot_2a_2 +
  theme(legend.position = "bottom") +
  labs(tag = "A")

fig_2_S2b <- wrap_plots(plot_2_S2b, ncol = 1) &
  theme(legend.position = "none")
fig_2_S2b[[1]][[1]] <- fig_2_S2b [[1]][[1]] +
  labs(tag = "C")

fig_2_S2c <- plot_2c_2 +
  theme(legend.position = "bottom") +
  labs(tag = "B")

fig_2_S2d <- plot_2d_3 +
  labs(tag = "D", title = "Forbes Width Overlap", subtitle = NULL)

# Supplementary Figure S4
fig_2_S4 <- (fig_2_S2a / fig_2_S2b | fig_2_S2c / fig_2_S2d) +
  plot_annotation(title = NULL)

export_fig(
  fig_2_S4,
  "Supplementary_Figure_S4"
)

## Figure 2 exp - S8
n_bams_names_fig_2_exp <- unname(names_dict[
  sub(paste0(".bed", "$"), "", names(n_bams_all_fig_2_exp))
])

plot_2_exp_S2b <- make_bars_1d(
  t_w_overlap_cent_fig_2_exp,
  t_w_overlap_telo_fig_2_exp,
  t_w_overlap_s_arms_fig_2_exp,
  n_bams_names_fig_2_exp,
  colors_ = list_colors_fig_2_exp
) +
  labs(y = "Coverage")

fig_2_exp_S2a <- plot_2_exp_a_2 +
  theme(legend.position = "bottom") +
  labs(tag = "A")

fig_2_exp_S2b <- wrap_plots(plot_2_exp_S2b, ncol = 1) +
  labs(tag = "C") &
  theme(legend.position = "none")

fig_2_exp_S2c <- plot_2_exp_c_2 +
  theme(legend.position = "bottom") +
  labs(tag = "B")

fig_2_exp_S2d <- plot_2_exp_d_3 +
  labs(tag = "D")

# "Supplementary Figure SX"
fig_2_exp_S4 <- (fig_2_exp_S2a / fig_2_exp_S2b | fig_2_exp_S2c / fig_2_exp_S2d) +
  plot_annotation(title = NULL)

## mm10 S1
### A

exclusion_set_fig_mm10 <- names_dict[unique_names_fig_mm10]

# All
num_total_fig_mm10 <- get_num_regions(n_bams_list_fig_mm10$all)

# Gaps
num_overlap_gaps_fig_mm10 <- get_num_regions(n_bams_list_fig_mm10$overlap_gaps)

# All
t_w_total_fig_mm10 <- get_c_w(n_bams_list_fig_mm10$all)

# Gaps
t_w_overlap_gaps_mm10 <- get_t_w(
  n_bams_list_fig_mm10$all,
  mm10_gaps_gr
)

gaps_counts_fig_mm10 <- data.frame(
  "e_s" = exclusion_set_fig_mm10,
  "total" = num_total_fig_mm10,
  "all_gaps" = num_overlap_gaps_fig_mm10,
  "not_gaps" = num_total_fig_mm10 - num_overlap_gaps_fig_mm10,

  # Keep column names unaltered
  check.names = FALSE
)
gaps_widths_fig_mm10 <- data.frame(
  "e_s" = exclusion_set_fig_mm10,
  "total" = t_w_total_fig_mm10,
  "all_gaps" = t_w_overlap_gaps_mm10,
  "not_gaps" = t_w_total_fig_mm10 - t_w_overlap_gaps_mm10,

  # Keep column names unaltered
  check.names = FALSE
)

plot_mm10_a_1 <- get_plot_1a(gaps_counts_fig_mm10, list_colors_fig_mm10) +
  labs(
    y = "Number of Regions"
  )
plot_mm10_a_2 <- get_plot_1a(gaps_widths_fig_mm10, list_colors_fig_mm10) +
  labs(
    y = "Coverage"
  )

### B

bed_list_mm10_b <- c("266_mm10og.bed", "mm10-blacklistv2.bed")

widths_mm10_b <- sapply(bed_list_mm10_b, function(l) {
  width(n_bams_list_fig_mm10$all[l][[1]])
})

df_mm10_b <- data.frame(
  list = c(),
  widths = c()
)
for (i in seq_along(widths_mm10_b)) {
  current_df <- data.frame(
    list = unname(
      names_dict[sub(paste0(".bed", "$"), "", names(widths_mm10_b)[i])]
    ),
    widths = widths_mm10_b[i]
  )
  names(current_df) <- c("list", "widths")
  df_mm10_b <- rbind(df_mm10_b, current_df)
}

df_mm10_b$list <- factor(
  df_mm10_b$list,
  levels = rev(unname(names_dict))
)

plot_mm10_b <- ggplot(
  df_mm10_b,
  aes(
    x = widths,
    y = list,
    fill = list
  )
) +
  ggridges::geom_density_ridges(
    color = NA
  ) +
  theme_minimal() +
  labs(
    title = NULL,
    y = NULL,
    x = expression(paste("Width (log", italic()[10], ")"))
  ) +
  scale_x_log10() +
  scale_fill_manual(
    name = NULL,
    values = list_colors_fig_mm10,
    breaks = NULL
  ) +
  theme(
    strip.text = element_blank(),
    strip.background = element_blank(),
    legend.position = "none"
  )

### C
bed_list_mm10_c <- c("266_mm10og.bed", "mm10-blacklistv2.bed")
gr_list_mm10_c <- n_bams_list_fig_mm10$all[bed_list_mm10_c]

overlaps_mm10_c <- lapply(seq_along(bed_list_mm10_c), function(i) {
  gr_overlaps_apply(
    gr_list_mm10_c[bed_list_mm10_c[bed_list_mm10_c != bed_list_mm10_c[[i]]]],
    n_bams_list_fig_mm10$all[[bed_list_mm10_c[[i]]]],
    io = FALSE
  )
})
names(overlaps_mm10_c) <- bed_list_mm10_c

counts_mm10_c <- lapply(
  overlaps_mm10_c, function(l_1) {
    lapply(l_1, function(l_2) {
      get_num_regions(l_2, unique_names_ = NULL)
    })
  }
)
widths_mm10_c <- lapply(
  overlaps_mm10_c, function(l_1) {
    lapply(l_1, function(l_2) {
      get_widths(l_2, unique_names_ = NULL)
    })
  }
)

plot_mm10_c_1 <- make_bars_1c(
  counts_mm10_c,
  mode = "counts",
  colors_ = list_colors_fig_mm10
) +
  labs(
    y = "Number of Regions"
  )
plot_mm10_c_2 <- make_bars_1c(
  widths_mm10_c,
  mode = "widths",
  colors_ = list_colors_fig_mm10
) +
  labs(
    y = "Width of Regions"
  )

### D
n_bams_names_fig_mm10 <- unname(names_dict[
  sub(paste0(".bed", "$"), "", names(n_bams_all_fig_mm10))
])

plot_mm10_d_1 <- make_bars_1d(
  num_overlap_cent_fig_mm10,
  num_overlap_telo_fig_mm10,
  num_overlap_s_arms_fig_mm10,
  n_bams_names_fig_mm10,
  list_colors_fig_mm10
) +
  labs(y = "Number of Regions")

plot_mm10_d_2 <- make_bars_1d(
  t_w_overlap_cent_fig_mm10,
  t_w_overlap_telo_fig_mm10,
  t_w_overlap_s_arms_fig_mm10,
  n_bams_names_fig_mm10,
  list_colors_fig_mm10
) +
  labs(y = "Coverage")

fig_mm10_a <- wrap_plots((
  plot_mm10_a_1 +
    theme(legend.position = "bottom")
), (
  plot_mm10_a_2 +
    theme(legend.position = "none")
), nrow = 2)
fig_mm10_a[[1]] <- fig_mm10_a[[1]] + labs(tag = "A")

fig_mm10_b <- plot_mm10_b + labs(tag = "D")

fig_mm10_c <- wrap_plots((
  plot_mm10_c_1 +
    theme(legend.position = "bottom")
), (
  plot_mm10_c_2 +
    theme(legend.position = "none")
), nrow = 2)
fig_mm10_c[[1]] <- fig_mm10_c[[1]] + labs(tag = "C")

fig_mm10_d <- wrap_plots((
  plot_mm10_d_1 &
    theme(
      legend.position = "none",
      legend.direction = "vertical"
    )
), (
  plot_mm10_d_2 &
    theme(legend.position = "none")
), nrow = 2)
fig_mm10_d[[1]][[1]] <- fig_mm10_d[[1]][[1]] + labs(tag = "B")

fig_mm10 <- (fig_mm10_a | fig_mm10_d | fig_mm10_c) +
  plot_annotation(title = NULL)

export_fig(
  fig_mm10,
  "Supplementary_Figure_S1"
)

## Aligners
### A

exclusion_set_fig_aligners <- names_dict[unique_names_fig_aligners]

# All
num_total_fig_aligners <- get_num_regions(n_bams_list_fig_aligners$all)

# Gaps
num_overlap_gaps_fig_aligners <- get_num_regions(n_bams_list_fig_aligners$overlap_gaps)

# All
t_w_total_fig_aligners <- get_c_w(n_bams_list_fig_aligners$all)

# Gaps
t_w_overlap_gaps <- get_t_w(
  n_bams_list_fig_aligners$all,
  gaps_gr
)

gaps_counts_fig_aligners <- data.frame(
  "e_s" = exclusion_set_fig_aligners,
  "total" = num_total_fig_aligners,
  "all_gaps" = num_overlap_gaps_fig_aligners,
  "not_gaps" = num_total_fig_aligners - num_overlap_gaps_fig_aligners,

  # Keep column names unaltered
  check.names = FALSE
)
gaps_widths_fig_aligners <- data.frame(
  "e_s" = exclusion_set_fig_aligners,
  "total" = t_w_total_fig_aligners,
  "all_gaps" = t_w_overlap_gaps,
  "not_gaps" = t_w_total_fig_aligners - t_w_overlap_gaps,

  # Keep column names unaltered
  check.names = FALSE
)

plot_aligners_a_1 <- get_plot_1a(
  gaps_counts_fig_aligners,
  list_colors_fig_aligners
) +
  labs(
    y = "Number of Regions"
  )
plot_aligners_a_2 <- get_plot_1a(
  gaps_widths_fig_aligners,
  list_colors_fig_aligners
) +
  labs(
    y = "Coverage"
  )

### B

bed_list_aligners_b <- names(n_bams_list_fig_aligners$all)

widths_aligners_b <- sapply(bed_list_aligners_b, function(l) {
  width(n_bams_list_fig_aligners$all[l][[1]])
})

df_aligners_b <- data.frame(
  list = c(),
  widths = c()
)
for (i in seq_along(widths_aligners_b)) {
  current_df <- data.frame(
    list = unname(
      names_dict[sub(paste0(".bed", "$"), "", names(widths_aligners_b)[i])]
    ),
    widths = widths_aligners_b[i]
  )
  names(current_df) <- c("list", "widths")
  df_aligners_b <- rbind(df_aligners_b, current_df)
}

df_aligners_b$list <- factor(
  df_aligners_b$list,
  levels = rev(unname(names_dict))
)

plot_aligners_b <- ggplot(
  df_aligners_b,
  aes(
    x = widths,
    y = list,
    fill = list
  )
) +
  ggridges::geom_density_ridges(
    color = NA
  ) +
  theme_minimal() +
  labs(
    title = NULL,
    y = NULL,
    x = expression(paste("Width (log", italic()[10], ")"))
  ) +
  scale_x_log10() +
  scale_fill_manual(
    name = NULL,
    values = list_colors_fig_aligners,
    breaks = NULL
  ) +
  theme(
    strip.text = element_blank(),
    strip.background = element_blank(),
    legend.position = "none"
  )

### C
bed_list_aligners_c <- names(n_bams_list_fig_aligners$all)
gr_list_aligners_c <- n_bams_list_fig_aligners$all[bed_list_aligners_c]
overlaps_aligners_c <- lapply(seq_along(bed_list_aligners_c), function(i) {
  gr_overlaps_apply(
    gr_list_aligners_c[bed_list_aligners_c[bed_list_aligners_c != bed_list_aligners_c[[i]]]],
    n_bams_list_fig_aligners$all[[bed_list_aligners_c[[i]]]],
    io = FALSE
  )
})
names(overlaps_aligners_c) <- bed_list_aligners_c

counts_aligners_c <- lapply(
  overlaps_aligners_c, function(l_1) {
    lapply(l_1, function(l_2) {
      get_num_regions(l_2, unique_names_ = NULL)
    })
  }
)
widths_aligners_c <- lapply(
  overlaps_aligners_c, function(l_1) {
    lapply(l_1, function(l_2) {
      get_widths(l_2, unique_names_ = NULL)
    })
  }
)

plot_aligners_c_1 <- make_bars_1c(
  counts_aligners_c,
  mode = "counts",
  list_colors_fig_aligners,
  focus_on = "hg38-blacklistv2.bed"
) +
  labs(
    y = "Number of Regions"
  )
plot_aligners_c_2 <- make_bars_1c(
  widths_aligners_c,
  mode = "widths",
  list_colors_fig_aligners,
  focus_on = "hg38-blacklistv2.bed"
) +
  labs(
    y = "Width of Regions"
  )

### D

cent_df <- num_overlap_cent_fig_aligners
telo_df <- num_overlap_telo_fig_aligners
s_arms_df <- num_overlap_s_arms_fig_aligners

n_bams_names_fig_aligners <- unname(names_dict[
  sub(paste0(".bed", "$"), "", names(n_bams_all_fig_aligners))
])

plot_aligners_d_1 <- make_bars_1d(
  num_overlap_cent_fig_aligners,
  num_overlap_telo_fig_aligners,
  num_overlap_s_arms_fig_aligners,
  n_bams_names_fig_aligners,
  list_colors_fig_aligners
) +
  labs(y = "Number of Regions")

plot_aligners_d_2 <- make_bars_1d(
  t_w_overlap_cent_fig_aligners,
  t_w_overlap_telo_fig_aligners,
  t_w_overlap_s_arms_fig_aligners,
  n_bams_names_fig_aligners,
  list_colors_fig_aligners
) +
  labs(y = "Coverage")

#plot_aligners_a_embeddings <- embeddings_plot_list[["figure_S5C.csv"]] +
#  labs(
#    title = "PCA of Embeddings"
#  ) +
#  scale_color_manual(
#    values = list_colors_fig_aligners[
#      ggplot_build(
#        embeddings_plot_list[["figure_S5C.csv"]]
#      )$data[[2]]$label
#    ]
#  ) +
#  theme(
#    legend.position = "none"
#  )

fig_aligners_a <- wrap_plots((
  plot_aligners_a_1 +
    theme(legend.position = "bottom")
), (
  plot_aligners_a_2 +
    theme(legend.position = "none")
), nrow = 2)
fig_aligners_a[[1]] <- fig_aligners_a[[1]] + labs(tag = "A")

fig_aligners_b <- plot_aligners_b + labs(tag = "D")

fig_aligners_c <- wrap_plots((
  plot_aligners_c_1 +
    theme(legend.position = "bottom")
), (
  plot_aligners_c_2 +
    theme(legend.position = "none")
), nrow = 2)
fig_aligners_c[[1]] <- fig_aligners_c[[1]] + labs(tag = "C")

fig_aligners_d <- wrap_plots((
  plot_aligners_d_1 &
    theme(
      legend.position = "none",
      legend.direction = "vertical"
    )
), (
  plot_aligners_d_2 &
    theme(legend.position = "none")
), nrow = 2)
fig_aligners_d[[1]][[1]] <- fig_aligners_d[[1]][[1]] + labs(tag = "B")

# Figure 3
fig_aligners <- (
  fig_aligners_a | fig_aligners_d | fig_aligners_c
) +
  plot_annotation(title = NULL)

export_fig(
  fig_aligners,
  "Figure_3"
)

plot_aln_mds_1 <- get_p_c_plots(
  j_c_distance_matrix_aln,
  "Jaccard Count Overlap",
  text_size,
  point_size,
  parameter_test = FALSE,
  colors_ = list_colors_fig_aligners
)
plot_aln_mds_2 <- get_p_c_plots(
  j_w_distance_matrix_aln,
  "Jaccard Width Overlap",
  text_size,
  point_size,
  parameter_test = FALSE,
  colors_ = list_colors_fig_aligners
)
plot_aln_mds_3 <- get_p_c_plots(
  f_w_distance_matrix_aln,
  "Forbes Width Overlap",
  text_size,
  point_size,
  parameter_test = FALSE,
  colors_ = list_colors_fig_aligners
)

plot_aln_den_1 <- get_dendro_plot(
  overlap_jaccard_count_hc_aln,
  "Hierarchical Clustering - Jaccard Count",
  "ward.D2"
)
plot_aln_den_2 <- get_dendro_plot(
  overlap_jaccard_width_hc_aln,
  "Hierarchical Clustering - Jaccard Width",
  "ward.D2"
)
plot_aln_den_3 <- get_dendro_plot(
  overlap_forbes_width_hc_aln,
  "Hierarchical Clustering - Forbes Width",
  "ward.D2"
)

fig_aligners_S2_a <- plot_aln_mds_1 +
  labs(
    title = "Jaccard Count Overlap",
    subtitle = NULL,
    tag = "A"
  )
fig_aligners_S2_b <- plot_aln_mds_3 +
  labs(
    title = "Forbes Width Overlap",
    subtitle = NULL,
    tag = "C"
  )
fig_aligners_S2_c <- plot_aln_den_1 +
  labs(
    title = NULL,
    subtitle = NULL,
    tag = "B"
  )
fig_aligners_S2_d <- plot_aln_den_3 +
  labs(
    title = NULL,
    subtitle = NULL,
    tag = "D"
  )

# Supplementary Figure S5
fig_aligners_S5 <- (
  fig_aligners_S2_a / fig_aligners_S2_b | fig_aligners_S2_c / fig_aligners_S2_d
) +
  plot_annotation(title = NULL)

export_fig(
  fig_aligners_S5,
  "Supplementary_Figure_S6"
)

#plot_fig_2_exp_a_embeddings <- fig_2_expanded_plot +
#  scale_color_manual(
#    values = list_colors_fig_2_exp[
#      ggplot_build(
#        fig_2_expanded_plot
#      )$data[[2]]$label
#    ]
#  ) +
#  theme(
#    legend.position = "none",
#  )


plot_fig_2_exp_mds_1 <- get_p_c_plots(
  j_c_distance_matrix_fig_2_exp,
  "Jaccard Count Overlap",
  text_size,
  point_size,
  parameter_test = FALSE,
  colors_ = list_colors_fig_2_exp
)

plot_fig_2_exp_mds_3 <- get_p_c_plots(
  f_w_distance_matrix_fig_2_exp,
  "Forbes Width Overlap",
  text_size,
  point_size,
  parameter_test = FALSE,
  colors_ = list_colors_fig_2_exp
)

plot_fig_2_exp_den_1 <- get_dendro_plot(
  overlap_jaccard_count_fig_2_exp_hc,
  "Hierarchical Clustering - Jaccard Count",
  "ward.D2"
)
plot_fig_2_exp_den_3 <- get_dendro_plot(
  overlap_forbes_width_fig_2_exp_hc,
  "Hierarchical Clustering - Forbes Width",
  "ward.D2"
)

fig_2_exp_S8_a <- plot_fig_2_exp_mds_1 +
  labs(
    subtitle = NULL,
    tag = "A"
  )

fig_2_exp_S8_b <- plot_fig_2_exp_mds_3 +
  labs(
    title = "Forbes Width Overlap",
    subtitle = NULL,
    tag = "C"
  )

fig_2_exp_S8_c <- plot_fig_2_exp_den_1 +
  labs(
    title = NULL,
    subtitle = NULL,
    tag = "B"
  )
fig_2_exp_S8_d <- plot_fig_2_exp_den_3 +
  labs(
    title = NULL,
    subtitle = NULL,
    tag = "D"
  )

# Supplementary Figure S8
fig_2_exp_S8 <- (
  fig_2_exp_S8_a / fig_2_exp_S8_b | fig_2_exp_S8_c / fig_2_exp_S8_d
) +
  plot_annotation(title = NULL)

export_fig(
  fig_2_exp_S8,
  "Supplementary_Figure_S11"
)

### S4

exclusion_set_fig_aligners_2 <- names_dict[unique_names_fig_aligners_2]

# All
num_total_fig_aligners_2 <- get_num_regions(n_bams_list_fig_aligners_2$all)

# Gaps
num_overlap_gaps_fig_aligners_2 <- get_num_regions(n_bams_list_fig_aligners_2$overlap_gaps)

# All
t_w_total_fig_aligners_2 <- get_c_w(n_bams_list_fig_aligners_2$all)

# Gaps
t_w_overlap_gaps <- get_t_w(
  n_bams_list_fig_aligners_2$all,
  gaps_gr
)

gaps_counts_fig_aligners_2 <- data.frame(
  "e_s" = exclusion_set_fig_aligners_2,
  "total" = num_total_fig_aligners_2,
  "all_gaps" = num_overlap_gaps_fig_aligners_2,
  "not_gaps" = num_total_fig_aligners_2 - num_overlap_gaps_fig_aligners_2,

  # Keep column names unaltered
  check.names = FALSE
)
gaps_widths_fig_aligners_2 <- data.frame(
  "e_s" = exclusion_set_fig_aligners_2,
  "total" = t_w_total_fig_aligners_2,
  "all_gaps" = t_w_overlap_gaps,
  "not_gaps" = t_w_total_fig_aligners_2 - t_w_overlap_gaps,

  # Keep column names unaltered
  check.names = FALSE
)

plot_aligners_2_a_1 <- get_plot_1a(
  gaps_counts_fig_aligners_2,
  list_colors_fig_aligners_2
) +
  labs(
    y = "Number of Regions"
  )
plot_aligners_2_a_2 <- get_plot_1a(
  gaps_widths_fig_aligners_2,
  list_colors_fig_aligners_2
) +
  labs(
    y = "Coverage"
  )

plot_aln_2_mds_1 <- get_p_c_plots(
  j_c_distance_matrix_aln_2,
  "Jaccard Count Overlap",
  text_size,
  point_size,
  parameter_test = FALSE,
  colors_ = list_colors_fig_aligners_2
)
plot_aln_2_mds_2 <- get_p_c_plots(
  j_w_distance_matrix_aln_2,
  "Jaccard Width Overlap",
  text_size,
  point_size,
  parameter_test = FALSE,
  colors_ = list_colors_fig_aligners_2
)
plot_aln_2_mds_3 <- get_p_c_plots(
  f_w_distance_matrix_aln_2,
  "Forbes Width Overlap",
  text_size,
  point_size,
  parameter_test = FALSE,
  colors_ = list_colors_fig_aligners_2
)
#plot_S6_b_embeddings <- embeddings_plot_list[["figure_S6B.csv"]] +
#  labs(
#    title = NULL,
#    subtitle = NULL
#  ) +
#  scale_color_manual(
#    values = list_colors_fig_aligners_2[
#      ggplot_build(
#        embeddings_plot_list[["figure_S6B.csv"]]
#      )$data[[2]]$label
#    ]
#  ) +
#  theme(
#    legend.position = "none"
#  )

fig_aln_S4_a <- (plot_aligners_2_a_1 + labs(tag = "A")) /
  plot_aligners_2_a_2 +
  theme(legend.position = "none")
fig_aln_S4_b <- (
  plot_aln_2_mds_1 +
    labs(tag = "B")
) /
  plot_aln_2_mds_3 +
    labs(title = "Forbes Width Overlap", subtitle = NULL)

# Supplementary Figure S6
fig_aln_S6 <- (fig_aln_S4_a | fig_aln_S4_b) +
  plot_annotation(title = NULL)

export_fig(
  fig_aln_S6,
  "Supplementary_Figure_S7"
)

## Merged
### A

exclusion_set_fig_merged <- names_dict[unique_names_fig_merged]

# All
num_total_fig_merged <- get_num_regions(n_bams_list_fig_merged$all)

# Gaps
num_overlap_gaps_fig_merged <- get_num_regions(n_bams_list_fig_merged$overlap_gaps)

# All
t_w_total_fig_merged <- get_c_w(n_bams_list_fig_merged$all)

# Gaps
t_w_overlap_gaps_merged <- get_t_w(
  n_bams_list_fig_merged$all,
  gaps_gr
)

gaps_counts_fig_merged <- data.frame(
  "e_s" = exclusion_set_fig_merged,
  "total" = num_total_fig_merged,
  "all_gaps" = num_overlap_gaps_fig_merged,
  "not_gaps" = num_total_fig_merged - num_overlap_gaps_fig_merged,

  # Keep column names unaltered
  check.names = FALSE
)
gaps_widths_fig_merged <- data.frame(
  "e_s" = exclusion_set_fig_merged,
  "total" = t_w_total_fig_merged,
  "all_gaps" = t_w_overlap_gaps_merged,
  "not_gaps" = t_w_total_fig_merged - t_w_overlap_gaps_merged,

  # Keep column names unaltered
  check.names = FALSE
)

plot_merged_a_1 <- get_plot_1a(gaps_counts_fig_merged, list_colors_fig_merged) +
  labs(
    y = "Number of Regions"
  )
plot_merged_a_2 <- get_plot_1a(gaps_widths_fig_merged, list_colors_fig_merged) +
  labs(
    y = "Coverage"
  )

### B

bed_list_merged_b <- names(n_bams_list_fig_merged$all)

widths_merged_b <- sapply(bed_list_merged_b, function(l) {
  width(n_bams_list_fig_merged$all[l][[1]])
})

df_merged_b <- data.frame(
  list = c(),
  widths = c()
)
for (i in seq_along(widths_merged_b)) {
  current_df <- data.frame(
    list = unname(
      names_dict[sub(paste0(".bed", "$"), "", names(widths_merged_b)[i])]
    ),
    widths = widths_merged_b[i]
  )
  names(current_df) <- c("list", "widths")
  df_merged_b <- rbind(df_merged_b, current_df)
}

df_merged_b$list <- factor(
  df_merged_b$list,
  levels = rev(unname(names_dict))
)

plot_merged_b <- ggplot(
  df_merged_b,
  aes(
    x = widths,
    y = list,
    fill = list
  )
) +
  ggridges::geom_density_ridges(
    color = NA
  ) +
  theme_minimal() +
  labs(
    title = NULL,
    y = NULL,
    x = expression(paste("Width (log", italic()[10], ")"))
  ) +
  scale_x_log10() +
  scale_fill_manual(
    name = NULL,
    values = list_colors_fig_merged,
    breaks = NULL
  ) +
  theme(
    strip.text = element_blank(),
    strip.background = element_blank(),
    legend.position = "none"
  )

### C
bed_list_merged_c <- names(n_bams_list_fig_merged$all)
gr_list_merged_c <- n_bams_list_fig_merged$all[bed_list_merged_c]

overlaps_merged_c <- lapply(seq_along(bed_list_merged_c), function(i) {
  gr_overlaps_apply(
    gr_list_merged_c[bed_list_merged_c[bed_list_merged_c != bed_list_merged_c[[i]]]],
    n_bams_list_fig_merged$all[[bed_list_merged_c[[i]]]],
    io = FALSE
  )
})
names(overlaps_merged_c) <- bed_list_merged_c

counts_merged_c <- lapply(
  overlaps_merged_c, function(l_1) {
    lapply(l_1, function(l_2) {
      get_num_regions(l_2, unique_names_ = NULL)
    })
  }
)
widths_merged_c <- lapply(
  overlaps_merged_c, function(l_1) {
    lapply(l_1, function(l_2) {
      get_widths(l_2, unique_names_ = NULL)
    })
  }
)

plot_merged_c_1 <- make_bars_1c(
  counts_merged_c,
  mode = "counts",
  colors_ = list_colors_fig_merged
) +
  labs(
    y = "Number of Regions"
  )
plot_merged_c_2 <- make_bars_1c(
  widths_merged_c,
  mode = "widths",
  colors_ = list_colors_fig_merged
) +
  labs(
    y = "Width of Regions"
  )

### D

n_bams_names_fig_merged <- unname(names_dict[
  sub(paste0(".bed", "$"), "", names(n_bams_all_fig_merged))
])

plot_merged_d_1 <- make_bars_1d(
  num_overlap_cent_fig_merged,
  num_overlap_telo_fig_merged,
  num_overlap_s_arms_fig_merged,
  n_bams_names_fig_merged,
  list_colors_fig_merged
) +
  labs(y = "Number of Regions")

plot_merged_d_2 <- make_bars_1d(
  t_w_overlap_cent_fig_merged,
  t_w_overlap_telo_fig_merged,
  t_w_overlap_s_arms_fig_merged,
  n_bams_names_fig_merged,
  list_colors_fig_merged
) +
  labs(y = "Coverage")

fig_merged_a <- wrap_plots((
  plot_merged_a_1 +
    theme(legend.position = "bottom")
), (
  plot_merged_a_2 +
    theme(legend.position = "none")
), nrow = 2)
fig_merged_a[[1]] <- fig_merged_a[[1]] + labs(tag = "A")

#fig_merged_b <- plot_merged_b + labs(tag = "D")

fig_merged_c <- wrap_plots((
  plot_merged_c_1 +
    theme(legend.position = "bottom")
), (
  plot_merged_c_2 +
    theme(legend.position = "none")
), nrow = 2)
fig_merged_c[[1]] <- fig_merged_c[[1]] + labs(tag = "C")

fig_merged_d <- wrap_plots((
  plot_merged_d_1 &
    theme(
      legend.position = "none",
      legend.direction = "vertical"
    )
), (
  plot_merged_d_2 &
    theme(legend.position = "none")
), nrow = 2)
fig_merged_d[[1]][[1]] <- fig_merged_d[[1]][[1]] + labs(tag = "B")

# Supplementary Figure S3
fig_merged <- (fig_merged_a | fig_merged_d | fig_merged_c) +
  plot_annotation(title = NULL)

export_fig(
  fig_merged,
  "Supplementary_Figure_S3"
)

## Parameters
### Heatmaps
ps_counts <- get_num_regions(n_bams_list_fig_parameters$all)
ps_sums <- get_c_x(n_bams_list_fig_parameters$all, sum)
ps_means <- get_c_x(n_bams_list_fig_parameters$all, mean)
ps_medians <- get_c_x(n_bams_list_fig_parameters$all, median)

ps_df <- data.frame(
  list = names_dict[sub(paste0(".bed", "$"), "", names(ps_counts))],
  "Number of Regions" = ps_counts,
  "Sum of Region Widths" = ps_sums,
  "Mean of Region Widths" = ps_means,
  "Median of Region Widths" = ps_medians,
  check.names = FALSE
)

labels <- gsub("STAR 36-BAMs ", "", ps_df$list)
labels <- gsub("20000", "20k", labels)
labels <- gsub("10000", "10k", labels)
labels <- gsub("1000", "1k", labels)
labels <- gsub("b", "", labels)
labels <- gsub(" k", " ", labels)
labels <- gsub("n", "", labels)

labels_split <- strsplit(labels, " ")
labels_split_matrix <- do.call(rbind, labels_split)
labels_df <- as.data.frame(labels_split_matrix, stringsAsFactors = FALSE)
rownames(labels_df) <- ps_df$list
colnames(labels_df) <- c("Bridge", "K-mer", "BAM Files")
labels_df$Bridge <- factor(
  labels_df$Bridge,
  levels = c("1k", "10k", "20k")
)
labels_df$`K-mer` <- factor(
  labels_df$`K-mer`,
  levels = c("36", "50", "100")
)
labels_df$`BAM Files` <- factor(
  labels_df$`BAM Files`,
  levels = c("10", "50", "100", "200", "300")
)

ps_df <- cbind(ps_df, labels_df)
ps_df <- na.exclude(ps_df)

text_size <- 1

ratio_1 <- 0.4
ratio_2 <- 0.5

ps_plots <- sapply(
  c(
    "Number of Regions",
    "Sum of Region Widths",
    "Mean of Region Widths",
    "Median of Region Widths"
  ),
  function(z) {
    # Number of files vs. bridge, separate for k-mer 36, 50, 100
    p1 <- get_ps_plots(
      ps_df,
      "K-mer",
      "BAM Files",
      "Bridge",
      z,
      text_size,
      ratio_1
    )

    # Number of files vs. k-mer, separate for bridge 1K, 10K, 20K
    p2 <- get_ps_plots(
      ps_df,
      "Bridge",
      "BAM Files",
      "K-mer",
      z,
      text_size,
      ratio_1
    )

    # Bridge vs. k-mer, separate for the different number of files
    p3 <- get_ps_plots(
      ps_df,
      "BAM Files",
      "Bridge",
      "K-mer",
      z,
      text_size,
      ratio_2
    )

    p1 <- clean_row_labels(p1, "K-mer: ")
    p1 <- wrap_plots(p1, nrow = 1)

    p2 <- clean_row_labels(p2, "Bridge: ")
    p2 <- wrap_plots(p2, nrow = 1)

    p3 <- clean_row_labels(p3, "BAM Files: ")
    p3 <- wrap_plots(p3, nrow = 1)

    out_plot <- wrap_plots(p1, p2, p3, ncol = 1) +
      plot_annotation(title = z)
    return(wrap_elements(panel = out_plot))
  },
  simplify = FALSE
)

# Supplementary Figure S7
fig_S7 <- wrap_plots(ps_plots) +
  plot_annotation(title = NULL)

export_fig(
  fig_S7,
  "Supplementary_Figure_S8"
)

### MDS
text_size <- 3
point_size <- 1

j_c_bam_plot <- get_para_mds_plots(
  in_matrix = j_c_distance_matrix_para,
  title_ = "Jaccard Count MDS",
  text_size = text_size,
  point_size = point_size,
  color_by = "BAM Files"
)
j_c_bridge_plot <- get_para_mds_plots(
  in_matrix = j_c_distance_matrix_para,
  title_ = "Jaccard Count MDS",
  text_size = text_size,
  point_size = point_size,
  color_by = "Bridge"
)
j_c_kmer_plot <- get_para_mds_plots(
  in_matrix = j_c_distance_matrix_para,
  title_ = "Jaccard Count MDS",
  text_size = text_size,
  point_size = point_size,
  color_by = "K-mer"
)

f_w_bam_plot <- get_para_mds_plots(
  in_matrix = f_w_distance_matrix_para,
  title_ = "Forbes Width MDS",
  text_size = text_size,
  point_size = point_size,
  color_by = "BAM Files"
)
f_w_bridge_plot <- get_para_mds_plots(
  in_matrix = f_w_distance_matrix_para,
  title_ = "Forbes Width MDS",
  text_size = text_size,
  point_size = point_size,
  color_by = "Bridge"
)
f_w_kmer_plot <- get_para_mds_plots(
  in_matrix = f_w_distance_matrix_para,
  title_ = "Forbes Width MDS",
  text_size = text_size,
  point_size = point_size,
  color_by = "K-mer"
)

#df_4 <- ggplot_build(
#  embeddings_plot_list[["figure_4H.csv"]]
#)$data[[1]][c("x", "y", "label")]

#e_bam_plot <- get_para_mds_plots_e(
#  in_df = df_4,
#  title_ = "Embeddings",
#  text_size = text_size,
#  point_size = point_size,
#  color_by = "BAM Files",
#  prefix_ = "STAR 36bp "
#)
#e_bridge_plot <- get_para_mds_plots_e(
#  in_df = df_4,
#  title_ = "Embeddings",
#  text_size = text_size,
#  point_size = point_size,
#  color_by = "Bridge",
#  prefix_ = "STAR 36bp "
#)
#e_kmer_plot <- get_para_mds_plots_e(
#  in_df = df_4,
#  title_ = "Embeddings",
#  text_size = text_size,
#  point_size = point_size,
#  color_by = "K-mer",
#  prefix_ = "STAR 36bp "
#)

f_w_rownames <- rownames(cmdscale(f_w_distance_matrix_para_gs))
f_w_labels <- ifelse(
  f_w_rownames %in% names(list_colors_fig_2),
  f_w_rownames,
  ""
)

f_w_labels_plot <-  get_para_mds_plots(
  in_matrix = f_w_distance_matrix_para_gs,
  title_ = "Forbes Width MDS",
  text_size = text_size,
  point_size = point_size,
  color_by = "label"
) +
  theme(
    legend.position = "none"
  ) +
  scale_color_manual(
    values = list_colors_fig_2
  ) +
  geom_text_repel(
    label = f_w_labels,
    color = "black",
    size = 3,
    max.overlaps = getOption(
      "ggrepel.max.overlaps",
      default = 100
    ),
    label.padding = 0.5
  )

#e_labels_plot <- embeddings_plot_before +
#  scale_color_manual(
#    values = list_colors_fig_2
#  ) +
#  theme(
#    legend.position = "none"
#  ) +
#  labs(
#    title = "Embeddings",
#    subtitle = "Principal Coordinate Analysis"
#  ) +
#  coord_cartesian(clip = "off")

j_c_rownames <- rownames(cmdscale(j_c_distance_matrix_para_gs))
j_c_labels <- ifelse(
  j_c_rownames %in% names(list_colors_fig_2),
  j_c_rownames,
  ""
)

j_c_labels_plot <-  get_para_mds_plots(
  in_matrix = j_c_distance_matrix_para_gs,
  title_ = "Jaccard Count MDS",
  text_size = text_size,
  point_size = point_size,
  color_by = "label"
) +
  theme(
    legend.position = "none"
  ) +
  scale_color_manual(
    values = list_colors_fig_2
  ) +
  geom_text_repel(
    label = j_c_labels,
    color = "black",
    size = 3,
    max.overlaps = getOption(
      "ggrepel.max.overlaps",
      default = 100
    ),
    label.padding = 0.5
  )

fig_para_mds <- wrap_plots(
  j_c_bam_plot +
    labs(tag = "A"),
  j_c_bridge_plot +
    labs(tag = "B"),
  j_c_kmer_plot +
    labs(tag = "C"),
  j_c_labels_plot +
    labs(tag = "D"),
  f_w_bam_plot +
    labs(tag = "E"),
  f_w_bridge_plot +
    labs(tag = "F"),
  f_w_kmer_plot +
    labs(tag = "G"),
  f_w_labels_plot +
    labs(tag = "H"),
  nrow = 2
) +
  plot_annotation(title = NULL) &
  labs(title = NULL, subtitle = NULL)

export_fig(
  fig_para_mds,
  "Supplementary_Figure_S9",
  x = 7,
  y = 3
)

## Biological Characterization

png_a <- file.path("images", "igv_snapshot.png")
png_b <- file.path("images", "igv_snapshot_RNA.png")

get_ratio <- function(image_path) {
  img <- image_read(image_path)
  image_info <- image_info(img)
  aspect_ratio <- image_info$width / image_info$height
  return(aspect_ratio)
}

ratio_a <- get_ratio(png_a)
ratio_b <- get_ratio(png_b)

igv_plot_a <- ggdraw(xlim = c(0, ratio_a), ylim = c(0, 1)) +
  draw_image(png_a, x = 0, y = 0, width = ratio_a, height = 1)
igv_plot_b <- ggdraw(xlim = c(0, ratio_b), ylim = c(0, 1)) +
  draw_image(png_b, x = 0, y = 0, width = ratio_b, height = 1)

fig_bio_a <- fig_bio_main &
  scale_fill_manual(
    values = list_colors_fig_2_exp
  )
fig_bio_a[[1]] <- fig_bio_a[[1]] + labs(tag = "A")

# Adjust `fig_bio_b_a` and `fig_bio_b_b` by removing extra space
fig_bio_b_a <- igv_plot_a +
  labs(tag = "B") +
  theme(
    plot.tag = element_text(face = "plain")
  )

fig_bio_b_b <- igv_plot_b +
  theme(
    plot.tag = element_text(face = "plain")
  )
fig_bio_b <- (fig_bio_b_a / fig_bio_b_b)

fig_bio_c <- dot_plot + labs(tag = "C")

# Assuming fig_bio_a, fig_bio_b, and fig_bio_c are your individual plots
fig_bio <- (
  free(fig_bio_a, side = "b") | fig_bio_c
) / free(fig_bio_b) +
  plot_layout(heights = c(0.5, 1))

export_fig(
  fig_bio,
  "Figure_4",
  x = 6,
  y = 5
)

## TF Heatmaps (summary_heatmaps.R)
export_fig(
  chip_final_plot,
  "Figure_5",
  x = 6,
  y = 4
)

# TF Boxplots (S10) (summary_heatmaps.R)
export_fig(
  tf_boxplots_all,
  "Supplementary_Figure_S10",
  scale = 3
)

## BAMs (S2)

png_bams <- file.path("images", "Read Length.png")

ratio_bams <- get_ratio(png_bams)

pie_plot <- ggdraw(xlim = c(0, ratio_bams), ylim = c(0, 1)) +
  draw_image(png_bams, x = 0, y = 0, width = ratio_bams, height = 1)

fig_bams_a <- pie_plot +
  theme(
    plot.tag = element_text(face = "plain")
  ) +
  labs(tag = "A")
fig_bams_b <- plot_mapped_reads + labs(tag = "B")
fig_bams_c <- plot_npd_counts + labs(tag = "C")
fig_bams_d <- plot_npd_reads + labs(tag = "D")

fig_bams <- (free(fig_bams_a) + fig_bams_b) /
  (fig_bams_c + fig_bams_d)

export_fig(
  fig_bams,
  "Supplementary_Figure_S2"
)

# S11 (RNA-seq_sponge.R and signal.R)

signal_plot[[1]][[1]] <- signal_plot[[1]][[1]] +
  labs(tag = "A")

fig_S11_a <- signal_plot

fig_S11_b <- (
  (rna_count_plot + labs(tag = "B")) / filtered_plot +
    theme(legend.box.margin = margin(0, 0, 0, -160))
)

fig_S11 <- wrap_plots(
  list(fig_S11_a, fig_S11_b),
  nrow = 1,
  widths = c(3, 1)
)

export_fig(
  fig_S11,
  "Supplementary_Figure_S12"
)
