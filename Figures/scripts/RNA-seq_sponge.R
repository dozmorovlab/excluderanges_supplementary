# Set up the environment
library(pander)
panderOptions('table.split.table', Inf)
set.seed(1)

library(tidyverse)
library(data.table)
library(ggplot2)

setwd(file.path(here::here(), "Figures"))

path_full <- file.path("data", "bulkRNASeq", "OUT_full", "star_rsem")
path_full_all <- file.path("data", "bulkRNASeq", "OUT_full_all", "star_rsem")

# File types to process
file_types <- c(
  "rsem.merged.gene_counts.tsv",
  "rsem.merged.gene_tpm.tsv",
  "rsem.merged.transcript_counts.tsv",
  "rsem.merged.transcript_tpm.tsv"
)

# Total count effect

# Function to read files and calculate column sums
process_files <- function(path, file_pattern) {
  files <- list.files(path, pattern = file_pattern, full.names = TRUE)

  # Read the file
  data <- fread(files)

  # Calculate column sums for all sample columns
  # (excluding gene/transcript ID columns)
  sample_cols <- names(data)[!grepl("_id", names(data))]
  col_sums <- colSums(data[, ..sample_cols])

  return(col_sums)
}

# Process gene counts
gene_counts_full <- process_files(
  path_full,
  "rsem.merged.gene_counts.tsv"
)
gene_counts_full_all <- process_files(
  path_full_all,
  "rsem.merged.gene_counts.tsv"
)

# Process gene TPM
gene_tpm_full <- process_files(path_full, "rsem.merged.gene_tpm.tsv")
gene_tpm_full_all <- process_files(path_full_all, "rsem.merged.gene_tpm.tsv")

# Create dataframes for statistical tests and plotting
gene_counts_df <- data.frame(
  full = gene_counts_full,
  full_all = gene_counts_full_all,
  group = c(rep("Cont", length(gene_counts_full)))
)

gene_tpm_df <- data.frame(
  full = gene_tpm_full,
  full_all = gene_tpm_full_all,
  group = c(rep("Cont", length(gene_counts_full)))
)

# Statistical tests
t_test_counts <- t.test(gene_counts_df$full, gene_counts_df$full_all)
wilcox_test_counts <- wilcox.test(gene_counts_df$full, gene_counts_df$full_all)
ks_test_counts <- ks.test(gene_counts_df$full, gene_counts_df$full_all)

t_test_tpm <- t.test(gene_tpm_df$full, gene_tpm_df$full_all)
wilcox_test_tpm <- wilcox.test(gene_tpm_df$full, gene_tpm_df$full_all)
ks_test_tpm <- ks.test(gene_tpm_df$full, gene_tpm_df$full_all)

# Visualization
# Prepare data for boxplot
boxplot_data_count <- data.frame(
  type = factor(rep(
    c("hg38", "hg38 & Sponge"), each = length(gene_counts_full)
  )),
  value = c(
    gene_counts_full,    # full Cont
    gene_counts_full_all# full_all Cont
  )
)
# Create boxplot
rna_count_plot <- ggplot(
  boxplot_data_count,
  aes(x = type, y = value, fill = type)
) +
  geom_boxplot() +
  labs(
    title = NULL,
    subtitle = "Count",
    x = NULL,
    y = "Total Sum"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 1),
    legend.position = "none"
  )

boxplot_data_tpm <- data.frame(
  type = factor(rep(
    c("hg38", "hg38 & Sponge"), each = length(gene_counts_full)
  )),
  value = c(
    gene_tpm_full,       # full TPM
    gene_tpm_full_all    # full_all TPM
  )
)

# Create boxplot
rna_tpm_plot <- ggplot(
  boxplot_data_tpm,
  aes(x = type, y = value, fill = type)
) +
  geom_boxplot() +
  labs(
    title = NULL,
    subtitle = "TPM",
    x = NULL,
    y = "Total Sum"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 1),
    legend.position = "none"
  )
