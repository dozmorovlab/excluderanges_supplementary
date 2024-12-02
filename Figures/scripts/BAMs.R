library(ggplot2)
library(dplyr)
library(tidyr)

bam_data_df_mr <- bam_data_df %>%
  mutate(Read_Type = ifelse(`Paired-End` == "True", "Paired End", "Single End"))

plot_mapped_reads <- ggplot(
  bam_data_df_mr,
  aes(x = `Mapped Reads`, y = Read_Type)
) +
  geom_violin(
    scale = "width",
    fill = "black"
  ) +
  scale_x_log10() +
  labs(
    title = NULL,
    x = "Number of Mapped Reads (log₁₀)",
    y = NULL
  ) +
  theme_minimal() +
  theme(legend.position = "none")

donor_data_df_npd <- donor_data_df %>%
  mutate(
    `Total BAM count` = `Single-End BAM count` + `Paired-End BAM count`
  ) %>%
  pivot_longer(
    cols = c(`Total BAM count`, `Single-End BAM count`, `Paired-End BAM count`),
    names_to = "BAM Type",
    values_to = "Count"
  )

plot_npd_counts <- ggplot(
  donor_data_df_npd,
  aes(x = Count, y = `BAM Type`)
) +
  geom_violin(
    scale = "width",
    fill = "black"
  ) +
  scale_x_log10() +
  labs(
    title = NULL,
    x = "Number per Donor (log₁₀)",
    y = NULL
  ) +
  theme_minimal()

donor_data_df_rt <- as_tibble(donor_data_df) %>%
  dplyr::select(Donor, `Total Mapped Reads`, `Total Unmapped Reads`) %>%
  pivot_longer(
    cols = c(`Total Mapped Reads`, `Total Unmapped Reads`),
    names_to = "Read Type",
    values_to = "Reads"
  )

plot_npd_reads <- ggplot(
  donor_data_df_rt,
  aes(x = Reads, y = `Read Type`)
) +
  geom_violin(
    scale = "width",
    fill = "Black"
  ) +
  scale_x_log10() +
  labs(
    title = NULL,
    x = "Number per Donor (log₁₀)",
    y = NULL
  ) +
  theme_minimal()