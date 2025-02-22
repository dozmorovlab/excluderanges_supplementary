library(ggplot2)
library(patchwork)
library(tidyverse)

setwd(file.path(here::here(), "Figures"))

lost_df <- read.csv(file.path("data", "WGS", "lost_counts.csv"))
gained_df <- read.csv(file.path("data", "WGS", "gained_counts.csv"))

lost_df <- lost_df[
  lost_df$File %in% c("lost_unique.vcf", "lost_shared.vcf"),
]
gained_df <- gained_df[
  gained_df$File %in% c("gained_unique.vcf", "gained_shared.vcf"),
]

# Add Category column to distinguish Lost vs. Gained
lost_df$Category <- "Lost"
gained_df$Category <- "Gained"

# Combine both dataframes
combined_df <- rbind(lost_df, gained_df)

# Reshape the dataframe into long format
counts_df <- reshape(
  combined_df,
  varying = c("SNPS", "INDELS"),
  v.names = "Count",
  timevar = "Variant",
  times = c("SNPS", "INDELS"),
  direction = "long"
)

# Rename Type column based on File names
counts_df$Type <- ifelse(grepl("unique", counts_df$File), "Unique", "Benchmark")

# Remove unnecessary columns
counts_df <- counts_df[, c("Category", "Variant", "Type", "Count")]

# Order
counts_df$Category <- factor(counts_df$Category, levels = c("Lost", "Gained"))
counts_df$Variant <- factor(counts_df$Variant, levels = c("SNPS", "INDELS"))

# Plot
filtered_plot <- ggplot(
  counts_df,
  aes(x = Category, y = Count, fill = Type)
) +
  geom_bar(stat = "identity", position = "stack") +
  geom_text(
    aes(label = Count),
    position = position_stack(vjust = 0.5), # Center within stacked bars
    size = 4
  ) +
  facet_wrap(~ Variant) +  # Separate SNPs and INDELs
  labs(
    title = NULL,
    x = NULL,
    y = "Count",
    fill = NULL
  ) +
  theme_minimal()
