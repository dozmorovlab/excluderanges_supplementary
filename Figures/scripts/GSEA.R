# Set up the environment
library(knitr)
library(pander)
panderOptions('table.split.table', Inf)
set.seed(1)
library(dplyr)
options(stringsAsFactors = FALSE)

# Libraries
library(readxl)
library(cowplot)
library(ggplot2)
library(stringr)
library(ComplexHeatmap)

# Settings

# Data
p_adj_cutoff   <- 0.05 # FDR cutoff
# degs_sheet     <- "1.BCM.15057EI_vs_BCM.15057" # Which worksheet contains differentially expressed genes
fileNameIn1 <- file.path("data", "biological_characterization_GSEA.xlsx")
fileNameOut1 <- file.path("data", "biological_characterization_GSEA.png")

msigdb_all     <- TRUE # Use all MSigDb categories (TRUE), or c("C2", "C5", "H") (FALSE)

max_GO_to_plot = 12   # Maximum number of GOs to plot
max_enrichment_length <- 50 # Maximum length of enrichment descriptions

# How many top enrichments to plot
ntop <- 50
# How long enrichment labels to plot
ntext <- 60
# Color palette for the heatmap, https://www.nceas.ucsb.edu/~frazier/RSpatialGuides/colorPaletteCheatsheet.pdf
col3 <- colorRampPalette(c('grey', 'red'))(20)
# col3 <- colorRampPalette(c('blue', 'gray', 'yellow'))(20)
# col3 <- colorRampPalette(c('green', 'black', 'red'))(20)
# col3 <- colorRamps::green2red(n = 20)
col3 <- RColorBrewer::brewer.pal(20, "RdBu")

# Sheets depending on what analyses were run
if (msigdb_all) {
  all_sheets <- excel_sheets(fileNameIn1)
} else {
  all_sheets <- c("Enrich.KEGG", "GSEA.KEGG", "Enrich.C2", "GSEA.C2", "Enrich.C5", "GSEA.C5", "Enrich.H", "GSEA.H")
}

# list to keep data from all files
res_list <- list()

for (sheet in all_sheets) {
  print("========================================================")
  print(paste(sheet, "enrichments"))
  res <- read_xlsx(fileNameIn1, sheet = sheet)
  res_list <- c(res_list, list(res[, c("Term", "P.value")]))
}
# Merge all results
res_all <- res_list %>% purrr::reduce(full_join, by = c("Term"))


# Make a matrix
mtx_to_plot <- res_all %>% dplyr::select(starts_with("P.value")) %>% as.data.frame()
rownames(mtx_to_plot) <- res_all$Term %>% str_trunc(., width = ntext) %>% make.unique()
colnames(mtx_to_plot) <- all_sheets
# Replace NAs with zeros
mtx_to_plot[is.na(mtx_to_plot)] <- 1
mtx_to_plot <- -log10(mtx_to_plot)
# Sort by average absoulte NES
mtx_to_plot <- mtx_to_plot[order(rowMeans(mtx_to_plot), decreasing = TRUE), ]
mtx_to_plot <- mtx_to_plot[1:ntop, ]
# Rename
rownames(mtx_to_plot)[rownames(mtx_to_plot) == "Non-alcoholic fatty liver disease (NAFLD)"] <- "Non-alcoholic fatty liver disease"
# Top most significant 
# Heatmap(mtx_to_plot, show_row_dend = FALSE, show_column_dend = FALSE,
#         column_title = paste("Top highest", ntop, "enrichments"))

# Load necessary libraries
library(ggplot2)
library(reshape2)

# Your matrix (mtx_to_plot) - make sure it is a data frame
#mtx_to_plot <- as.matrix(mtx_to_plot[1:10, ])

# Convert the matrix into a long format
mtx_long <- melt(as.matrix(mtx_to_plot[1:12, ]))

# Rename columns for easier use in ggplot
colnames(mtx_long) <- c("Pathway", "Condition", "Value")
unique(mtx_long$Pathway)
mtx_long$Pathway <- factor(mtx_long$Pathway, levels = rev(mtx_long$Pathway[1:12]))
mtx_long$Condition <- factor(mtx_long$Condition, levels =  c("GitHub Blacklist", "Generated Blacklist", "Kundaje Unified", "High Signal", "Low Mappability", "HS + LM", "GreyListChIP", "Nordin CUT&RUN"))

#c("GitHub Blacklist", "Generated Blacklist", "Kundaje Unified", "High Signal", "Low Mappability", "HS + LM", "GreyListChIP", "Nordin CUT&RUN")

#c("hg38 GitHub List", "hg38 Generated List", "hg38 Kundaje Unified", "High Signal + Centromere", "Low Mappability", "HS + LM + CM", "GreyListChIP STAR 1k", "hg38 Nordin CUT&RUN")

# Create the dot plot
dot_plot <- ggplot(mtx_long, aes(x = Condition, y = Pathway)) +
  geom_point(aes(size = Value, color = Value)) +
  scale_size(range = c(2, 10)) +
  scale_color_gradientn(colors = c("white", "red"),
                        values = c(0, 0.1, 1)) +  # Pushes the gradient more towards red
  theme_minimal() +
  labs(
    title = NULL,
    x = NULL,
    y = NULL,
    size = "-log₁₀(FDR)",
    color = NULL
  ) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))
