library(openxlsx)
library(dplyr)
library(ggplot2)
library(patchwork)
library(rtracklayer)
library(GenomicRanges)
library(biomaRt)
library(jsonlite)

bed_dir <- file.path("data", "figure_2_exp")

# https://www.gencodegenes.org/human/
gff3_path <- file.path("data", "gencode.v46.annotation.gff3.gz")

include_chr <- c(
  "chr1", "chr2", "chr3", "chr4", "chr5",
  "chr6", "chr7", "chr8", "chr9", "chr10",
  "chr11", "chr12", "chr13", "chr14", "chr15",
  "chr16", "chr17", "chr18", "chr19", "chr20",
  "chr21", "chr22", "chrX", "chrY"
)

# Load BED files
bed_list <- list.files(bed_dir)
gr_list <- sapply(bed_list, function(b_l) {
  import(file.path(bed_dir, b_l), format = "BED")
}, simplify = FALSE)
compressed_gr_list <- GRangesList(gr_list)

# Load GFF3
gff3_gr <- import(gff3_path, "GFF3")

# Filter
gr_list <- sapply(gr_list, function(g_l) {
  keepSeqlevels(
    g_l,
    include_chr[include_chr %in% seqnames(g_l)],
    pruning.mode = "tidy"
  )
})
gff3_gr <- keepSeqlevels(
  gff3_gr,
  include_chr[include_chr %in% seqnames(gff3_gr)],
  pruning.mode = "tidy"
)

genes <- gff3_gr[gff3_gr$type == "gene"]
gene_list <- split(genes, mcols(genes)$gene_id)
exons <- gff3_gr[gff3_gr$type == "exon"]
gene_exon_list <- split(exons, mcols(exons)$gene_id)
#colnames(mcols(gff3_gr))

# Get descriptions
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

gene_bm <- getBM(
  attributes = c(
    "ensembl_gene_id_version",
    "description",
    "interpro",
    "interpro_short_description",
    "interpro_description",
    "gene_biotype"
  ),
  filters = "ensembl_gene_id_version",
  values = genes$gene_id,
  mart = ensembl
)

gene_bm$description <- gsub(" \\[(.*?)\\]", "", gene_bm$description)

#write(toJSON(gene_bm, pretty = TRUE), "./gene_bm.json")

# python call
#system("python ./bio.py")

#linkage_df <- as.data.frame(fromJSON(readLines(
#  "./gene_bm_linkage_matrix_all.json",
#  warn = FALSE
#)))

#linkage_merge <- linkage_df[, c("V1", "V2")] + 1
#linkage_merge[
#  linkage_merge <= nrow(linkage_df) + 1
#] <- -linkage_merge[linkage_merge <= nrow(linkage_df) + 1]
#linkage_merge[
#  linkage_merge > 0
#] <- linkage_merge[linkage_merge > 0] -
#  (length(linkage_merge[[1]]) + 1)

#linkage_height <- linkage_df$V3

#hc <- list()
#hc$merge <- as.matrix(linkage_merge)
#hc$height <- linkage_height
#hc$order <- 1:(nrow(linkage_df) + 1)
#hc$labels <- gene_bm$ensembl_gene_id_version
#class(hc) <- "hclust"
##plot(hc)

#cutree_hc <- cutree(hc, h = 1)
#cutree_hc_table <- table(cutree_hc)
#cutree_hc_table <- cutree_hc_table[
#  order(cutree_hc_table, decreasing = TRUE)
#]
#cutree_hc_dict <- split(names(cutree_hc), cutree_hc)
#clusters_all <- cutree_hc_dict[names(cutree_hc_table)]

# Filter to include only protein-coding
#protein_coding_accessions <- gene_bm$ensembl_gene_id_version[
#  gene_bm$gene_biotype == "protein_coding"
#]
#clusters_protein_coding <- sapply(clusters_all, function(c) {
#  out_v <- c[c %in% protein_coding_accessions]
#  if(length(out_v) > 0) {
#    return(out_v)
#  }
#}, simplify = FALSE)
#clusters_protein_coding <- clusters_protein_coding[
#  sapply(clusters_protein_coding, function(c) length(c) > 0)
#]

# Sort decreasing
#clusters_protein_coding <- clusters_protein_coding[
#  order(
#    sapply(clusters_protein_coding, function(c) length(c)),
#    decreasing = TRUE
#  )
#]

# Get descriptions
#clusters_protein_coding_desc <- sapply(clusters_protein_coding, function(c) {
#  gene_bm$description[gene_bm$ensembl_gene_id_version %in% c]
#})

# Flatten
#clusters_protein_coding_flat <- do.call(c, lapply(
#  names(clusters_protein_coding),
#  function(n) {
#    setNames(
#      rep(n, length(clusters_protein_coding[[n]])),
#      clusters_protein_coding[[n]]
#    )
#  }
#))

gene_desc <- setNames(gene_bm$description, gene_bm$ensembl_gene_id_version)
gene_interpro <- setNames(gene_bm$interpro, gene_bm$ensembl_gene_id_version)

interpro_df <- unique(
  gene_bm[c(
    "interpro",
    "interpro_short_description",
    "interpro_description"
  )]
)
rownames(interpro_df) <- interpro_df$interpro

exon_widths <- sum(width(gene_exon_list))
main_df <- data.frame(
  Name = genes$gene_name,
  ID = genes$ID,
  Description = gene_desc[genes$gene_id],
  `Gene Width` = sum(width(gene_list))[genes$gene_id],
  `Sum of Exon Widths` = exon_widths[genes$gene_id],
  `Number of Exons` =  elementNROWS(gene_exon_list)[genes$gene_id],
  Class = genes$gene_type,
  `InterPro Accession` = gene_interpro[genes$gene_id],
  #`Description Cluster` = clusters_protein_coding_flat[genes$gene_id],
  check.names = FALSE
)

write.csv(main_df[main_df$Class == "protein_coding", ], "main_df.csv")

gr <- gr_list[[1]]
main_results_list <- lapply(gr_list, function(gr) {

  # Number of exclusion regions affecting a gene
  genes_pairs <- findOverlapPairs(
    query = genes,
    subject = gr
  )
  i_genes <- pintersect(genes_pairs)
  i_genes_split <- split(i_genes, mcols(i_genes)$gene_id)
  i_genes_widths <- setNames(sum(width(i_genes_split)), names(i_genes_split))
  n_ex_regions_per_gene <- table(i_genes$gene_id)
  no_hit_genes <- setNames(
    rep(0, length(setdiff(genes$gene_id, names(n_ex_regions_per_gene)))),
    setdiff(genes$gene_id, names(n_ex_regions_per_gene))
  )
  n_ex_regions_per_gene <- c(
    n_ex_regions_per_gene,
    no_hit_genes
  )

  # Number of exons overlapped with an exclusion region by at least 1bp
  all_pairs <- findOverlapPairs(
    query = gff3_gr,
    subject = gr
  )
  i_all <- pintersect(all_pairs)
  i_all_exons <- i_all[i_all$type == "exon"]
  i_all_exons_per_gene <- split(i_all_exons, mcols(i_all_exons)$gene_id)
  i_exons_widths <- setNames(sum(width(i_all_exons_per_gene)), names(i_all_exons_per_gene))
  n_exons_overlap_ex_per_gene <- elementNROWS(width(i_all_exons_per_gene))
  no_hit_n_exons_overlap <- setNames(
    rep(0, length(setdiff(genes$gene_id, names(n_exons_overlap_ex_per_gene)))),
    setdiff(genes$gene_id, names(n_exons_overlap_ex_per_gene))
  )
  n_exons_overlap_ex <- c(
    n_exons_overlap_ex_per_gene,
    no_hit_n_exons_overlap
  )

  # Width of bases covered by exclusion regions
  no_hit_w_bases_ex <- setNames(
    rep(0, length(setdiff(genes$gene_id, names(i_genes_widths)))),
    setdiff(genes$gene_id, names(i_genes_widths))
  )
  w_bases_ex <- c(i_genes_widths, no_hit_w_bases_ex)
  #sum(duplicated(names(i_genes_widths)))

  # Proportion of the gene length affected by exclusion regions
  prop_gene_ex <- w_bases_ex /
    setNames(width(genes), genes$gene_id)[names(w_bases_ex)]

  # Sum of widths of exons affected by exclusion regions
  no_hit_w_bases_exons_ex <- setNames(
    rep(0, length(setdiff(genes$gene_id, names(i_exons_widths)))),
    setdiff(genes$gene_id, names(i_exons_widths))
  )
  w_bases_exons_ex <- c(i_exons_widths, no_hit_w_bases_exons_ex)
  prop_gene_exons_ex <- w_bases_exons_ex /
    exon_widths[names(w_bases_exons_ex)]
  #sum(duplicated(names(i_exons_widths)))

  # Add fields to out_df
  out_df <- main_df[main_df$ID, ]
  out_df$`Number of Exclusion Regions Overlapping Gene` <-
    n_ex_regions_per_gene[main_df$ID]
  out_df$`Number of Exons Overlapping an Exclusion Region` <-
    n_exons_overlap_ex[main_df$ID]
  out_df$`Width of Gene Covered by Exclusion Regions` <-
    w_bases_ex[main_df$ID]
  out_df$`Proportion of Gene Width Covered by Exclusion Regions` <-
    prop_gene_ex[main_df$ID]
  out_df$`Sum of Exon Widths Covered by Exclusion Regions` <-
    w_bases_exons_ex[main_df$ID]
  out_df$`Proportion of Exon Widths Covered by Exclusion Regions` <-
    prop_gene_exons_ex[main_df$ID]
  return(out_df)
})

names(main_results_list) <- unname(names_dict[
  sub(paste0(".bed", "$"), "", names(main_results_list))
])

interpro_results_list <- lapply(gr_list, function(gr) {

  # Number of exclusion regions affecting a gene
  genes_pairs <- findOverlapPairs(
    query = genes,
    subject = gr
  )
  i_genes <- pintersect(genes_pairs)

  # Interpro
  affected_genes <- unique(i_genes$gene_id)
  interpro_counts <- table(gene_interpro[affected_genes])
  interpro_counts_sorted <- sort(interpro_counts, decreasing = TRUE)
  interpro_counts_sorted <- interpro_counts_sorted[
    names(interpro_counts_sorted) != ""
  ]
  out_df <- data.frame(
    `Number of Genes Affected` = c(interpro_counts_sorted),
    `InterPro Accession` = names(interpro_counts_sorted),
    `Short Description` = interpro_df[names(interpro_counts_sorted), ]$interpro_short_description,
    `Full Description` = interpro_df[names(interpro_counts_sorted), ]$interpro_description,
    check.names = FALSE
  )

  return(out_df)
})

names(interpro_results_list) <- unname(names_dict[
  sub(paste0(".bed", "$"), "", names(interpro_results_list))
])

# Create a new workbook
main_wb <- createWorkbook()
interpro_wb <- createWorkbook()

# Loop through the list and add each data frame to a new sheet
for (name in names(main_results_list)) {
  addWorksheet(main_wb, substr(name, 1, 31))
  writeData(main_wb, sheet = substr(name, 1, 31), main_results_list[[name]])
}
for (name in names(interpro_results_list)) {
  addWorksheet(interpro_wb, substr(name, 1, 31))
  writeData(
    interpro_wb,
    sheet = substr(name, 1, 31),
    interpro_results_list[[name]]
  )
}

# Save the workbook
saveWorkbook(
  main_wb,
  file = "biological_characterization.xlsx",
  overwrite = TRUE
)
saveWorkbook(
  interpro_wb,
  file = "interpro_analysis.xlsx",
  overwrite = TRUE
)

# Subset based on ex overlap
class_list_main <- lapply(main_results_list, function(df) {
  df[df$`Number of Exclusion Regions Overlapping Gene` > 0, ]
})

# Counts
## class
count_class <- function(df, name_) {
  df %>%
    count(Class) %>%
    mutate(BED = name_)
}
class_counts_main_df <- bind_rows(mapply(
  count_class,
  class_list_main,
  names(class_list_main),
  SIMPLIFY = FALSE
))
class_levels_main <- sapply(unique(class_counts_main_df$Class), function(c) {
  max(class_counts_main_df$n[class_counts_main_df$Class == c])
}, simplify = FALSE)
class_levels_main <- unlist(class_levels_main)
class_levels_main <- class_levels_main[
  order(class_levels_main, decreasing = TRUE)
]
class_counts_main_df$Class <- factor(
  class_counts_main_df$Class,
  levels = names(class_levels_main)
)

## InterPro
count_ip <- function(df, name_) {
  df %>%
    count(`InterPro Accession`) %>%
    mutate(BED = name_)
}
class_counts_ip_df <- bind_rows(mapply(
  count_ip,
  class_list_main,
  names(class_list_main),
  SIMPLIFY = FALSE
))
ip_list <- setNames(interpro_df$interpro_description, interpro_df$interpro)
class_counts_ip_df$`Class` <- ip_list[
  class_counts_ip_df$`InterPro Accession`
]

# Remove NA
class_counts_ip_df <- class_counts_ip_df[! is.na(class_counts_ip_df$Class), ]

order_by <- "GitHub Blacklist"
c <- class_counts_ip_df$Class[[1]]
class_levels_ip <- sapply(unique(class_counts_ip_df$Class), function(c) {
  max(
    class_counts_ip_df$n[
      (class_counts_ip_df$Class == c) &
        (class_counts_ip_df$BED == order_by)
    ]
  )
}, simplify = FALSE)
class_levels_ip <- unlist(class_levels_ip)
class_levels_ip <- class_levels_ip[
  order(class_levels_ip, decreasing = TRUE)
]
class_counts_ip_df$Class <- factor(
  class_counts_ip_df$Class,
  levels = names(class_levels_ip)
)

# Widths
## class
width_class <- function(df, name_) {
  df %>%
    group_by(Class) %>%
    summarise(
      n = sum(
        `Width of Gene Covered by Exclusion Regions`,
        na.rm = TRUE
      )
    ) %>%
    mutate(BED = name_)
}
class_widths_main_df <- bind_rows(mapply(
  width_class,
  class_list_main,
  names(class_list_main),
  SIMPLIFY = FALSE
))
class_widths_main_df$Class <- factor(
  class_widths_main_df$Class,
  levels = names(class_levels_main)
)

## InterPro
width_ip <- function(df, name_) {
  df %>%
    group_by(`InterPro Accession`) %>%
    summarise(
      n = sum(
        `Width of Gene Covered by Exclusion Regions`,
        na.rm = TRUE
      )
    ) %>%
    mutate(BED = name_)
}
class_widths_ip_df <- bind_rows(mapply(
  width_ip,
  class_list_main,
  names(class_list_main),
  SIMPLIFY = FALSE
))
class_widths_ip_df$`Class` <- ip_list[
  class_widths_ip_df$`InterPro Accession`
]
class_widths_ip_df <- class_widths_ip_df[! is.na(class_widths_ip_df$Class), ]

class_widths_ip_df$Class <- factor(
  class_widths_ip_df$Class,
  levels = names(class_levels_ip)
)

# Order by width
## Widths
class_widths_ip_by_width_df <- bind_rows(mapply(
  width_ip,
  class_list_main,
  names(class_list_main),
  SIMPLIFY = FALSE
))
class_widths_ip_by_width_df$`Class` <- ip_list[
  class_widths_ip_by_width_df$`InterPro Accession`
]
class_widths_ip_by_width_df <- class_widths_ip_by_width_df[! is.na(class_widths_ip_by_width_df$Class), ]

order_by <- "GitHub Blacklist"
class_levels_ip <- sapply(unique(class_widths_ip_by_width_df$Class), function(c) {
  max(
    class_widths_ip_by_width_df$n[
      (class_widths_ip_by_width_df$Class == c) &
        (class_widths_ip_by_width_df$BED == order_by)
    ]
  )
}, simplify = FALSE)
class_levels_ip <- unlist(class_levels_ip)
class_levels_ip <- class_levels_ip[
  order(class_levels_ip, decreasing = TRUE)
]
class_widths_ip_by_width_df$Class <- factor(
  class_widths_ip_by_width_df$Class,
  levels = names(class_levels_ip)
)

## Counts
class_counts_ip_by_width_df <- bind_rows(mapply(
  count_ip,
  class_list_main,
  names(class_list_main),
  SIMPLIFY = FALSE
))
class_counts_ip_by_width_df$`Class` <- ip_list[
  class_counts_ip_by_width_df$`InterPro Accession`
]

class_counts_ip_by_width_df <- class_counts_ip_by_width_df[! is.na(class_counts_ip_by_width_df$Class), ]

order_by <- "GitHub Blacklist"
c <- class_counts_ip_by_width_df$Class[[1]]

class_counts_ip_by_width_df$Class <- factor(
  class_counts_ip_by_width_df$Class,
  levels = names(class_levels_ip)
)

# Reorder
list_order <- c(
  "GitHub Blacklist",
  "Generated Blacklist",
  "Kundaje Unified",
  "High Signal",
  "Low Mappability",
  "HS + LM",
  "Nordin CUT&RUN",
  "GreyListChIP"
)
class_counts_main_df$BED <- factor(
  class_counts_main_df$BED,
  levels = list_order
)
class_widths_main_df$BED <- factor(
  class_widths_main_df$BED,
  levels = list_order
)

df <- class_counts_main_df
get_bars <- function(
  df,
  y_lab = NULL,
  cap_at = 2,
  other_lab = "Other",
  show_other = TRUE
) {

  if (length(unique(df$Class)) > cap_at) {
    # Get the top 'cap_at' unique classes
    top_classes <- as.character(levels(df$Class)[1:cap_at])

    # Update the Class column: keep top classes, others become "Other"
    df <- df %>%
      mutate(
        Class = ifelse(
          Class %in% top_classes,
          as.character(Class),
          other_lab
        )
      ) %>%
      group_by(Class, BED) %>%
      summarise(n = sum(n)) %>%
      ungroup()
    if (!show_other) {
      df <- df[df$Class != other_lab, ]
    }
  }

  df$Class <- gsub("_", " ", df$Class)
  class_levels <- gsub("_", " ", c(top_classes, other_lab))
  # Convert Class column to a factor with the specified levels
  df$Class <- factor(df$Class, levels = class_levels)

  out_plot <- ggplot(
    df,
    aes(
      x = factor(Class, levels = rev(levels(factor(Class)))),
      y = n,
      fill = BED,
      group = factor(BED, levels = rev(levels(factor(BED))))
    )
  ) +
    geom_bar(
      stat = "identity",
      position = position_dodge()
    ) +
    labs(
      x = NULL,
      y = y_lab,
      fill = NULL
    ) +
    theme_minimal() +
    theme(
      strip.text = element_blank(),
      strip.background = element_blank(),
      legend.position = "bottom"
    ) +
    coord_flip()

  return(out_plot)
}

class_counts_main_plot <- get_bars(
  class_counts_main_df,
  "Number of\nAffected Genes",
  2,
  "other transcripts"
) +
scale_y_continuous(labels = scales::label_number(scale = 1e-3, suffix = "k"))
class_widths_main_plot <- get_bars(
  class_widths_main_df,
  "Coverage of\nAffected Genes (Mbp)",
  2,
  "other transcripts"
) +
  scale_y_continuous(labels = scales::label_number(scale = 1e-6))

fig_bio_main <- class_counts_main_plot +
  (
    class_widths_main_plot +
      labs(x = NULL) +
      theme(
        axis.text.y = element_blank()
      )
  ) +
  plot_layout(guides = "collect") &
  theme(
    legend.position = "bottom",
    legend.direction = "horizontal"
    #axis.text.x = element_text(angle = 90)
  )
