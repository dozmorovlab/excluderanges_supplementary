library(rtracklayer)
library(GenomicRanges)
library(ggplot2)
library(reshape2)
library(viridis)
library(patchwork)
library(dplyr)
library(openxlsx)

data_folder <- file.path("data", "peaks_data")

# Functions
make_forbes_dist <- function(in_matrix, original_values = TRUE) {
  epsilon <- 1e-10 # Define a small constant epsilon to avoid division by zero

  # Convert Forbes coefficients to distances
  dist_matrix <- 1 / (in_matrix + epsilon)

  # Set diagonal to 0 since distance from an element to itself is 0
  diag(dist_matrix) <- 0

  # Create a full distance matrix instead of dist object
  full_dist_matrix <- as.matrix(dist_matrix)

  # Perform hierarchical clustering
  hc <- hclust(dist(full_dist_matrix), method = "ward.D2")

  # Get the order of the rows
  ordered_rows <- hc$order

  # Create the ordered matrix based on the clustering
  if (original_values) {
    ordered_matrix <- in_matrix[ordered_rows, ordered_rows]
  } else {
    ordered_matrix <- full_dist_matrix[ordered_rows, ordered_rows]  # Use full distance matrix
  }

  return(ordered_matrix)  # Return the ordered matrix
}
make_heatmap <- function(
  in_matrix,
  title_ = NULL,
  scale_ = NULL,
  max_value_scale = Inf,
  color_trans = "identity" #"identity" "log10" "sqrt" "reverse"
) {
  df <- melt(in_matrix)
  # Reverse the levels of Var2
  df$Var2 <- factor(df$Var2, levels = rev(unique(df$Var2)))

  out_plot <- ggplot(
    df,
    aes(x = Var1, y = Var2, fill = value)
  ) +
    geom_tile() +
    scale_fill_viridis_c(
      option = "magma",
      limits = c(0, max_value_scale),
      trans = color_trans,
      na.value = "white"
    ) +
    labs(
      x = NULL,
      y = NULL,
      fill = scale_
    ) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    coord_fixed(ratio = 1)

  return(out_plot)
}
make_heatmap_delta <- function(
  in_matrix,
  title_ = NULL,
  scale_ = NULL,
  max_value_scale = 100,
  color_trans = "identity" #"identity" "log10" "sqrt" "reverse"
) {
  df <- melt(in_matrix)
  # Reverse the levels of Var2
  df$Var2 <- factor(df$Var2, levels = rev(unique(df$Var2)))

  out_plot <- ggplot(
    df,
    aes(x = Var1, y = Var2, fill = value)
  ) +
    geom_tile() +
    scale_fill_gradient2(
      low = "red",
      mid = "white",
      high = "blue",
      midpoint = 0,
      limits = c(-max_value_scale, max_value_scale),
      na.value = "white",
      trans = "pseudo_log",
      breaks = sort(c(
        -round(10 ^ seq(
          from = 0,
          to = log10(max_value_scale),
          length.out = 4
        )),
        round(10 ^ seq(
          from = 0,
          to = log10(max_value_scale),
          length.out = 4
        ))
      )),
      labels = scales::label_number()
    ) +
    labs(
      x = NULL,
      y = NULL,
      fill = scale_
    ) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    coord_fixed(ratio = 1)

  return(out_plot)
}

# Load peaks
peak_gr_list <- sapply(
  list.files(file.path(".", data_folder)), function(bed_path) {
    bed_path <- file.path(".", data_folder, bed_path)
    return(import(bed_path, format = "narrowPeak"))
  }
)

# Load exclusion lists
ex_gr_list <- make_excludeList(
  path_to_sets = file.path("data", "ex_data"),
  remove_gaps = FALSE,
  path_to_gaps = NULL
)

# Filter chr
include_chr <- c(
  "chr1", "chr2", "chr3", "chr4", "chr5",
  "chr6", "chr7", "chr8", "chr9", "chr10",
  "chr11", "chr12", "chr13", "chr14", "chr15",
  "chr16", "chr17", "chr18", "chr19", "chr20",
  "chr21", "chr22", "chrX", "chrY"
)
peak_gr_list <- lapply(peak_gr_list, function(gr) {
  gr <- keepSeqlevels(
    gr,
    include_chr[include_chr %in% seqnames(gr)],
    pruning.mode = "tidy"
  )
  return(gr)
})
ex_gr_list <- lapply(ex_gr_list, function(gr) {
  gr <- keepSeqlevels(
    gr,
    include_chr[include_chr %in% seqnames(gr)],
    pruning.mode = "tidy"
  )
  return(gr)
})

genome_size <- get_genome_size(
  genome_id = "hg38",
  only_autosomal = FALSE,
  include_chr = include_chr
)

# Exclude peaks
ex_peaks_list <- sapply(ex_gr_list, function(ex_gr) {
  sapply(peak_gr_list, function(p_gr) {
    GenomicRanges::setdiff(p_gr, ex_gr)
  }, simplify = FALSE)
}, simplify = FALSE)

all_list <- c(list(Original = peak_gr_list), ex_peaks_list)

#saveRDS(matricies, "matricies.rds")

if (file.exists("matricies.rds")) {
  matricies <- readRDS("matricies.rds")
} else {
  matricies <- sapply(all_list, function(gr_list) {
    peak_targets <- gsub(
      "-human_local_peaks.narrowPeak",
      "",
      names(all_list[[1]])
    )
    f_w_matrix <- overlap_matrix(
      LIST = gr_list,
      Names = peak_targets,
      fun = forbes_coeff_width,
      Vend = 1,
      genome_size = genome_size,
      only_autosomal = FALSE
    )
    return(f_w_matrix)
  }, simplify = FALSE)
}

# Get delta
matricies_delta <- sapply(matricies, function(m) {
  m - matricies[["Original"]]
}, simplify = FALSE)

# Get row order
#ordered_og <- make_forbes_dist(
#  matricies[["Original"]],
#  original_values = TRUE
#)
#row_order <- rownames(ordered_og)

max_matrix <- max(unlist(lapply(matricies, function(m) {
  max(abs(m), na.rm = TRUE)
})))

plots <- sapply(seq_along(matricies), function(n) {
  m <- matricies[[n]]
  m_means <- unlist(sapply(rownames(m), function(r) {
    summary(m[, r])[["Mean"]]
  }, simplify = FALSE))
  row_order <- names(m_means)[order(m_means)]
  out_plot <- make_heatmap(
    m[row_order, row_order],
    title_ = NULL,
    scale_ = "Forbes Width",
    max_value_scale = max_matrix
    #color_trans = "none"
  ) + labs(title = names(matricies)[[n]])
  return(out_plot)
}, simplify = FALSE)

max_delta <- max(unlist(lapply(matricies, function(m) {
  max(abs(m - matricies[["Original"]]), na.rm = TRUE)
})))
n <- 1
plots_delta <- sapply(seq_along(matricies), function(n) {
  m <- matricies[[n]]
  m <- m - matricies[["Original"]]
  m_means <- unlist(sapply(rownames(m), function(r) {
    summary(m[, r])[["Median"]]
  }, simplify = FALSE))
  row_order <- names(m_means)[order(m_means)]
  current_name <- names(matricies)[[n]]
  title_ <- names_dict[
    sub(
      paste0(".bed", "$"),
      "",
      current_name
    )
  ]
  out_plot <- make_heatmap_delta(
    m[row_order, row_order],
    title_ = NULL,
    scale_ = "Forbes Width\nChange",
    max_value_scale = max_delta
    #color_trans = "none"
  ) + labs(subtitle = title_)
  return(out_plot)
}, simplify = FALSE)
names(plots_delta) <- names(matricies)

# Order by overall mean
delta_order <- names(matricies)[
  order(unlist(sapply(matricies, function(m) {
    mean(m - matricies[["Original"]], na.rm = TRUE)
  }, simplify = FALSE)), decreasing = FALSE)
]
plots_delta <- plots_delta[delta_order]

# Exclude Original
plots_delta <- plots_delta[names(plots_delta) != "Original"]

# Summary Stats
sum_stats_delta <- sapply(
  matricies_delta[names(matricies_delta) != "Original"],
  function(m) {
    # Get upper triangle
    upper_tri_values <- m[upper.tri(m)]
    out_stats <- c(
      summary(upper_tri_values),
      SD = sd(upper_tri_values)
    )
  }, simplify = FALSE
)

sum_stats_delta <- sum_stats_delta[
  delta_order[delta_order != "Original"]
]
delta_summary_table <- as.data.frame(t(
  as.data.frame(sum_stats_delta, check.names = FALSE)
), check.names = FALSE)[c("Mean", "Median", "SD")]

# Format names
rownames(delta_summary_table) <- names_dict[
  sub(
    paste0(".bed", "$"),
    "",
    rownames(delta_summary_table)
  )
]
names(plots_delta) <- names_dict[
  sub(
    paste0(".bed", "$"),
    "",
    names(plots_delta)
  )
]
names(matricies) <- names_dict[
  sub(
    paste0(".bed", "$"),
    "",
    names(matricies)
  )
]
names(matricies_delta) <- names_dict[
  sub(
    paste0(".bed", "$"),
    "",
    names(matricies_delta)
  )
]

long_summary_delta_df <- do.call(
  rbind,
  lapply(names(sum_stats_delta), function(l) {

    l_df <- data.frame(
      List = l,
      Min = sum_stats_delta[[l]][["Min."]],
      Q1 = sum_stats_delta[[l]][["1st Qu."]],
      Median = sum_stats_delta[[l]][["Median"]],
      Mean = sum_stats_delta[[l]][["Mean"]],
      Q3 = sum_stats_delta[[l]][["3rd Qu."]],
      Max = sum_stats_delta[[l]][["Max."]],
      SD = sum_stats_delta[[l]][["SD"]],
      stringsAsFactors = FALSE
    )
    return(l_df)
  })
)
long_summary_delta_df$List <- factor(
  long_summary_delta_df$List,
  levels = long_summary_delta_df$List
)


sum_stats_delta_plot <- ggplot(
  long_summary_delta_df,
  aes(
    ymin = Min,
    lower = Q1,
    middle = Median,
    upper = Q3,
    ymax = Max,
    x = List,
    group = List
  )
) +
  # Horizontal reference line at y = 0
  geom_hline(
    yintercept = 0,
    color = "grey",
    linetype = "dashed"
  ) +
  # Box plot based on summary statistics
  geom_boxplot(stat = "identity") +
  # Mean points in red
  geom_point(
    aes(y = Mean),
    color = "red",
    size = 3,
    shape = 18
  ) +
  # Error bars for mean ± SD
  geom_errorbar(
    aes(ymin = Mean - SD, ymax = Mean + SD),
    width = 0.2,
    color = "red"
  ) +
  labs(
    title = "Summary Statistics",
    subtitle = "Upper Triangle of Transcription Factors",
    x = "List",
    y = "Forbes Width Change"
  ) +
  theme_minimal() +
  coord_flip()

#wrap_plots(plots) +
#  plot_layout(guides = "collect")

heatmaps_a <- wrap_plots(
  c(
    list(summary = tableGrob(
      round(delta_summary_table, 2),
      theme = ttheme_minimal(
        core = list(
          fg_params = list(fontface = "plain", fontsize = 8),
          padding = unit(c(5, 3), "pt")
        ),
        colhead = list(
          fg_params = list(fontface = "plain", fontsize = 8),
          padding = unit(c(5, 3), "pt")
        ),
        rowhead = list(
          fg_params = list(fontface = "plain", fontsize = 8),
          padding = unit(c(5, 3), "pt")
        )
      )
    )),
    plots_delta
  ),
  widths = c(1, 1, 1),
  heights = c(1, 1, 1),
  ncol = 3,
  nrow = 3
) & theme(
  legend.position = "bottom",
  legend.text = element_text(angle = 90, size = 3),
  legend.key.width = unit(0.4, "cm"),
  legend.key.height = unit(0.2, "cm"),
  axis.text.x = element_text(size = 4),
  axis.text.y = element_text(size = 4)
)

in_m <- matricies_delta
get_tf_boxplots <- function(in_m, xlsx_path = NULL, y_lab = "Forbes Width") {
  # Summary Stats
  sum_stats <- sapply(in_m, function(m) {
    sapply(rownames(m), function(tf) {
      q <- as.numeric(na.omit(m[tf, ]))
      query_summary <- c(
        summary(q),
        SD = sd(q)
      )
      return(query_summary)
    }, simplify = FALSE)

  }, simplify = FALSE)

  long_summary_df <- do.call(rbind, lapply(names(sum_stats), function(l) {

    l_df <- do.call(rbind, lapply(names(sum_stats[[l]]), function(tf) {

      tf_df <- data.frame(
        List = l,
        TF = tf,
        Min = sum_stats[[l]][[tf]][["Min."]],
        Q1 = sum_stats[[l]][[tf]][["1st Qu."]],
        Median = sum_stats[[l]][[tf]][["Median"]],
        Mean = sum_stats[[l]][[tf]][["Mean"]],
        Q3 = sum_stats[[l]][[tf]][["3rd Qu."]],
        Max = sum_stats[[l]][[tf]][["Max."]],
        SD = sum_stats[[l]][[tf]][["SD"]],
        stringsAsFactors = FALSE
      )

      return(tf_df)
    }))
    return(l_df)
  }))

  long_summary_df <- long_summary_df[
    order(long_summary_df$Mean, decreasing = FALSE),
  ]

  # Save xlsx
  if (!is.null(xlsx_path)) {
    wb <- createWorkbook()
    lapply(unique(long_summary_df$List), function(l) {
      addWorksheet(wb, substr(l, 1, 31))
      l_df <- long_summary_df[long_summary_df$List == l, ]
      writeData(wb, substr(l, 1, 31), l_df[, names(l_df) != "List"])
    })
    saveWorkbook(wb, xlsx_path, overwrite = TRUE)
  }

  # Plot TF
  tf_plots <- sapply(unique(long_summary_df$List), function(l) {
    df <- long_summary_df[long_summary_df$List == l, ]

    # Order by Mean
    df$TF <- factor(df$TF, levels = df$TF[order(df$Mean, decreasing = TRUE)])

    out_plot <- ggplot(
      df,
      aes(
        ymin = Min,
        lower = Q1,
        middle = Median,
        upper = Q3,
        ymax = Max,
        x = TF,
        group = TF
      )
    ) +
      geom_hline(
        yintercept = 0,
        color = "grey",
        linetype = "dashed"
      ) +
      geom_boxplot(stat = "identity") +
      # Add mean points
      geom_point(
        aes(y = Mean),
        color = "red",
        size = 3,
        shape = 18
      ) +
      labs(
        title = l,
        x = "Gene",
        y = y_lab
      ) +
      theme_minimal() +
      coord_flip()
    return(out_plot)
  }, simplify = FALSE)

  return(wrap_plots(tf_plots))
}

in_m <- matricies_delta
get_list_boxplots <- function(in_m, csv_path = NULL, y_lab = "Forbes Width") {
  # Summary Stats
  sum_stats <- sapply(in_m, function(m) {
    # Get upper triangle
    upper_tri_values <- m[upper.tri(m)]
    out_stats <- c(
      summary(upper_tri_values),
      SD = sd(upper_tri_values)
    )
  }, simplify = FALSE)

  long_summary_df <- do.call(rbind, lapply(names(sum_stats), function(l) {

    l_df <- data.frame(
      List = l,
      Min = sum_stats[[l]][["Min."]],
      Q1 = sum_stats[[l]][["1st Qu."]],
      Median = sum_stats[[l]][["Median"]],
      Mean = sum_stats[[l]][["Mean"]],
      Q3 = sum_stats[[l]][["3rd Qu."]],
      Max = sum_stats[[l]][["Max."]],
      SD = sum_stats[[l]][["SD"]],
      stringsAsFactors = FALSE
    )
    return(l_df)
  }))

  long_summary_df <- long_summary_df[
    order(long_summary_df$Mean, decreasing = FALSE),
  ]

  # Save csv
  if (!is.null(csv_path)) {
    write.csv(long_summary_df, csv_path)
  }

  # Plot

  df <- long_summary_df

  # Order by Mean
  df$List <- factor(
    df$List,
    levels = df$List[order(df$Mean, decreasing = TRUE)]
  )

  out_plot <- ggplot(
    df,
    aes(
      ymin = Min,
      lower = Q1,
      middle = Median,
      upper = Q3,
      ymax = Max,
      x = List,
      group = List
    )
  ) +
    geom_hline(
      yintercept = 0,
      color = "grey",
      linetype = "dashed"
    ) +
    geom_boxplot(stat = "identity") +
    # Add mean points
    geom_point(
      aes(y = Mean),
      color = "red",
      size = 3,
      shape = 18
    ) +
    labs(
      title = "Summary Statistics",
      subtitle = "Upper Triangle of Transcription Factors",
      x = "List",
      y = y_lab
    ) +
    theme_minimal() +
    coord_flip()
  return(out_plot)
}

in_m <- matricies_delta
get_tf_boxplots_2 <- function(
  in_m,
  xlsx_path = NULL,
  y_lab = "Forbes Width Change",
  tf_list = c("UBTF", "CTCF", "REST")
) {
  # Summary Stats
  sum_stats <- sapply(in_m, function(m) {
    sapply(rownames(m), function(tf) {
      q <- as.numeric(na.omit(m[tf, ]))
      query_summary <- c(
        summary(q),
        SD = sd(q)
      )
      return(query_summary)
    }, simplify = FALSE)

  }, simplify = FALSE)

  long_summary_df <- do.call(rbind, lapply(names(sum_stats), function(l) {

    l_df <- do.call(rbind, lapply(names(sum_stats[[l]]), function(tf) {

      tf_df <- data.frame(
        List = l,
        TF = tf,
        Min = sum_stats[[l]][[tf]][["Min."]],
        Q1 = sum_stats[[l]][[tf]][["1st Qu."]],
        Median = sum_stats[[l]][[tf]][["Median"]],
        Mean = sum_stats[[l]][[tf]][["Mean"]],
        Q3 = sum_stats[[l]][[tf]][["3rd Qu."]],
        Max = sum_stats[[l]][[tf]][["Max."]],
        SD = sum_stats[[l]][[tf]][["SD"]],
        stringsAsFactors = FALSE
      )

      return(tf_df)
    }))
    return(l_df)
  }))

  long_summary_df <- long_summary_df[
    order(long_summary_df$Mean, decreasing = FALSE),
  ]

  long_summary_df <- long_summary_df[long_summary_df$TF %in% tf_list, ]

  # Save xlsx
  if (!is.null(xlsx_path)) {
    wb <- createWorkbook()
    lapply(unique(long_summary_df$TF), function(tf) {
      addWorksheet(wb, substr(tf, 1, 31))
      tf_df <- long_summary_df[long_summary_df$TF == tf, ]
      writeData(wb, substr(tf, 1, 31), tf_df[, names(tf_df) != "TF"])
    })
    saveWorkbook(wb, xlsx_path, overwrite = TRUE)
  }

  # Plot TF
  tf_plots <- sapply(tf_list, function(tf) {
    df <- long_summary_df[long_summary_df$TF == tf, ]

    # Order by Mean
    df$List <- factor(df$List, levels = df$List[
      order(df$Mean, decreasing = TRUE)
    ])

    out_plot <- ggplot(
      df,
      aes(
        ymin = Min,
        lower = Q1,
        middle = Median,
        upper = Q3,
        ymax = Max,
        x = List,
        group = List
      )
    ) +
      geom_hline(
        yintercept = 0,
        color = "grey",
        linetype = "dashed"
      ) +
      geom_boxplot(stat = "identity") +
      # Add mean points
      geom_point(
        aes(y = Mean),
        color = "red",
        size = 3,
        shape = 18
      ) +
      labs(
        subtitle = tf,
        x = NULL,
        y = NULL
      ) +
      theme_minimal() +
      coord_flip()
    return(out_plot)
  }, simplify = FALSE)

  # Add y label
  tf_plots[[length(tf_plots)]] <- tf_plots[[length(tf_plots)]] + labs(y = y_lab)

  # Return output
  out_plots <- wrap_plots(tf_plots, ncol = 1)
  return(out_plots)
}


#get_tf_boxplots(
#  matricies,
#  xlsx_path = "tf_raw_forbes.xlsx",
#  y_lab = "Forbes Width"
#)
#get_tf_boxplots(
#  matricies_delta,
#  xlsx_path = "tf_change_forbes.xlsx",
#  y_lab = "Forbes Width Change"
#)

heatmaps_b <- get_tf_boxplots_2(
  matricies_delta[names(matricies_delta) != "Original"],
  xlsx_path = "tf_change_forbes_CTCF_REST_UBTF.xlsx",
  y_lab = "Forbes Width Change",
  tf_list = c("UBTF", "REST", "CTCF")
)

extract_legend <- function(p) {
  g <- ggplotGrob(p)
  legend <- g$grobs[which(sapply(g$grobs, function(x) x$name) == "guide-box")]
  return(legend)
}

# Format
legend <- extract_legend(heatmaps_a[[2]])
heatmaps_a <- heatmaps_a & theme(
  legend.position = "none"
)
heatmaps_a[[1]] <- (
  (
    heatmaps_a[[1]] +
      labs(tag = "A") +
      theme(plot.tag.position = c(0.00, 0.98))
  ) / legend
) +
  plot_layout(heights = c(3, 1))
heatmaps_b[[1]] <- heatmaps_b[[1]] +
  labs(tag = "B") +
  theme(plot.tag.position = c(0.00, 0.98))
chip_final_plot <- (heatmaps_a | heatmaps_b) + plot_layout(
  ncol = 2, widths = c(2, 1)
)

#ggsave(
#  filename = "delta_heatmaps.png",
#  plot = final_plot,
#  width = 6,
#  height = 4,
#  dpi = 300,
#  scale = 2
#)

#get_list_boxplots(
#  matricies_delta,
#  csv_path = "list_all_tf_delta_forbes.csv",
#  y_lab = "Forbes Width Change"
#)

#summary_df <- as.data.frame(do.call(rbind, sum_stats))

#make_datatable <- function(df) {
#  return(
#    DT::datatable(
#      df,
#      options = list(pageLength = nrow(df)),
#      rownames = TRUE
#    )
#  )
#}

#make_datatable(summary_df)