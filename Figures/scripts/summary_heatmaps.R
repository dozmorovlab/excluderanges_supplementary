library(rtracklayer)
library(GenomicRanges)
library(ggplot2)
library(reshape2)
library(viridis)
library(patchwork)
library(dplyr)
library(openxlsx)
library(gridExtra)

data_folder <- file.path("data", "tf_correlations")

setwd(file.path(here::here(), "Figures"))
source(file.path("scripts", "names.R"))

# Functions
make_heatmap <- function(
  in_matrix,
  title_ = NULL,
  scale_ = NULL,
  max_value_scale = Inf,
  color_trans = "identity" #"identity" "log10" "sqrt" "reverse"
) {
  df <- melt(as.matrix(in_matrix))
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
  df <- melt(as.matrix(in_matrix))
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
      #trans = "pseudo_log",
      breaks = sort(c(
        -seq(
          from = 0,
          to = max_value_scale,
          length.out = 4
        ),
        seq(
          from = 0,
          to = max_value_scale,
          length.out = 4
        )
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
get_tf_boxplots <- function(in_m, xlsx_path = NULL, y_lab = "Delta") {
  # Summary Stats
  sum_stats <- sapply(in_m, function(m) {
    m <- as.matrix(m)
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
get_list_boxplots <- function(in_m, csv_path = NULL, y_lab = "Delta") {
  # Summary Stats
  sum_stats <- sapply(in_m, function(m) {
    m <- as.matrix(m)
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
get_tf_boxplots_2 <- function(
  in_m,
  xlsx_path = NULL,
  y_lab = "Delta",
  tf_list = c("UBTF", "CTCF", "REST"),
  interesting_tf = c()
) {
  # Summary Stats
  sum_stats <- sapply(in_m, function(m) {
    m <- as.matrix(m)
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

    i_color <- "black"
    if (tf %in% interesting_tf) {
      i_color <- "red"
    }

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
      theme(plot.subtitle = element_text(color = i_color)) +
      coord_flip()
    return(out_plot)
  }, simplify = FALSE)

  # Add y label
  tf_plots[[length(tf_plots)]] <- tf_plots[[length(tf_plots)]] + labs(y = y_lab)

  # Return output
  return(tf_plots)
}

# Load corr
tf_corr_all <- sapply(
  list.files(file.path(".", data_folder)), function(tsv_path) {
    tsv_path <- file.path(".", data_folder, tsv_path)
    df <- read.delim(
      tsv_path,
      header = TRUE,
      sep = "\t",
      skip = 1,
      check.names = FALSE,
      row.names = 1
    )
    rownames(df) <- gsub(
      "'|-human_(sponge_)*Aligned.sortedByCoord.out.autosomes.bam",
      "",
      rownames(df)
    )
    colnames(df) <- gsub(
      "'|-human_(sponge_)*Aligned.sortedByCoord.out.autosomes.bam",
      "",
      colnames(df)
    )
    # reorder columns
    df <- df[, order(colnames(df))]
    # reorder rows
    df <- df[order(rownames(df)), ]
    return(df)
  },
  simplify = FALSE
)

# Separate into groups
tf_corr <- list(
  pearson = tf_corr_all[grep("pearson", names(tf_corr_all))],
  spearman = tf_corr_all[grep("spearman", names(tf_corr_all))]
)
tf_corr <- lapply(tf_corr, function(in_list) {
  list(
    with_outliers = in_list[
      grep("no_outliers", names(in_list), invert = TRUE)
    ],
    no_outliers = in_list[
      grep("no_outliers", names(in_list))
    ]
  )
})
tf_corr <- lapply(tf_corr, function(list_a) {
  lapply(list_a, function(in_list) {
    list(
      no_sponge = in_list[
        grep("sponge", names(in_list), invert = TRUE)
      ],
      with_sponge = in_list[
        grep("sponge", names(in_list))
      ]
    )
  })
})
tf_corr <- lapply(tf_corr, function(list_a) {
  lapply(list_a, function(list_b) {
    lapply(list_b, function(in_list) {
      names(in_list) <- gsub(
        "all\\.|\\.*sponge|\\.*no_outliers|\\.*pearson|\\.*spearman|\\.*matrix\\.tsv",
        "",
        names(in_list)
      )
      names(in_list) <- names_dict[names(in_list)]
      return(in_list)
    })
  })
})


# Get delta
names(tf_corr$pearson$with_outliers$with_sponge$`NA`)
tf_corr[[1]][[1]]
tf_corr_delta <- lapply(tf_corr, function(list_a) {
  lapply(list_a, function(list_b){
    og_hg38 <- list_b$no_sponge$`NA`
    sponge <- list_b$with_sponge$`NA` - og_hg38
    others <- lapply(list_b$no_sponge, function(df) {
      return(df - og_hg38)
    })
    return(c(others, list(Sponge = sponge)))
  })
})

# Choose dataset
matricies <- tf_corr_delta$pearson$with_outliers

max_matrix <- max(unlist(lapply(matricies, function(m) {
  max(abs(m), na.rm = TRUE)
})))

n <- 1
m <- matricies$`GitHub Blacklist`
clustered_order <- rownames(m)[hclust(
  dist(as.matrix(m)),
  method = "complete"
)$order]
plots_delta <- sapply(seq_along(matricies), function(n) {
  m <- matricies[[n]]
  row_order <- clustered_order
  current_name <- names(matricies)[[n]]
  title_ <- current_name
  out_plot <- make_heatmap_delta(
    m[row_order, row_order],
    title_ = NULL,
    scale_ = "Delta",
    max_value_scale = max_matrix
    #color_trans = "none"
  ) + labs(subtitle = title_)
  return(out_plot)
}, simplify = FALSE)
names(plots_delta) <- names(matricies)

# Order by overall mean
delta_order <- names(matricies)[
  order(unlist(sapply(matricies, function(m) {
    mean(as.matrix(m), na.rm = TRUE)
  }, simplify = FALSE)), decreasing = FALSE)
]
plots_delta <- plots_delta[delta_order]

# Exclude Original
plots_delta <- plots_delta[!is.na(names(plots_delta))]

# Summary Stats
m <- matricies[[1]]
sum_stats_delta <- sapply(
  matricies[!is.na(names(matricies))],
  function(m) {
    m <- as.matrix(m)
    # Get upper triangle
    upper_tri_values <- m[upper.tri(m)]
    out_stats <- c(
      summary(upper_tri_values),
      SD = sd(upper_tri_values)
    )
  }, simplify = FALSE
)

sum_stats_delta <- sum_stats_delta[
  delta_order[!is.na(delta_order)]
]
delta_summary_table <- as.data.frame(t(
  as.data.frame(sum_stats_delta, check.names = FALSE)
), check.names = FALSE)[c("Mean", "Median", "SD")]

# Format names


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
    title = "Upper Triangle of Transcription Factors",
    subtitle = NULL,
    x = "List",
    y = "Delta"
  ) +
  theme_minimal() +
  coord_flip()

#wrap_plots(plots) +
#  plot_layout(guides = "collect")

table_a <- tableGrob(
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
)

heatmaps_a <- wrap_plots(
  plots_delta,
  guides = "collect",
  widths = c(1, 1, 1),
  heights = c(1, 1, 1),
  ncol = 3,
  nrow = 3
) & theme(
  legend.position = "left",
  legend.title = element_text(size = 5),
  legend.text = element_text(angle = 0, size = 3),
  legend.key.width = unit(2, "pt"),
  legend.key.height = unit(75, "pt"),
  axis.text.x = element_text(size = 4),
  axis.text.y = element_text(size = 4),
  legend.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt")
)

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
  matricies[!is.na(names(matricies))],
  xlsx_path = "tf_change_forbes_CTCF_REST_UBTF.xlsx",
  y_lab = "Delta",
  tf_list = clustered_order[1:3]#c("UBTF", "REST", "CTCF")
)
heatmaps_b <- wrap_plots(heatmaps_b, ncol = 1)

interesting_tf <- c(
  "ATF2", "ELF1", "NFATC1", "NFIC", "PBX3",
  "RUNX3", "SIX5", "STAT5A", "ZEB1"
)

tf_boxplots_all <- get_tf_boxplots_2(
  matricies[!is.na(names(matricies))],
  xlsx_path = "tf_change_forbes_all_tf.xlsx",
  y_lab = "Delta",
  tf_list = c(names(matricies[[1]])),
  interesting_tf = interesting_tf
)
tf_boxplots_all <- wrap_plots(tf_boxplots_all)

extract_legend <- function(p) {
  g <- ggplotGrob(p)
  legend <- g$grobs[which(sapply(g$grobs, function(x) x$name) == "guide-box")]
  return(legend)
}

# Format
legend <- extract_legend(heatmaps_a[[2]])
legend_plot <- wrap_plots(
  ggplot() + legend,
  ncol = 1
)

heatmaps_a <- heatmaps_a & theme(
  legend.position = "none"
)
heatmaps_a[[1]] <- (
  heatmaps_a[[1]] +
    labs(tag = "A") +
    theme(plot.tag.position = c(0.00, 0.98))
) #+
  #plot_layout(heights = c(3, 1))
heatmaps_b[[1]] <- heatmaps_b[[1]] +
  labs(tag = "B") +
  theme(plot.tag.position = c(0.00, 0.98))
chip_final_plot <- (heatmaps_a | legend_plot| heatmaps_b) +
  plot_layout(
    ncol = 3,
    widths = c(2, 0, 1)
  )

#ggsave(
#  filename = file.path("images", "delta_heatmaps_pearson.png"),
#  plot = chip_final_plot,
#  width = 6,
#  height = 4,
#  dpi = 300,
#  scale = 2
#)

#ggsave(
#  filename = file.path("images", "tf_boxplots_all.png"),
#  plot = tf_boxplots_all,
#  width = 6,
#  height = 4,
#  dpi = 300,
#  scale = 3
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