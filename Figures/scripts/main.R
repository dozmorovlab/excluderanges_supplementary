library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(GenomicRanges)
library(ggplot2)
library(ggpattern)
library(openxlsx)
library(IRanges)
library(dplyr)
library(ggvenn)
library(patchwork)
library(ggrepel)
library(ggdendro)

# functions

# Figure export
export_fig <- function(
  fig, name,
  dir = file.path("figures"),
  x = 7,
  y = 5,
  dpi = 1200,
  scale = 2
) {
  ggsave(
    filename = file.path(dir, paste0(name, ".svg")),
    plot = fig,
    width = x,
    height = y,
    units = "in",
    dpi = dpi,
    scale = scale
  )
  ggsave(
    filename = file.path(dir, paste0(name, ".eps")),
    plot = fig,
    width = x,
    height = y,
    units = "in",
    dpi = dpi,
    scale = scale
  )
  ggsave(
    filename = file.path(dir, paste0(name, ".png")),
    plot = fig,
    width = x,
    height = y,
    units = "in",
    dpi = dpi,
    scale = scale,
    bg = "white"
  )
}

gr_overlaps <- function(gr_a, gr_b, io = TRUE) {

  a_inside_b <- subsetByOverlaps(
    gr_a,
    gr_b,
    type = "within"
  )

  a_any_b <- subsetByOverlaps(
    gr_a,
    gr_b,
    type = "any"
  )

  a_overlap_b <- subsetByOverlaps(
    a_any_b,
    a_inside_b,
    invert = TRUE
  )

  a_outside_b <- subsetByOverlaps(
    gr_a,
    gr_b,
    invert = TRUE
  )

  out_list <- list(
    inside = a_inside_b,
    trans = a_overlap_b,
    outside = a_outside_b,
    any = a_any_b
  )

  if (!io) {

    b_inside_a <- subsetByOverlaps(
      gr_b,
      gr_a,
      type = "within"
    )

    b_any_a <- subsetByOverlaps(
      gr_b,
      gr_a,
      type = "any"
    )

    b_overlap_a <- subsetByOverlaps(
      b_any_a,
      b_inside_a,
      invert = TRUE
    )

    b_outside_a <- subsetByOverlaps(
      gr_b,
      gr_a,
      invert = TRUE
    )

    tw_a <- IRanges::setdiff(gr_a, gr_b)
    tw_c <- IRanges::intersect(gr_a, gr_b)
    tw_b <- IRanges::setdiff(gr_b, gr_a)

    out_list <- list(
      inside = a_outside_b,
      trans = c(a_overlap_b, b_overlap_a),
      trans_a = a_any_b,
      trans_b = b_any_a,
      outside = b_outside_a,
      any = a_any_b,
      tw_a = tw_a,
      tw_c = tw_c,
      tw_b = tw_b
    )
  }

  return(out_list)
}

gr_overlaps_apply_2 <- function(gr_list, gr_b, io = TRUE) {
  # Apply the gr_overlaps function to each range in gr_list
  result_list <- lapply(gr_list, function(gr_in) gr_overlaps(gr_in, gr_b, io))

  # Base components for output list using sapply
  out_list <- list(
    inside = sapply(result_list, `[[`, "inside", simplify = FALSE),
    trans  = sapply(result_list, `[[`, "trans", simplify = FALSE),
    outside = sapply(result_list, `[[`, "outside", simplify = FALSE),
    any = sapply(result_list, `[[`, "any", simplify = FALSE)
  )

  # Add additional components if io is FALSE
  if (!io) {
    extra_components <- c("trans_a", "trans_b", "tw_a", "tw_b", "tw_c")

    for (comp in extra_components) {
      out_list[[comp]] <- sapply(result_list, `[[`, comp, simplify = FALSE)
    }
  }

  return(out_list)
}

gr_overlaps_apply <- function(gr_list, gr_b, io = TRUE) {
  inside_list <- list()
  outside_list <- list()
  trans_list <- list()
  trans_a_list <- list()
  trans_b_list <- list()
  tw_a_list <- list()
  tw_b_list <- list()
  tw_c_list <- list()

  # Apply the gr_overlaps function to each range in gr_all
  result_list <- lapply(gr_list, function(gr_in) gr_overlaps(gr_in, gr_b, io))

  # Extract the inside, trans, and outside lists from the result_list
  inside_list <- lapply(result_list, function(res) res$inside)
  trans_list <- lapply(result_list, function(res) res$trans)
  outside_list <- lapply(result_list, function(res) res$outside)
  any_list <- lapply(result_list, function(res) res$any)

  # Combine the lists
  out_list <- list(
    inside = unlist(inside_list, recursive = FALSE),
    trans = unlist(trans_list, recursive = FALSE),
    outside = unlist(outside_list, recursive = FALSE),
    any = unlist(any_list, recursive = FALSE)
  )

  if (!io) {
    trans_a_list <- lapply(result_list, function(res) res$trans_a)
    trans_b_list <- lapply(result_list, function(res) res$trans_b)
    tw_a_list <- lapply(result_list, function(res) res$tw_a)
    tw_b_list <- lapply(result_list, function(res) res$tw_b)
    tw_c_list <- lapply(result_list, function(res) res$tw_c)

    # Combine the lists
    out_list <- list(
      inside = unlist(inside_list, recursive = FALSE),
      trans = unlist(trans_list, recursive = FALSE),
      trans_a = unlist(trans_a_list, recursive = FALSE),
      trans_b = unlist(trans_b_list, recursive = FALSE),
      outside = unlist(outside_list, recursive = FALSE),
      any = unlist(any_list, recursive = FALSE),
      tw_a = unlist(tw_a_list, recursive = FALSE),
      tw_b = unlist(tw_b_list, recursive = FALSE),
      tw_c = unlist(tw_c_list, recursive = FALSE)
    )
  }

  return(out_list)
}

get_num_regions <- function(
  n_bams,
  get_sd = FALSE,
  unique_names_ = NULL
) {
  if (is.null(unique_names_)) {
    unique_names_ <- names(n_bams)
    names_ <- unique_names_
  }
  out_data <- sapply(
    unique_names_,
    function(n) {
      # Get all sets created with same nBams
      temp_set <- n_bams[which(names_ == n)]
      if (get_sd) {
        sd_ <- sapply(temp_set, length) |>
          sd() |>
          round(0)
      }
      total_length <- sapply(temp_set, length) |> sum()
      average_length <- (total_length / length(temp_set))
      if (get_sd) {
        return(paste(average_length), "+/-", sd_)
      } else {
        return(average_length)
      }
    }
  )
  return(out_data)
}

get_percent <- function(values, totals) {
  return(sprintf("%.1f%%", values / totals * 100))
}

get_c_w <- function(
  n_bams,
  get_sd = FALSE,
  unique_names_ = NULL
) {
  if (is.null(unique_names_)) {
    unique_names_ <- names(n_bams)
    names_ <- unique_names_
  }
  out_data <- sapply(
    unique_names_,
    function(n) {
      # Get all sets created with same nBams
      temp_set <- n_bams[which(names_ == n)]
      sd_ <- sapply(
        temp_set,
        function(set_) {
          set_ |>
            width() |>
            sum()
        }
      ) |>
        sd() |>
        round(0)
      total_width <- sapply(
        temp_set,
        function(set_) {
          set_ |>
            width() |>
            sum()
        }
      ) |>
        sum()
      average_width <- (total_width / length(temp_set)) |>
        round(0)
      if (get_sd) {
        return(paste(average_width), "+/-", sd_)
      } else {
        return(average_width)
      }
    }
  )
  return(out_data)
}

get_c_x <- function(
  n_bams,
  func,
  get_sd = FALSE,
  unique_names_ = NULL
) {
  if (is.null(unique_names_)) {
    unique_names_ <- names(n_bams)
    names_ <- unique_names_
  }
  out_data <- sapply(
    unique_names_,
    function(n) {
      # Get all sets created with same nBams
      temp_set <- n_bams[which(names_ == n)]
      sd_ <- sapply(
        temp_set,
        function(set_) {
          set_ |>
            width() |>
            func()
        }
      ) |>
        sd() |>
        round(0)
      total_width <- sapply(
        temp_set,
        function(set_) {
          set_ |>
            width() |>
            func()
        }
      ) |>
        func()
      average_width <- (total_width / length(temp_set)) |>
        round(0)
      if (get_sd) {
        return(paste(average_width), "+/-", sd_)
      } else {
        return(average_width)
      }
    }
  )
  return(out_data)
}

get_t_w <- function(
  n_bams,
  gap,
  get_sd = FALSE,
  unique_names_ = NULL
) {
  if (is.null(unique_names_)) {
    unique_names_ <- names(n_bams)
    names_ <- unique_names_
  }
  out_data <- sapply(
    unique_names_,
    function(n) {
      # Get all sets created with same nBams
      temp_set <- n_bams[which(names_ == n)]
      total_width <- sapply(
        temp_set,
        function(set_) {
          return(
            sum(
              width(GenomicRanges::intersect(
                set_,
                gap,
                ignore.strand = TRUE
              ))
            )
          )
        }
      ) |>
        sum()
      #average_width <- (total_width / length(temp_set)) |>
      #  round(0)
      return(total_width)
    }
  )
  return(out_data)
}

make_boxplots <- function(data) {
  # Convert character columns to numeric
  data <- data %>%
    mutate(
      across(
        c(Min, Q1, Median, Mean, Q3, Max), ~as.numeric(gsub(",", "", .))
      )
    )

  out_plot <- ggplot(data, aes(x = reorder(Exclusion.set, Mean), y = Median)) +
    geom_boxplot(
      aes(
        y = Median,
        ymin = Min,
        lower = Q1,
        middle = Median,
        upper = Q3,
        ymax = Max
      ),
      position = position_dodge(width = 0.8),
      stat = "identity"
    ) +
    stat_summary(
      aes(y = Mean),
      fun = mean,
      geom = "point",
      position = position_dodge(width = 0.8),
      color = "red",
      size = 3,
      show.legend = TRUE
    ) +  # Add mean point
    scale_y_log10() +  # Set the y-axis to log scale
    labs(title = NULL, y = "log10", x = "Exclusion Set") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1)
    ) +  # Rotates x-axis labels for better visibility
    scale_color_manual(
      name = "Summary",
      values = c("red"),
      labels = c("Mean")
    )  # Legend for the mean point

  return(out_plot)
}

make_datatable <- function(df) {
  return(
    DT::datatable(
      df,
      options = list(pageLength = nrow(df)),
      rownames = FALSE
    )
  )
}
# Get number of grA regions intersecting grB regions
intersection_count <- function(grA, grB) {
  intersect_ <- subsetByOverlaps(grA, grB) |> length()
  return(intersect_)
}

get_widths <- function(
  n_bams,
  get_sd = FALSE,
  unique_names_ = NULL
) {
  if (is.null(unique_names_)) {
    unique_names_ <- names(n_bams)
    names_ <- unique_names_
  }
  out_data <- sapply(
    unique_names_,
    function(n) {
      # Get all sets created with same nBams
      temp_set <- n_bams[which(names_ == n)]
      if (get_sd) {
        sd_ <- sapply(
          temp_set,
          function(set_) {
            set_ |>
              width() |>
              sum()
          }
        ) |>
          sd() |>
          round(0)
      }
      total_width <- sapply(
        temp_set,
        function(set_) {
          set_ |>
            width() |>
            sum()
        }
      ) |>
        sum()
      average_width <- (total_width / length(temp_set)) |>
        round(0)
      if (get_sd) {
        return(paste(average_width), "+/-", sd_)
      } else {
        return(average_width)
      }
    }
  )
  return(out_data)
}

# Plots

get_plot_1a <- function(gaps_df, colors_ = list_colors_fig_1) {
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
    levels = c("Gaps", setdiff(names(colors_), "Gaps"))
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
      values = colors_,
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

make_bars_1c <- function(
  in_list,
  mode = "counts",
  colors_ = list_colors_fig_1,
  focus_on = NULL
) {

  in_combos <- combn(c(names(in_list)[[1]], names(in_list[[1]][["inside"]])), 2)
  combined_df <- data.frame(
    list = c(),
    Ex = c(),
    Shared = c(),
    pair = c()
  )

  # Loop through pairs
  i <- 1
  combo <- split(in_combos, col(in_combos))[[1]]
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
    if (is.null(focus_on)) {
      combined_df <- rbind(combined_df, current_df)
      i <- i + 1
    } else if (focus_on %in% c(bed_a_old_name, bed_b_new_name)) {
      combined_df <- rbind(combined_df, current_df)
      i <- i + 1
    }
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
  combined_df_long$group <- factor(
    combined_df_long$group,
    levels = names(colors_)
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
      values = colors_,
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

make_bars_1d <- function(
  cent_df,
  telo_df,
  s_arms_df,
  n_bams_names,
  colors_ = list_colors_fig_1
) {

  df_1d <- as.data.frame(
    t(rbind(cent_df, telo_df, s_arms_df))
  )

  # Correct Names
  lists_1d <- n_bams_names
  rownames(df_1d) <- lists_1d

  colnames(df_1d) <- list(
    "cent_df" = "Centromeres",
    "telo_df" = "Telomeres",
    "s_arms_df" = "Short Arms"
  )[colnames(df_1d)]

  df_1d <- as.data.frame(t(df_1d))

  df_1d$feature <- factor(
    c("Centromeres", "Telomeres", "Short Arms"),
    levels = c("Centromeres", "Telomeres", "Short Arms")
  )

  df_1d_long <- df_1d %>%
    tidyr::pivot_longer(
      cols = unname(lists_1d),
      names_to = "list",
      values_to = "value"
    )
  df_1d_long$list <- factor(
    df_1d_long$list,
    levels = rev(unname(names_dict))
  )
  out_plot_list <- list()
  for (gap in c("Centromeres", "Telomeres", "Short Arms")) {
    plot <- ggplot(
      df_1d_long[df_1d_long$feature == gap, ],
      aes(
        x = feature,
        y = value,
        fill = list,
        group = list
      )
    ) +
      geom_bar(
        stat = "identity",
        position = position_dodge()
      ) +
      scale_fill_manual(
        values = colors_,
        name = NULL,
        guide = guide_legend(
          override.aes = list(pattern = "none"),
          reverse = TRUE
        )
      ) +
      labs(
        x = NULL,
        y = NULL
      ) +
      theme_minimal() +
      theme(
        strip.text = element_blank(),
        strip.background = element_blank(),
        legend.position = "bottom"
      ) +
      coord_flip()
    out_plot_list[[length(out_plot_list) + 1]] <- plot
  }
  out_plot <- wrap_plots(out_plot_list, nrow = 3) +
    plot_layout(guides = "collect") &
    theme(legend.position = "bottom")
  return(out_plot)
}

get_p_c_plots <- function(
  in_matrix,
  title_,
  text_size = text_size,
  point_size = point_size,
  parameter_test = FALSE,
  shape_ = 13,
  colors_ = list_colors_fig_2
) {

  # Principal Coordinates
  loc <- cmdscale(in_matrix)
  loc_df <- data.frame(x = loc[, 1], y = -loc[, 2], label = rownames(loc))

  # Labels
  p1 <- ggplot(
    loc_df,
    aes(x = x, y = y, label = label)
  ) +
    geom_point(
      size = point_size,
      aes(color = label)
    ) +
    scale_color_manual(values = colors_[loc_df$label]) +
    geom_text_repel(
      #min.segment.length = 0,
      size = text_size,
      max.overlaps = 99999
    ) +
    theme_minimal() +
    theme(
      legend.position = "none"
    ) +
    labs(
      title = title_,
      subtitle = "Principal Coordinates",
      x = "Principal Coordinate 1",
      y = "Principal Coordinate 2"
    )
  return(p1)
}

get_para_mds_plots <- function(
  in_matrix,
  title_,
  text_size = 3,
  point_size = 2,
  color_by = "BAM Files",
  prefix_ = "STAR 36-BAMs "
) {

  # Principal Coordinates
  loc <- cmdscale(in_matrix)
  loc_df <- data.frame(x = loc[, 1], y = -loc[, 2], label = rownames(loc))

  para <- rownames(loc)[grep("000", rownames(loc), invert = FALSE)]
  not_para <- rownames(loc)[grep("000", rownames(loc), invert = TRUE)]

  labels <- gsub(prefix_, "", para)
  labels <- gsub("20000", "20k", labels)
  labels <- gsub("10000", "10k", labels)
  labels <- gsub("1000", "1k", labels)
  labels <- gsub("b", "", labels)
  labels <- gsub(" k", " ", labels)
  labels <- gsub("n", "", labels)

  labels_split <- strsplit(labels, " ")
  labels_split_matrix <- do.call(rbind, labels_split)
  labels_df <- as.data.frame(labels_split_matrix, stringsAsFactors = FALSE)
  rownames(labels_df) <- para
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
  not_para_df <- data.frame(matrix(ncol = 3, nrow = length(not_para)), row.names = not_para)
  colnames(not_para_df) <- colnames(labels_df)
  labels_df <- rbind(labels_df, not_para_df)
  loc_df <- cbind(loc_df, labels_df)
  #loc_df <- na.exclude(loc_df)

  colors_ <- ggsci::pal_lancet(
    palette = "lanonc",
    alpha = 1
  )(
    length(unique(loc_df[[color_by]]))
  )
  names(colors_) <- unique(loc_df[[color_by]])

  # Labels
  p1 <- ggplot(
    loc_df,
    aes(
      x = x,
      y = y,
      color = .data[[color_by]]
    )
  ) +
    geom_point(
      size = point_size
    ) +
    theme_minimal() +
    theme(
      legend.position = "right"
    ) +
    labs(
      title = title_,
      subtitle = "Principal Coordinates",
      x = "Principal Coordinate 1",
      y = "Principal Coordinate 2",
      color = color_by
    ) +
    scale_color_discrete(
      type = colors_
    )
  return(p1)
}

get_para_mds_plots_e <- function(
  in_df,
  title_,
  text_size = 3,
  point_size = 2,
  color_by = "BAM Files",
  prefix_ = "hg38 STAR 36bp "
) {

  # Principal Coordinates
  loc_df <- in_df

  labels <- gsub(prefix_, "", in_df$label)
  labels <- gsub("20000", "20k", labels)
  labels <- gsub("10000", "10k", labels)
  labels <- gsub("1000", "1k", labels)
  labels <- gsub("b", "", labels)
  labels <- gsub(" k", " ", labels)
  labels <- gsub("n", "", labels)

  labels_split <- strsplit(labels, " ")
  labels_split_matrix <- do.call(rbind, labels_split)
  labels_df <- as.data.frame(labels_split_matrix, stringsAsFactors = FALSE)
  rownames(labels_df) <- in_df$label
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

  loc_df <- cbind(loc_df, labels_df)
  loc_df <- na.exclude(loc_df)

  colors_ <- ggsci::pal_lancet(
    palette = "lanonc",
    alpha = 1
  )(
    length(unique(loc_df[[color_by]]))
  )
  names(colors_) <- unique(loc_df[[color_by]])

  # Labels
  p1 <- ggplot(
    loc_df,
    aes(
      x = x,
      y = y,
      color = .data[[color_by]]
    )
  ) +
    geom_point(
      size = point_size
    ) +
    theme_minimal() +
    theme(
      legend.position = "right"
    ) +
    labs(
      title = title_,
      subtitle = "Principal Coordinate Analysis",
      x = "Principal Coordinate 1",
      y = "Principal Coordinate 2",
      color = color_by
    ) +
    scale_color_discrete(
      type = colors_
    )
  return(p1)
}

get_para_mds_plots_labels <- function(
  in_matrix,
  title_,
  text_size = 3,
  point_size = 2,
  color_by = "BAM Files",
  prefix_ = "STAR 36-BAMs ",
  colors_ = list_colors_fig_2
) {

  # Principal Coordinates
  loc <- cmdscale(in_matrix)
  loc_df <- data.frame(x = loc[, 1], y = -loc[, 2], label = rownames(loc))

  labels <- gsub(prefix_, "", rownames(loc))
  labels <- gsub("20000", "20k", labels)
  labels <- gsub("10000", "10k", labels)
  labels <- gsub("1000", "1k", labels)
  labels <- gsub("b", "", labels)
  labels <- gsub(" k", " ", labels)
  labels <- gsub("n", "", labels)

  labels_split <- strsplit(labels, " ")
  labels_split_matrix <- do.call(rbind, labels_split)
  labels_df <- as.data.frame(labels_split_matrix, stringsAsFactors = FALSE)
  rownames(labels_df) <- rownames(loc)
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

  loc_df <- cbind(loc_df, labels_df)
  loc_df <- na.exclude(loc_df)

  colors_ <- ggsci::pal_lancet(
    palette = "lanonc",
    alpha = 1
  )(
    length(unique(loc_df[[color_by]]))
  )
  names(colors_) <- unique(loc_df[[color_by]])

  # Labels
  p1 <- ggplot(
    loc_df,
    aes(
      x = x,
      y = y,
      color = .data[[color_by]]
    )
  ) +
    geom_point(
      size = point_size
    ) +
    theme_minimal() +
    theme(
      legend.position = "right"
    ) +
    labs(
      title = title_,
      subtitle = "Principal Coordinates",
      x = "Principal Coordinate 1",
      y = "Principal Coordinate 2",
      color = color_by
    ) +
    scale_color_discrete(
      type = colors_
    )
  return(p1)
}

get_para_mds_plots_e_labels <- function(
  in_df,
  title_,
  text_size = 3,
  point_size = 2,
  color_by = "BAM Files",
  prefix_ = "hg38 STAR 36bp ",
  list_colors = list_colors_fig_2
) {

  # Principal Coordinates
  loc_df <- in_df

  labels <- gsub(prefix_, "", in_df$label)
  labels <- gsub("20000", "20k", labels)
  labels <- gsub("10000", "10k", labels)
  labels <- gsub("1000", "1k", labels)
  labels <- gsub("b", "", labels)
  labels <- gsub(" k", " ", labels)
  labels <- gsub("n", "", labels)

  labels_split <- strsplit(labels, " ")
  labels_split_matrix <- do.call(rbind, labels_split)
  labels_df <- as.data.frame(labels_split_matrix, stringsAsFactors = FALSE)
  rownames(labels_df) <- in_df$label
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

  loc_df <- cbind(loc_df, labels_df)
  loc_df <- na.exclude(loc_df)

  # Labels
  p1 <- ggplot(
    loc_df,
    aes(
      x = x,
      y = y,
      color = labels
    )
  ) +
    geom_point(
      size = point_size
    ) +
    theme_minimal() +
    theme(
      legend.position = "right"
    ) +
    labs(
      title = title_,
      subtitle = "Principal Coordinate Analysis",
      x = "Principal Coordinate 1",
      y = "Principal Coordinate 2",
    ) +
    scale_color_manual(values = list_colors[loc_df$label])
  return(p1)
}

get_ps_plots <- function(
  df,
  out_col,
  x_col,
  y_col,
  z_col,
  text_size = 4,
  ratio_ = 1
) {
  out_plots <- sapply(levels(df[[out_col]]), function(x) {
    sub_df <- df[df[[out_col]] == x, ]
    ggplot(
      sub_df,
      aes(
        x = .data[[x_col]],
        y = .data[[y_col]],
        fill = .data[[z_col]],
        label = .data[[z_col]]
      )
    ) +
      geom_tile(color = "white") +
      coord_fixed() +
      theme_minimal() +
      labs(
        subtitle = paste0(out_col, ": ", x),
        x = x_col,
        y = y_col
      ) +
      theme(
        legend.position = "none"
      ) +
      #geom_text(
      #  aes(
      #    label = formatC(
      #      sub_df[[z_col]],
      #      format = "G",
      #      digits = 2
      #    )
      #  ),
      #  color = "white",
      #  size = text_size
      #) +
      scale_fill_continuous(trans = "reverse") +
      coord_fixed(ratio = ratio_)

  }, simplify = FALSE)
  return(out_plots)
}

clean_row_labels <- function(
  plot_list,
  gsub_title = ""
) {
  out_list <- c(
    list(plot_list[[1]]),
    sapply(plot_list[-1], function(p) {
      o_p <- p + labs(
        subtitle = gsub(gsub_title, "", p$labels$subtitle),
        x = NULL,
        y = NULL
      ) +
        theme(
          axis.text.y = element_blank()
        )
      return(o_p)
    }, simplify = FALSE)
  )
  return(out_list)
}

get_dendro_plot <- function(hc, title_, subtitle_) {
  # Create dendrogram data using ggdendro
  hc_dendro <- as.dendrogram(hc)
  hc_data <- dendro_data(hc_dendro)

  # Create a ggplot dendrogram object with labels
  out_plot <- ggplot(hc_data$segments) +
    geom_segment(
      aes(
        x = x,
        y = y,
        xend = xend,
        yend = yend
      ),
    ) +
    labs(
      title = title_,
      subtitle = subtitle_,
      x = NULL, y = NULL
    ) +
    scale_x_continuous(
      breaks = hc_data$labels$x,
      labels = hc_data$labels$label
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_blank()
    ) +
    coord_flip()
  return(out_plot)
}

#load gaps
# Read in gaps as grange objects
centromeres <- bed_to_GRanges(path_to_centromeres, merge_gaps = TRUE)
telomeres <- bed_to_GRanges(path_to_telomeres, merge_gaps = TRUE)
short_arms <- bed_to_GRanges(path_to_shortarms, merge_gaps = TRUE)
gaps_gr <- reduce(c(centromeres, telomeres, short_arms), ignore.strand = TRUE)

mm10_centromeres <- bed_to_GRanges(mm10_path_to_centromeres, merge_gaps = TRUE)
mm10_telomeres <- bed_to_GRanges(mm10_path_to_telomeres, merge_gaps = TRUE)
mm10_short_arms <- bed_to_GRanges(mm10_path_to_shortarms, merge_gaps = TRUE)
mm10_gaps_gr <- reduce(c(mm10_centromeres, mm10_telomeres, mm10_short_arms), ignore.strand = TRUE)

#filter gaps
# Filter
centromeres <- keepSeqlevels(
  centromeres,
  include_chr[include_chr %in% seqnames(centromeres)],
  pruning.mode = "tidy"
)
telomeres <- keepSeqlevels(
  telomeres,
  include_chr[include_chr %in% seqnames(telomeres)],
  pruning.mode = "tidy"
)
short_arms <- keepSeqlevels(
  short_arms,
  include_chr[include_chr %in% seqnames(short_arms)],
  pruning.mode = "tidy"
)

mm10_centromeres <- keepSeqlevels(
  mm10_centromeres,
  mm10_include_chr[mm10_include_chr %in% seqnames(mm10_centromeres)],
  pruning.mode = "tidy"
)
mm10_telomeres <- keepSeqlevels(
  mm10_telomeres,
  mm10_include_chr[mm10_include_chr %in% seqnames(mm10_telomeres)],
  pruning.mode = "tidy"
)
mm10_short_arms <- keepSeqlevels(
  mm10_short_arms,
  mm10_include_chr[mm10_include_chr %in% seqnames(mm10_short_arms)],
  pruning.mode = "tidy"
)
# load bio
get_genome_knownGenes(
  exclude_chr = c(),
  genome = "hg38",
  gene_data = TxDb.Hsapiens.UCSC.hg38.knownGene
)
get_genome_knownGenes(
  exclude_chr = c(),
  genome = "mm10",
  gene_data = TxDb.Mmusculus.UCSC.mm10.knownGene
)

genes_gr <- bed_to_GRanges(genes_file, merge_gaps = TRUE)
exons_gr <- bed_to_GRanges(exons_file, merge_gaps = TRUE)
cds_gr <- bed_to_GRanges(cds_file, merge_gaps = TRUE)
bio_gr <- reduce(c(genes_gr, exons_gr, cds_gr), ignore.strand = TRUE)

mm10_genes_gr <- bed_to_GRanges(mm10_genes_file, merge_gaps = TRUE)
mm10_exons_gr <- bed_to_GRanges(mm10_exons_file, merge_gaps = TRUE)
mm10_cds_gr <- bed_to_GRanges(mm10_cds_file, merge_gaps = TRUE)
mm10_bio_gr <- reduce(c(mm10_genes_gr, mm10_exons_gr, mm10_cds_gr), ignore.strand = TRUE)
# filter bio, eval=FALSE
mm10_gaps_gr <- keepSeqlevels(
  mm10_gaps_gr,
  mm10_include_chr[mm10_include_chr %in% seqnames(mm10_gaps_gr)],
  pruning.mode = "tidy"
)
mm10_genes_gr <- keepSeqlevels(
  mm10_genes_gr,
  mm10_include_chr[mm10_include_chr %in% seqnames(mm10_genes_gr)],
  pruning.mode = "tidy"
)
exons_gr <- keepSeqlevels(
  mm10_exons_gr,
  mm10_include_chr[mm10_include_chr %in% seqnames(mm10_exons_gr)],
  pruning.mode = "tidy"
)
mm10_cds_gr <- keepSeqlevels(
  mm10_cds_gr,
  mm10_include_chr[mm10_include_chr %in% seqnames(mm10_cds_gr)],
  pruning.mode = "tidy"
)

# make vars
n_bams_all_fig_1 <- make_excludeList(
  path_to_sets = figure_1_dir,
  remove_gaps = FALSE,
  path_to_gaps = NULL
)
n_bams_all_fig_2 <- make_excludeList(
  path_to_sets = figure_2_dir,
  remove_gaps = FALSE,
  path_to_gaps = NULL
)
n_bams_all_fig_2_exp <- make_excludeList(
  path_to_sets = figure_2_exp_dir,
  remove_gaps = FALSE,
  path_to_gaps = NULL
)
n_bams_all_fig_mm10 <- make_excludeList(
  path_to_sets = figure_mm10_dir,
  remove_gaps = FALSE,
  path_to_gaps = NULL
)
n_bams_all_fig_aligners <- make_excludeList(
  path_to_sets = figure_aligners_dir,
  remove_gaps = FALSE,
  path_to_gaps = NULL
)
n_bams_all_fig_aligners_2 <- make_excludeList(
  path_to_sets = figure_aligners_2_dir,
  remove_gaps = FALSE,
  path_to_gaps = NULL
)
n_bams_all_fig_parameters <- make_excludeList(
  path_to_sets = figure_parameters_dir,
  remove_gaps = FALSE,
  path_to_gaps = NULL
)
n_bams_all_fig_merged <- make_excludeList(
  path_to_sets = figure_merged_dir,
  remove_gaps = FALSE,
  path_to_gaps = NULL
)

# filter chrs
# Exclude chrs
n_bams_all_fig_1 <- lapply(n_bams_all_fig_1, function(gr) {
  gr <- keepSeqlevels(
    gr,
    include_chr[include_chr %in% seqnames(gr)],
    pruning.mode = "tidy"
  )
  return(gr)
})
n_bams_all_fig_2 <- lapply(n_bams_all_fig_2, function(gr) {
  gr <- keepSeqlevels(
    gr,
    include_chr[include_chr %in% seqnames(gr)],
    pruning.mode = "tidy"
  )
  return(gr)
})
n_bams_all_fig_2_exp <- lapply(n_bams_all_fig_2_exp, function(gr) {
  gr <- keepSeqlevels(
    gr,
    include_chr[include_chr %in% seqnames(gr)],
    pruning.mode = "tidy"
  )
  return(gr)
})
n_bams_all_fig_mm10 <- lapply(n_bams_all_fig_mm10, function(gr) {
  gr <- keepSeqlevels(
    gr,
    mm10_include_chr[mm10_include_chr %in% seqnames(gr)],
    pruning.mode = "tidy"
  )
  return(gr)
})
n_bams_all_fig_aligners <- lapply(n_bams_all_fig_aligners, function(gr) {
  gr <- keepSeqlevels(
    gr,
    include_chr[include_chr %in% seqnames(gr)],
    pruning.mode = "tidy"
  )
  return(gr)
})
n_bams_all_fig_aligners_2 <- lapply(n_bams_all_fig_aligners_2, function(gr) {
  gr <- keepSeqlevels(
    gr,
    include_chr[include_chr %in% seqnames(gr)],
    pruning.mode = "tidy"
  )
  return(gr)
})
n_bams_all_fig_parameters <- lapply(n_bams_all_fig_parameters, function(gr) {
  gr <- keepSeqlevels(
    gr,
    include_chr[include_chr %in% seqnames(gr)],
    pruning.mode = "tidy"
  )
  return(gr)
})
n_bams_all_fig_merged <- lapply(n_bams_all_fig_merged, function(gr) {
  gr <- keepSeqlevels(
    gr,
    include_chr[include_chr %in% seqnames(gr)],
    pruning.mode = "tidy"
  )
  return(gr)
})

# make vars fig 1
gaps_overlaps_fig_1 <- gr_overlaps_apply(n_bams_all_fig_1, gaps_gr)
cent_overlaps_fig_1 <- gr_overlaps_apply(n_bams_all_fig_1, centromeres)
telo_overlaps_fig_1 <- gr_overlaps_apply(n_bams_all_fig_1, telomeres)
s_arms_overlaps_fig_1 <- gr_overlaps_apply(n_bams_all_fig_1, short_arms)
#genes_overlaps_fig_1 <- gr_overlaps_apply(n_bams_all_fig_1, genes_gr)
#exons_overlaps_fig_1 <- gr_overlaps_apply(n_bams_all_fig_1, exons_gr)
#cds_overlaps_fig_1 <- gr_overlaps_apply(n_bams_all_fig_1, cds_gr)
#bio_overlaps_fig_1 <- gr_overlaps_apply(n_bams_all_fig_1, bio_gr)
#r make vars fig 2}
gaps_overlaps_fig_2 <- gr_overlaps_apply(n_bams_all_fig_2, gaps_gr)
cent_overlaps_fig_2 <- gr_overlaps_apply(n_bams_all_fig_2, centromeres)
telo_overlaps_fig_2 <- gr_overlaps_apply(n_bams_all_fig_2, telomeres)
s_arms_overlaps_fig_2 <- gr_overlaps_apply(n_bams_all_fig_2, short_arms)
#genes_overlaps_fig_2 <- gr_overlaps_apply(n_bams_all_fig_2, genes_gr)
#exons_overlaps_fig_2 <- gr_overlaps_apply(n_bams_all_fig_2, exons_gr)
#cds_overlaps_fig_2 <- gr_overlaps_apply(n_bams_all_fig_2, cds_gr)
#bio_overlaps_fig_2 <- gr_overlaps_apply(n_bams_all_fig_2, bio_gr)
#r make vars fig mm10}
gaps_overlaps_fig_2_exp <- gr_overlaps_apply(n_bams_all_fig_2_exp, gaps_gr)
cent_overlaps_fig_2_exp <- gr_overlaps_apply(n_bams_all_fig_2_exp, centromeres)
telo_overlaps_fig_2_exp <- gr_overlaps_apply(n_bams_all_fig_2_exp, telomeres)
s_arms_overlaps_fig_2_exp <- gr_overlaps_apply(n_bams_all_fig_2_exp, short_arms)
#genes_overlaps_fig_2_exp <- gr_overlaps_apply(n_bams_all_fig_2_exp, genes_gr)
#exons_overlaps_fig_2_exp <- gr_overlaps_apply(n_bams_all_fig_2_exp, exons_gr)
#cds_overlaps_fig_2_exp <- gr_overlaps_apply(n_bams_all_fig_2_exp, cds_gr)
#bio_overlaps_fig_2_exp <- gr_overlaps_apply(n_bams_all_fig_2_exp, bio_gr)
#r make vars fig mm10}
gaps_overlaps_fig_mm10 <- gr_overlaps_apply(n_bams_all_fig_mm10, mm10_gaps_gr)
cent_overlaps_fig_mm10 <- gr_overlaps_apply(n_bams_all_fig_mm10, mm10_centromeres)
telo_overlaps_fig_mm10 <- gr_overlaps_apply(n_bams_all_fig_mm10, mm10_telomeres)
s_arms_overlaps_fig_mm10 <- gr_overlaps_apply(n_bams_all_fig_mm10, mm10_short_arms)
#genes_overlaps_fig_mm10 <- gr_overlaps_apply(n_bams_all_fig_mm10, mm10_genes_gr)
#exons_overlaps_fig_mm10 <- gr_overlaps_apply(n_bams_all_fig_mm10, mm10_exons_gr)
#cds_overlaps_fig_mm10 <- gr_overlaps_apply(n_bams_all_fig_mm10, mm10_cds_gr)
#bio_overlaps_fig_mm10 <- gr_overlaps_apply(n_bams_all_fig_mm10, mm10_bio_gr)
#r make vars fig aligners}
gaps_overlaps_fig_aligners <- gr_overlaps_apply(n_bams_all_fig_aligners, gaps_gr)
cent_overlaps_fig_aligners <- gr_overlaps_apply(n_bams_all_fig_aligners, centromeres)
telo_overlaps_fig_aligners <- gr_overlaps_apply(n_bams_all_fig_aligners, telomeres)
s_arms_overlaps_fig_aligners <- gr_overlaps_apply(n_bams_all_fig_aligners, short_arms)
#genes_overlaps_fig_aligners <- gr_overlaps_apply(n_bams_all_fig_aligners, genes_gr)
#exons_overlaps_fig_aligners <- gr_overlaps_apply(n_bams_all_fig_aligners, exons_gr)
#cds_overlaps_fig_aligners <- gr_overlaps_apply(n_bams_all_fig_aligners, cds_gr)
#bio_overlaps_fig_aligners <- gr_overlaps_apply(n_bams_all_fig_aligners, bio_gr)
#r make vars fig aligners 2}
gaps_overlaps_fig_aligners_2 <- gr_overlaps_apply(n_bams_all_fig_aligners_2, gaps_gr)
cent_overlaps_fig_aligners_2 <- gr_overlaps_apply(n_bams_all_fig_aligners_2, centromeres)
telo_overlaps_fig_aligners_2 <- gr_overlaps_apply(n_bams_all_fig_aligners_2, telomeres)
s_arms_overlaps_fig_aligners_2 <- gr_overlaps_apply(n_bams_all_fig_aligners_2, short_arms)
#genes_overlaps_fig_aligners_2 <- gr_overlaps_apply(n_bams_all_fig_aligners_2, genes_gr)
#exons_overlaps_fig_aligners_2 <- gr_overlaps_apply(n_bams_all_fig_aligners_2, exons_gr)
#cds_overlaps_fig_aligners_2 <- gr_overlaps_apply(n_bams_all_fig_aligners_2, cds_gr)
#bio_overlaps_fig_aligners_2 <- gr_overlaps_apply(n_bams_all_fig_aligners_2, bio_gr)
#r make vars fig parameters}
gaps_overlaps_fig_parameters <- gr_overlaps_apply(n_bams_all_fig_parameters, gaps_gr)
cent_overlaps_fig_parameters <- gr_overlaps_apply(n_bams_all_fig_parameters, centromeres)
telo_overlaps_fig_parameters <- gr_overlaps_apply(n_bams_all_fig_parameters, telomeres)
s_arms_overlaps_fig_parameters <- gr_overlaps_apply(n_bams_all_fig_parameters, short_arms)
#genes_overlaps_fig_parameters <- gr_overlaps_apply(n_bams_all_fig_parameters, genes_gr)
#exons_overlaps_fig_parameters <- gr_overlaps_apply(n_bams_all_fig_parameters, exons_gr)
#cds_overlaps_fig_parameters <- gr_overlaps_apply(n_bams_all_fig_parameters, cds_gr)
#bio_overlaps_fig_parameters <- gr_overlaps_apply(n_bams_all_fig_parameters, bio_gr)
#r make vars fig merged}
gaps_overlaps_fig_merged <- gr_overlaps_apply(n_bams_all_fig_merged, gaps_gr)
cent_overlaps_fig_merged <- gr_overlaps_apply(n_bams_all_fig_merged, centromeres)
telo_overlaps_fig_merged <- gr_overlaps_apply(n_bams_all_fig_merged, telomeres)
s_arms_overlaps_fig_merged <- gr_overlaps_apply(n_bams_all_fig_merged, short_arms)
#genes_overlaps_fig_merged <- gr_overlaps_apply(n_bams_all_fig_merged, genes_gr)
#exons_overlaps_fig_merged <- gr_overlaps_apply(n_bams_all_fig_merged, exons_gr)
#cds_overlaps_fig_merged <- gr_overlaps_apply(n_bams_all_fig_merged, cds_gr)
#bio_overlaps_fig_merged <- gr_overlaps_apply(n_bams_all_fig_merged, bio_gr)
#r make lists}
n_bams_list_fig_1 <- list(
  all = n_bams_all_fig_1,

  # Gaps
  inside_gaps = gaps_overlaps_fig_1$inside,
  trans_gaps = gaps_overlaps_fig_1$trans,
  outside_gaps = gaps_overlaps_fig_1$outside,
  overlap_gaps = gaps_overlaps_fig_1$any,

  # Centromeres
  inside_cent = cent_overlaps_fig_1$inside,
  trans_cent = cent_overlaps_fig_1$trans,
  outside_cent = cent_overlaps_fig_1$outside,
  overlap_cent = cent_overlaps_fig_1$any,

  # Telomeres
  inside_telo = telo_overlaps_fig_1$inside,
  trans_telo = telo_overlaps_fig_1$trans,
  outside_telo = telo_overlaps_fig_1$outside,
  overlap_telo = telo_overlaps_fig_1$any,

  # Short Arms
  inside_s_arms = s_arms_overlaps_fig_1$inside,
  trans_s_arms = s_arms_overlaps_fig_1$trans,
  outside_s_arms = s_arms_overlaps_fig_1$outside,
  overlap_s_arms = s_arms_overlaps_fig_1$any

  # Bio
  #inside_bio = bio_overlaps_fig_1$inside,
  #trans_bio = bio_overlaps_fig_1$trans,
  #outside_bio = bio_overlaps_fig_1$outside,
  #overlap_bio = bio_overlaps_fig_1$any,

  # Genes
  #inside_genes = genes_overlaps_fig_1$inside,
  #trans_genes = genes_overlaps_fig_1$trans,
  #outside_genes = genes_overlaps_fig_1$outside,
  #overlap_genes = genes_overlaps_fig_1$any,

  # Exons
  #inside_exons = exons_overlaps_fig_1$inside,
  #trans_exons = exons_overlaps_fig_1$trans,
  #outside_exons = exons_overlaps_fig_1$outside,
  #overlap_exons = exons_overlaps_fig_1$any,

  # CDS
  #inside_cds = cds_overlaps_fig_1$inside,
  #trans_cds = cds_overlaps_fig_1$trans,
  #outside_cds = cds_overlaps_fig_1$outside,
  #overlap_cds = cds_overlaps_fig_1$any
)

# Make names
names_fig_1 <- make_names(
  path_to_sets = figure_1_dir,
  remove_prefix = FALSE,
  remove_version = FALSE
)

n_bams_list_fig_2 <- list(
  all = n_bams_all_fig_2,

  # Gaps
  inside_gaps = gaps_overlaps_fig_2$inside,
  trans_gaps = gaps_overlaps_fig_2$trans,
  outside_gaps = gaps_overlaps_fig_2$outside,
  overlap_gaps = gaps_overlaps_fig_2$any,

  # Centromeres
  inside_cent = cent_overlaps_fig_2$inside,
  trans_cent = cent_overlaps_fig_2$trans,
  outside_cent = cent_overlaps_fig_2$outside,
  overlap_cent = cent_overlaps_fig_2$any,

  # Telomeres
  inside_telo = telo_overlaps_fig_2$inside,
  trans_telo = telo_overlaps_fig_2$trans,
  outside_telo = telo_overlaps_fig_2$outside,
  overlap_telo = telo_overlaps_fig_2$any,

  # Short Arms
  inside_s_arms = s_arms_overlaps_fig_2$inside,
  trans_s_arms = s_arms_overlaps_fig_2$trans,
  outside_s_arms = s_arms_overlaps_fig_2$outside,
  overlap_s_arms = s_arms_overlaps_fig_2$any

  # Bio
  #inside_bio = bio_overlaps_fig_2$inside,
  #trans_bio = bio_overlaps_fig_2$trans,
  #outside_bio = bio_overlaps_fig_2$outside,
  #overlap_bio = bio_overlaps_fig_2$any,

  # Genes
  #inside_genes = genes_overlaps_fig_2$inside,
  #trans_genes = genes_overlaps_fig_2$trans,
  #outside_genes = genes_overlaps_fig_2$outside,
  #overlap_genes = genes_overlaps_fig_2$any,

  # Exons
  #inside_exons = exons_overlaps_fig_2$inside,
  #trans_exons = exons_overlaps_fig_2$trans,
  #outside_exons = exons_overlaps_fig_2$outside,
  #overlap_exons = exons_overlaps_fig_2$any,

  # CDS
  #inside_cds = cds_overlaps_fig_2$inside,
  #trans_cds = cds_overlaps_fig_2$trans,
  #outside_cds = cds_overlaps_fig_2$outside,
  #overlap_cds = cds_overlaps_fig_2$any
)

# Make names
names_fig_2 <- make_names(
  path_to_sets = figure_2_dir,
  remove_prefix = FALSE,
  remove_version = FALSE
)
n_bams_list_fig_2_exp <- list(
  all = n_bams_all_fig_2_exp,

  # Gaps
  inside_gaps = gaps_overlaps_fig_2_exp$inside,
  trans_gaps = gaps_overlaps_fig_2_exp$trans,
  outside_gaps = gaps_overlaps_fig_2_exp$outside,
  overlap_gaps = gaps_overlaps_fig_2_exp$any,

  # Centromeres
  inside_cent = cent_overlaps_fig_2_exp$inside,
  trans_cent = cent_overlaps_fig_2_exp$trans,
  outside_cent = cent_overlaps_fig_2_exp$outside,
  overlap_cent = cent_overlaps_fig_2_exp$any,

  # Telomeres
  inside_telo = telo_overlaps_fig_2_exp$inside,
  trans_telo = telo_overlaps_fig_2_exp$trans,
  outside_telo = telo_overlaps_fig_2_exp$outside,
  overlap_telo = telo_overlaps_fig_2_exp$any,

  # Short Arms
  inside_s_arms = s_arms_overlaps_fig_2_exp$inside,
  trans_s_arms = s_arms_overlaps_fig_2_exp$trans,
  outside_s_arms = s_arms_overlaps_fig_2_exp$outside,
  overlap_s_arms = s_arms_overlaps_fig_2_exp$any

  # Bio
  #inside_bio = bio_overlaps_fig_2_exp$inside,
  #trans_bio = bio_overlaps_fig_2_exp$trans,
  #outside_bio = bio_overlaps_fig_2_exp$outside,
  #overlap_bio = bio_overlaps_fig_2_exp$any,

  # Genes
  #inside_genes = genes_overlaps_fig_2_exp$inside,
  #trans_genes = genes_overlaps_fig_2_exp$trans,
  #outside_genes = genes_overlaps_fig_2_exp$outside,
  #overlap_genes = genes_overlaps_fig_2_exp$any,

  # Exons
  #inside_exons = exons_overlaps_fig_2_exp$inside,
  #trans_exons = exons_overlaps_fig_2_exp$trans,
  #outside_exons = exons_overlaps_fig_2_exp$outside,
  #overlap_exons = exons_overlaps_fig_2_exp$any,

  # CDS
  #inside_cds = cds_overlaps_fig_2_exp$inside,
  #trans_cds = cds_overlaps_fig_2_exp$trans,
  #outside_cds = cds_overlaps_fig_2_exp$outside,
  #overlap_cds = cds_overlaps_fig_2_exp$any
)

names_fig_2_exp <- make_names(
  path_to_sets = figure_2_exp_dir,
  remove_prefix = FALSE,
  remove_version = FALSE
)

n_bams_list_fig_mm10 <- list(
  all = n_bams_all_fig_mm10,

  # Gaps
  inside_gaps = gaps_overlaps_fig_mm10$inside,
  trans_gaps = gaps_overlaps_fig_mm10$trans,
  outside_gaps = gaps_overlaps_fig_mm10$outside,
  overlap_gaps = gaps_overlaps_fig_mm10$any,

  # Centromeres
  inside_cent = cent_overlaps_fig_mm10$inside,
  trans_cent = cent_overlaps_fig_mm10$trans,
  outside_cent = cent_overlaps_fig_mm10$outside,
  overlap_cent = cent_overlaps_fig_mm10$any,

  # Telomeres
  inside_telo = telo_overlaps_fig_mm10$inside,
  trans_telo = telo_overlaps_fig_mm10$trans,
  outside_telo = telo_overlaps_fig_mm10$outside,
  overlap_telo = telo_overlaps_fig_mm10$any,

  # Short Arms
  inside_s_arms = s_arms_overlaps_fig_mm10$inside,
  trans_s_arms = s_arms_overlaps_fig_mm10$trans,
  outside_s_arms = s_arms_overlaps_fig_mm10$outside,
  overlap_s_arms = s_arms_overlaps_fig_mm10$any

  # Bio
  #inside_bio = bio_overlaps_fig_mm10$inside,
  #trans_bio = bio_overlaps_fig_mm10$trans,
  #outside_bio = bio_overlaps_fig_mm10$outside,
  #overlap_bio = bio_overlaps_fig_mm10$any,

  # Genes
  #inside_genes = genes_overlaps_fig_mm10$inside,
  #trans_genes = genes_overlaps_fig_mm10$trans,
  #outside_genes = genes_overlaps_fig_mm10$outside,
  #overlap_genes = genes_overlaps_fig_mm10$any,

  # Exons
  #inside_exons = exons_overlaps_fig_mm10$inside,
  #trans_exons = exons_overlaps_fig_mm10$trans,
  #outside_exons = exons_overlaps_fig_mm10$outside,
  #overlap_exons = exons_overlaps_fig_mm10$any,

  # CDS
  #inside_cds = cds_overlaps_fig_mm10$inside,
  #trans_cds = cds_overlaps_fig_mm10$trans,
  #outside_cds = cds_overlaps_fig_mm10$outside,
  #overlap_cds = cds_overlaps_fig_mm10$any
)

# Make names
names_fig_mm10 <- make_names(
  path_to_sets = figure_mm10_dir,
  remove_prefix = FALSE,
  remove_version = FALSE
)

n_bams_list_fig_aligners <- list(
  all = n_bams_all_fig_aligners,

  # Gaps
  inside_gaps = gaps_overlaps_fig_aligners$inside,
  trans_gaps = gaps_overlaps_fig_aligners$trans,
  outside_gaps = gaps_overlaps_fig_aligners$outside,
  overlap_gaps = gaps_overlaps_fig_aligners$any,

  # Centromeres
  inside_cent = cent_overlaps_fig_aligners$inside,
  trans_cent = cent_overlaps_fig_aligners$trans,
  outside_cent = cent_overlaps_fig_aligners$outside,
  overlap_cent = cent_overlaps_fig_aligners$any,

  # Telomeres
  inside_telo = telo_overlaps_fig_aligners$inside,
  trans_telo = telo_overlaps_fig_aligners$trans,
  outside_telo = telo_overlaps_fig_aligners$outside,
  overlap_telo = telo_overlaps_fig_aligners$any,

  # Short Arms
  inside_s_arms = s_arms_overlaps_fig_aligners$inside,
  trans_s_arms = s_arms_overlaps_fig_aligners$trans,
  outside_s_arms = s_arms_overlaps_fig_aligners$outside,
  overlap_s_arms = s_arms_overlaps_fig_aligners$any

  # Bio
  #inside_bio = bio_overlaps_fig_aligners$inside,
  #trans_bio = bio_overlaps_fig_aligners$trans,
  #outside_bio = bio_overlaps_fig_aligners$outside,
  #overlap_bio = bio_overlaps_fig_aligners$any,

  # Genes
  #inside_genes = genes_overlaps_fig_aligners$inside,
  #trans_genes = genes_overlaps_fig_aligners$trans,
  #outside_genes = genes_overlaps_fig_aligners$outside,
  #overlap_genes = genes_overlaps_fig_aligners$any,

  # Exons
  #inside_exons = exons_overlaps_fig_aligners$inside,
  #trans_exons = exons_overlaps_fig_aligners$trans,
  #outside_exons = exons_overlaps_fig_aligners$outside,
  #overlap_exons = exons_overlaps_fig_aligners$any,

  # CDS
  #inside_cds = cds_overlaps_fig_aligners$inside,
  #trans_cds = cds_overlaps_fig_aligners$trans,
  #outside_cds = cds_overlaps_fig_aligners$outside,
  #overlap_cds = cds_overlaps_fig_aligners$any
)

# Make names
names_fig_aligners <- make_names(
  path_to_sets = figure_aligners_dir,
  remove_prefix = FALSE,
  remove_version = FALSE
)

n_bams_list_fig_aligners_2 <- list(
  all = n_bams_all_fig_aligners_2,

  # Gaps
  inside_gaps = gaps_overlaps_fig_aligners_2$inside,
  trans_gaps = gaps_overlaps_fig_aligners_2$trans,
  outside_gaps = gaps_overlaps_fig_aligners_2$outside,
  overlap_gaps = gaps_overlaps_fig_aligners_2$any,

  # Centromeres
  inside_cent = cent_overlaps_fig_aligners_2$inside,
  trans_cent = cent_overlaps_fig_aligners_2$trans,
  outside_cent = cent_overlaps_fig_aligners_2$outside,
  overlap_cent = cent_overlaps_fig_aligners_2$any,

  # Telomeres
  inside_telo = telo_overlaps_fig_aligners_2$inside,
  trans_telo = telo_overlaps_fig_aligners_2$trans,
  outside_telo = telo_overlaps_fig_aligners_2$outside,
  overlap_telo = telo_overlaps_fig_aligners_2$any,

  # Short Arms
  inside_s_arms = s_arms_overlaps_fig_aligners_2$inside,
  trans_s_arms = s_arms_overlaps_fig_aligners_2$trans,
  outside_s_arms = s_arms_overlaps_fig_aligners_2$outside,
  overlap_s_arms = s_arms_overlaps_fig_aligners_2$any

  # Bio
  #inside_bio = bio_overlaps_fig_aligners_2$inside,
  #trans_bio = bio_overlaps_fig_aligners_2$trans,
  #outside_bio = bio_overlaps_fig_aligners_2$outside,
  #overlap_bio = bio_overlaps_fig_aligners_2$any,

  # Genes
  #inside_genes = genes_overlaps_fig_aligners_2$inside,
  #trans_genes = genes_overlaps_fig_aligners_2$trans,
  #outside_genes = genes_overlaps_fig_aligners_2$outside,
  #overlap_genes = genes_overlaps_fig_aligners_2$any,

  # Exons
  #inside_exons = exons_overlaps_fig_aligners_2$inside,
  #trans_exons = exons_overlaps_fig_aligners_2$trans,
  #outside_exons = exons_overlaps_fig_aligners_2$outside,
  #overlap_exons = exons_overlaps_fig_aligners_2$any,

  # CDS
  #inside_cds = cds_overlaps_fig_aligners_2$inside,
  #trans_cds = cds_overlaps_fig_aligners_2$trans,
  #outside_cds = cds_overlaps_fig_aligners_2$outside,
  #overlap_cds = cds_overlaps_fig_aligners_2$any
)

# Make names
names_fig_aligners_2 <- make_names(
  path_to_sets = figure_aligners_2_dir,
  remove_prefix = FALSE,
  remove_version = FALSE
)

n_bams_list_fig_parameters <- list(
  all = n_bams_all_fig_parameters,

  # Gaps
  inside_gaps = gaps_overlaps_fig_parameters$inside,
  trans_gaps = gaps_overlaps_fig_parameters$trans,
  outside_gaps = gaps_overlaps_fig_parameters$outside,
  overlap_gaps = gaps_overlaps_fig_parameters$any,

  # Centromeres
  inside_cent = cent_overlaps_fig_parameters$inside,
  trans_cent = cent_overlaps_fig_parameters$trans,
  outside_cent = cent_overlaps_fig_parameters$outside,
  overlap_cent = cent_overlaps_fig_parameters$any,

  # Telomeres
  inside_telo = telo_overlaps_fig_parameters$inside,
  trans_telo = telo_overlaps_fig_parameters$trans,
  outside_telo = telo_overlaps_fig_parameters$outside,
  overlap_telo = telo_overlaps_fig_parameters$any,

  # Short Arms
  inside_s_arms = s_arms_overlaps_fig_parameters$inside,
  trans_s_arms = s_arms_overlaps_fig_parameters$trans,
  outside_s_arms = s_arms_overlaps_fig_parameters$outside,
  overlap_s_arms = s_arms_overlaps_fig_parameters$any

  # Bio
  #inside_bio = bio_overlaps_fig_parameters$inside,
  #trans_bio = bio_overlaps_fig_parameters$trans,
  #outside_bio = bio_overlaps_fig_parameters$outside,
  #overlap_bio = bio_overlaps_fig_parameters$any,

  # Genes
  #inside_genes = genes_overlaps_fig_parameters$inside,
  #trans_genes = genes_overlaps_fig_parameters$trans,
  #outside_genes = genes_overlaps_fig_parameters$outside,
  #overlap_genes = genes_overlaps_fig_parameters$any,

  # Exons
  #inside_exons = exons_overlaps_fig_parameters$inside,
  #trans_exons = exons_overlaps_fig_parameters$trans,
  #outside_exons = exons_overlaps_fig_parameters$outside,
  #overlap_exons = exons_overlaps_fig_parameters$any,

  # CDS
  #inside_cds = cds_overlaps_fig_parameters$inside,
  #trans_cds = cds_overlaps_fig_parameters$trans,
  #outside_cds = cds_overlaps_fig_parameters$outside,
  #overlap_cds = cds_overlaps_fig_parameters$any
)

# Make names
names_fig_parameters <- make_names(
  path_to_sets = figure_parameters_dir,
  remove_prefix = FALSE,
  remove_version = FALSE
)

n_bams_list_fig_merged <- list(
  all = n_bams_all_fig_merged,

  # Gaps
  inside_gaps = gaps_overlaps_fig_merged$inside,
  trans_gaps = gaps_overlaps_fig_merged$trans,
  outside_gaps = gaps_overlaps_fig_merged$outside,
  overlap_gaps = gaps_overlaps_fig_merged$any,

  # Centromeres
  inside_cent = cent_overlaps_fig_merged$inside,
  trans_cent = cent_overlaps_fig_merged$trans,
  outside_cent = cent_overlaps_fig_merged$outside,
  overlap_cent = cent_overlaps_fig_merged$any,

  # Telomeres
  inside_telo = telo_overlaps_fig_merged$inside,
  trans_telo = telo_overlaps_fig_merged$trans,
  outside_telo = telo_overlaps_fig_merged$outside,
  overlap_telo = telo_overlaps_fig_merged$any,

  # Short Arms
  inside_s_arms = s_arms_overlaps_fig_merged$inside,
  trans_s_arms = s_arms_overlaps_fig_merged$trans,
  outside_s_arms = s_arms_overlaps_fig_merged$outside,
  overlap_s_arms = s_arms_overlaps_fig_merged$any

  # Bio
  #inside_bio = bio_overlaps_fig_merged$inside,
  #trans_bio = bio_overlaps_fig_merged$trans,
  #outside_bio = bio_overlaps_fig_merged$outside,
  #overlap_bio = bio_overlaps_fig_merged$any,

  # Genes
  #inside_genes = genes_overlaps_fig_merged$inside,
  #trans_genes = genes_overlaps_fig_merged$trans,
  #outside_genes = genes_overlaps_fig_merged$outside,
  #overlap_genes = genes_overlaps_fig_merged$any,

  # Exons
  #inside_exons = exons_overlaps_fig_merged$inside,
  #trans_exons = exons_overlaps_fig_merged$trans,
  #outside_exons = exons_overlaps_fig_merged$outside,
  #overlap_exons = exons_overlaps_fig_merged$any,

  # CDS
  #inside_cds = cds_overlaps_fig_merged$inside,
  #trans_cds = cds_overlaps_fig_merged$trans,
  #outside_cds = cds_overlaps_fig_merged$outside,
  #overlap_cds = cds_overlaps_fig_merged$any
)

# Make names
names_fig_merged <- make_names(
  path_to_sets = figure_merged_dir,
  remove_prefix = FALSE,
  remove_version = FALSE
)

# Get unique names
unique_names_fig_1 <- names_fig_1 |> unique()
unique_names_fig_2 <- names_fig_2 |> unique()
unique_names_fig_2_exp <- names_fig_2_exp |> unique()
unique_names_fig_mm10 <- names_fig_mm10 |> unique()
unique_names_fig_aligners <- names_fig_aligners |> unique()
unique_names_fig_aligners_2 <- names_fig_aligners_2 |> unique()
unique_names_fig_parameters <- names_fig_parameters |> unique()
unique_names_fig_merged <- names_fig_merged |> unique()

if (!dir.exists("results")) {
  # Create the folder only if it doesn't exist
  dir.create("results")
}
# calculations fig 1

# All
num_total_fig_1 <- get_num_regions(n_bams_list_fig_1$all)

# Gaps
num_overlap_gaps_fig_1 <- get_num_regions(n_bams_list_fig_1$overlap_gaps)

# Centromeres
num_overlap_cent_fig_1 <- get_num_regions(n_bams_list_fig_1$overlap_cent)

# Telomeres
num_overlap_telo_fig_1 <- get_num_regions(n_bams_list_fig_1$overlap_telo)

# Short Arms
num_overlap_s_arms_fig_1 <- get_num_regions(n_bams_list_fig_1$overlap_s_arms)

# Bio
#num_overlap_bio_fig_1 <- get_num_regions(n_bams_list_fig_1$overlap_bio)

# Genes
#num_overlap_genes_fig_1 <- get_num_regions(n_bams_list_fig_1$overlap_genes)

# Exons
#num_overlap_exons_fig_1 <- get_num_regions(n_bams_list_fig_1$overlap_exons)

# CDS
#num_overlap_cds_fig_1 <- get_num_regions(n_bams_list_fig_1$overlap_cds)

# Centromeres
cent_tw <- centromeres |> width() |> sum()
t_w_overlap_cent_fig_1 <- get_t_w(
  n_bams_list_fig_1$all,
  centromeres
)

# Telomeres
telo_tw <- telomeres |> width() |> sum()
t_w_overlap_telo_fig_1 <- get_t_w(
  n_bams_list_fig_1$all,
  telomeres
)

# Short Arms
s_arms_tw <- short_arms |> width() |> sum()
t_w_overlap_s_arms_fig_1 <- get_t_w(
  n_bams_list_fig_1$all,
  short_arms
)

# Bio
#bio_tw <- bio_gr |> width() |> sum()
#t_w_overlap_bio_fig_1 <- get_t_w(
#  n_bams_list_fig_1$all,
#  bio_gr
#)

# Genes
#genes_tw <- genes_gr |> width() |> sum()
#t_w_overlap_genes_fig_1 <- get_t_w(
#  n_bams_list_fig_1$all,
#  genes_gr
#)

# Exons
#exons_tw <- exons_gr |> width() |> sum()
#t_w_overlap_exons_fig_1 <- get_t_w(
#  n_bams_list_fig_1$all,
#  exons_gr
#)

# CDS
#cds_tw <- cds_gr |> width() |> sum()
#t_w_overlap_cds_fig_1 <- get_t_w(
#  n_bams_list_fig_1$all,
#  cds_gr
#)
#r calculations fig 2}
# All
num_total_fig_2 <- get_num_regions(n_bams_list_fig_2$all)

# Gaps
num_overlap_gaps_fig_2 <- get_num_regions(n_bams_list_fig_2$overlap_gaps)

# Centromeres
num_overlap_cent_fig_2 <- get_num_regions(n_bams_list_fig_2$overlap_cent)

# Telomeres
num_overlap_telo_fig_2 <- get_num_regions(n_bams_list_fig_2$overlap_telo)

# Short Arms
num_overlap_s_arms_fig_2 <- get_num_regions(n_bams_list_fig_2$overlap_s_arms)

# Bio
#num_overlap_bio_fig_2 <- get_num_regions(n_bams_list_fig_2$overlap_bio)

# Genes
#num_overlap_genes_fig_2 <- get_num_regions(n_bams_list_fig_2$overlap_genes)

# Exons
#num_overlap_exons_fig_2 <- get_num_regions(n_bams_list_fig_2$overlap_exons)

# CDS
#num_overlap_cds_fig_2 <- get_num_regions(n_bams_list_fig_2$overlap_cds)

# Centromeres
cent_tw <- centromeres |> width() |> sum()
t_w_overlap_cent_fig_2 <- get_t_w(
  n_bams_list_fig_2$all,
  centromeres
)

# Telomeres
telo_tw <- telomeres |> width() |> sum()
t_w_overlap_telo_fig_2 <- get_t_w(
  n_bams_list_fig_2$all,
  telomeres
)

# Short Arms
s_arms_tw <- short_arms |> width() |> sum()
t_w_overlap_s_arms_fig_2 <- get_t_w(
  n_bams_list_fig_2$all,
  short_arms
)

# Bio
#bio_tw <- bio_gr |> width() |> sum()
#t_w_overlap_bio_fig_2 <- get_t_w(
#  n_bams_list_fig_2$all,
#  bio_gr
#)

# Genes
#genes_tw <- genes_gr |> width() |> sum()
#t_w_overlap_genes_fig_2 <- get_t_w(
#  n_bams_list_fig_2$all,
#  genes_gr
#)

# Exons
#exons_tw <- exons_gr |> width() |> sum()
#t_w_overlap_exons_fig_2 <- get_t_w(
#  n_bams_list_fig_2$all,
#  exons_gr
#)

# CDS
#cds_tw <- cds_gr |> width() |> sum()
#t_w_overlap_cds_fig_2 <- get_t_w(
#  n_bams_list_fig_2$all,
#  cds_gr
#)
#r calculations fig mm10}

# All
num_total_fig_2_exp <- get_num_regions(n_bams_list_fig_2_exp$all)

# Gaps
num_overlap_gaps_fig_2_exp <- get_num_regions(n_bams_list_fig_2_exp$overlap_gaps)

# Centromeres
num_overlap_cent_fig_2_exp <- get_num_regions(n_bams_list_fig_2_exp$overlap_cent)

# Telomeres
num_overlap_telo_fig_2_exp <- get_num_regions(n_bams_list_fig_2_exp$overlap_telo)

# Short Arms
num_overlap_s_arms_fig_2_exp <- get_num_regions(n_bams_list_fig_2_exp$overlap_s_arms)

# Bio
#num_overlap_bio_fig_2_exp <- get_num_regions(n_bams_list_fig_2_exp$overlap_bio)

# Genes
#num_overlap_genes_fig_2_exp <- get_num_regions(n_bams_list_fig_2_exp$overlap_genes)

# Exons
#num_overlap_exons_fig_2_exp <- get_num_regions(n_bams_list_fig_2_exp$overlap_exons)

# CDS
#num_overlap_cds_fig_2_exp <- get_num_regions(n_bams_list_fig_2_exp$overlap_cds)

# Centromeres
cent_tw <- centromeres |> width() |> sum()
t_w_overlap_cent_fig_2_exp <- get_t_w(
  n_bams_list_fig_2_exp$all,
  centromeres
)

# Telomeres
telo_tw <- telomeres |> width() |> sum()
t_w_overlap_telo_fig_2_exp <- get_t_w(
  n_bams_list_fig_2_exp$all,
  telomeres
)

# Short Arms
s_arms_tw <- short_arms |> width() |> sum()
t_w_overlap_s_arms_fig_2_exp <- get_t_w(
  n_bams_list_fig_2_exp$all,
  short_arms
)

# Bio
#bio_tw <- bio_gr |> width() |> sum()
#t_w_overlap_bio_fig_2_exp <- get_t_w(
#  n_bams_list_fig_2_exp$all,
#  bio_gr
#)

# Genes
#genes_tw <- genes_gr |> width() |> sum()
#t_w_overlap_genes_fig_2_exp <- get_t_w(
#  n_bams_list_fig_2_exp$all,
#  genes_gr
#)

# Exons
#exons_tw <- exons_gr |> width() |> sum()
#t_w_overlap_exons_fig_2_exp <- get_t_w(
#  n_bams_list_fig_2_exp$all,
#  exons_gr
#)

# CDS
#cds_tw <- cds_gr |> width() |> sum()
#t_w_overlap_cds_fig_2_exp <- get_t_w(
#  n_bams_list_fig_2_exp$all,
#  cds_gr
#)
#r calculations fig mm10}

# All
num_total_fig_mm10 <- get_num_regions(n_bams_list_fig_mm10$all)

# Gaps
num_overlap_gaps_fig_mm10 <- get_num_regions(n_bams_list_fig_mm10$overlap_gaps)

# Centromeres
num_overlap_cent_fig_mm10 <- get_num_regions(n_bams_list_fig_mm10$overlap_cent)

# Telomeres
num_overlap_telo_fig_mm10 <- get_num_regions(n_bams_list_fig_mm10$overlap_telo)

# Short Arms
num_overlap_s_arms_fig_mm10 <- get_num_regions(n_bams_list_fig_mm10$overlap_s_arms)

# Bio
#num_overlap_bio_fig_mm10 <- get_num_regions(n_bams_list_fig_mm10$overlap_bio)

# Genes
#num_overlap_genes_fig_mm10 <- get_num_regions(n_bams_list_fig_mm10$overlap_genes)

# Exons
#num_overlap_exons_fig_mm10 <- get_num_regions(n_bams_list_fig_mm10$overlap_exons)

# CDS
#num_overlap_cds_fig_mm10 <- get_num_regions(n_bams_list_fig_mm10$overlap_cds)

# Centromeres
cent_tw <- mm10_centromeres |> width() |> sum()
t_w_overlap_cent_fig_mm10 <- get_t_w(
  n_bams_list_fig_mm10$all,
  mm10_centromeres
)

# Telomeres
telo_tw <- mm10_telomeres |> width() |> sum()
t_w_overlap_telo_fig_mm10 <- get_t_w(
  n_bams_list_fig_mm10$all,
  mm10_telomeres
)

# Short Arms
s_arms_tw <- mm10_short_arms |> width() |> sum()
t_w_overlap_s_arms_fig_mm10 <- get_t_w(
  n_bams_list_fig_mm10$all,
  mm10_short_arms
)

# Bio
#bio_tw <- mm10_bio_gr |> width() |> sum()
#t_w_overlap_bio_fig_mm10 <- get_t_w(
#  n_bams_list_fig_mm10$all,
#  mm10_bio_gr
#)

# Genes
#genes_tw <- mm10_genes_gr |> width() |> sum()
#t_w_overlap_genes_fig_mm10 <- get_t_w(
#  n_bams_list_fig_mm10$all,
#  mm10_genes_gr
#)

# Exons
#exons_tw <- mm10_exons_gr |> width() |> sum()
#t_w_overlap_exons_fig_mm10 <- get_t_w(
#  n_bams_list_fig_mm10$all,
#  mm10_exons_gr
#)

# CDS
#cds_tw <- mm10_cds_gr |> width() |> sum()
#t_w_overlap_cds_fig_mm10 <- get_t_w(
#  n_bams_list_fig_mm10$all,
#  mm10_cds_gr
#)
#r calculations fig aligners}

# All
num_total_fig_aligners <- get_num_regions(n_bams_list_fig_aligners$all)

# Gaps
num_overlap_gaps_fig_aligners <- get_num_regions(n_bams_list_fig_aligners$overlap_gaps)

# Centromeres
num_overlap_cent_fig_aligners <- get_num_regions(n_bams_list_fig_aligners$overlap_cent)

# Telomeres
num_overlap_telo_fig_aligners <- get_num_regions(n_bams_list_fig_aligners$overlap_telo)

# Short Arms
num_overlap_s_arms_fig_aligners <- get_num_regions(n_bams_list_fig_aligners$overlap_s_arms)

# Bio
#num_overlap_bio_fig_aligners <- get_num_regions(n_bams_list_fig_aligners$overlap_bio)

# Genes
#num_overlap_genes_fig_aligners <- get_num_regions(n_bams_list_fig_aligners$overlap_genes)

# Exons
#num_overlap_exons_fig_aligners <- get_num_regions(n_bams_list_fig_aligners$overlap_exons)

# CDS
#num_overlap_cds_fig_aligners <- get_num_regions(n_bams_list_fig_aligners$overlap_cds)

# Centromeres
cent_tw <- centromeres |> width() |> sum()
t_w_overlap_cent_fig_aligners <- get_t_w(
  n_bams_list_fig_aligners$all,
  centromeres
)

# Telomeres
telo_tw <- telomeres |> width() |> sum()
t_w_overlap_telo_fig_aligners <- get_t_w(
  n_bams_list_fig_aligners$all,
  telomeres
)

# Short Arms
s_arms_tw <- short_arms |> width() |> sum()
t_w_overlap_s_arms_fig_aligners <- get_t_w(
  n_bams_list_fig_aligners$all,
  short_arms
)

# Bio
#bio_tw <- bio_gr |> width() |> sum()
#t_w_overlap_bio_fig_aligners <- get_t_w(
#  n_bams_list_fig_aligners$all,
#  bio_gr
#)

# Genes
#genes_tw <- genes_gr |> width() |> sum()
#t_w_overlap_genes_fig_aligners <- get_t_w(
#  n_bams_list_fig_aligners$all,
#  genes_gr
#)

# Exons
#exons_tw <- exons_gr |> width() |> sum()
#t_w_overlap_exons_fig_aligners <- get_t_w(
#  n_bams_list_fig_aligners$all,
#  exons_gr
#)

# CDS
#cds_tw <- cds_gr |> width() |> sum()
#t_w_overlap_cds_fig_aligners <- get_t_w(
#  n_bams_list_fig_aligners$all,
#  cds_gr
#)
#r calculations fig aligners 2}

# All
num_total_fig_aligners_2 <- get_num_regions(n_bams_list_fig_aligners_2$all)

# Gaps
num_overlap_gaps_fig_aligners_2 <- get_num_regions(n_bams_list_fig_aligners_2$overlap_gaps)

# Centromeres
num_overlap_cent_fig_aligners_2 <- get_num_regions(n_bams_list_fig_aligners_2$overlap_cent)

# Telomeres
num_overlap_telo_fig_aligners_2 <- get_num_regions(n_bams_list_fig_aligners_2$overlap_telo)

# Short Arms
num_overlap_s_arms_fig_aligners_2 <- get_num_regions(n_bams_list_fig_aligners_2$overlap_s_arms)

# Bio
#num_overlap_bio_fig_aligners_2 <- get_num_regions(n_bams_list_fig_aligners_2$overlap_bio)

# Genes
#num_overlap_genes_fig_aligners_2 <- get_num_regions(n_bams_list_fig_aligners_2$overlap_genes)

# Exons
#num_overlap_exons_fig_aligners_2 <- get_num_regions(n_bams_list_fig_aligners_2$overlap_exons)

# CDS
#num_overlap_cds_fig_aligners_2 <- get_num_regions(n_bams_list_fig_aligners_2$overlap_cds)

# Centromeres
cent_tw <- centromeres |> width() |> sum()
t_w_overlap_cent_fig_aligners_2 <- get_t_w(
  n_bams_list_fig_aligners_2$all,
  centromeres
)

# Telomeres
telo_tw <- telomeres |> width() |> sum()
t_w_overlap_telo_fig_aligners_2 <- get_t_w(
  n_bams_list_fig_aligners_2$all,
  telomeres
)

# Short Arms
s_arms_tw <- short_arms |> width() |> sum()
t_w_overlap_s_arms_fig_aligners_2 <- get_t_w(
  n_bams_list_fig_aligners_2$all,
  short_arms
)

# Bio
#bio_tw <- bio_gr |> width() |> sum()
#t_w_overlap_bio_fig_aligners_2 <- get_t_w(
#  n_bams_list_fig_aligners_2$all,
#  bio_gr
#)

# Genes
#genes_tw <- genes_gr |> width() |> sum()
#t_w_overlap_genes_fig_aligners_2 <- get_t_w(
#  n_bams_list_fig_aligners_2$all,
#  genes_gr
#)

# Exons
#exons_tw <- exons_gr |> width() |> sum()
#t_w_overlap_exons_fig_aligners_2 <- get_t_w(
#  n_bams_list_fig_aligners_2$all,
#  exons_gr
#)

# CDS
#cds_tw <- cds_gr |> width() |> sum()
#t_w_overlap_cds_fig_aligners_2 <- get_t_w(
#  n_bams_list_fig_aligners_2$all,
#  cds_gr
#)
#r calculations fig merged}

# All
num_total_fig_merged <- get_num_regions(n_bams_list_fig_merged$all)

# Gaps
num_overlap_gaps_fig_merged <- get_num_regions(n_bams_list_fig_merged$overlap_gaps)

# Centromeres
num_overlap_cent_fig_merged <- get_num_regions(n_bams_list_fig_merged$overlap_cent)

# Telomeres
num_overlap_telo_fig_merged <- get_num_regions(n_bams_list_fig_merged$overlap_telo)

# Short Arms
num_overlap_s_arms_fig_merged <- get_num_regions(n_bams_list_fig_merged$overlap_s_arms)

# Bio
#num_overlap_bio_fig_merged <- get_num_regions(n_bams_list_fig_merged$overlap_bio)

# Genes
#num_overlap_genes_fig_merged <- get_num_regions(n_bams_list_fig_merged$overlap_genes)

# Exons
#num_overlap_exons_fig_merged <- get_num_regions(n_bams_list_fig_merged$overlap_exons)

# CDS
#num_overlap_cds_fig_merged <- get_num_regions(n_bams_list_fig_merged$overlap_cds)

# Centromeres
cent_tw <- centromeres |> width() |> sum()
t_w_overlap_cent_fig_merged <- get_t_w(
  n_bams_list_fig_merged$all,
  centromeres
)

# Telomeres
telo_tw <- telomeres |> width() |> sum()
t_w_overlap_telo_fig_merged <- get_t_w(
  n_bams_list_fig_merged$all,
  telomeres
)

# Short Arms
s_arms_tw <- short_arms |> width() |> sum()
t_w_overlap_s_arms_fig_merged <- get_t_w(
  n_bams_list_fig_merged$all,
  short_arms
)

# Bio
#bio_tw <- bio_gr |> width() |> sum()
#t_w_overlap_bio_fig_merged <- get_t_w(
#  n_bams_list_fig_merged$all,
#  bio_gr
#)

# Genes
#genes_tw <- genes_gr |> width() |> sum()
#t_w_overlap_genes_fig_merged <- get_t_w(
#  n_bams_list_fig_merged$all,
#  genes_gr
#)

# Exons
#exons_tw <- exons_gr |> width() |> sum()
#t_w_overlap_exons_fig_merged <- get_t_w(
#  n_bams_list_fig_merged$all,
#  exons_gr
#)

# CDS
#cds_tw <- cds_gr |> width() |> sum()
#t_w_overlap_cds_fig_merged <- get_t_w(
#  n_bams_list_fig_merged$all,
#  cds_gr
#)
#r calculations fig parameters}
# All
num_total_fig_parameters <- get_num_regions(n_bams_list_fig_parameters$all)

# Gaps
num_overlap_gaps_fig_parameters <- get_num_regions(n_bams_list_fig_parameters$overlap_gaps)

# Centromeres
num_overlap_cent_fig_parameters <- get_num_regions(n_bams_list_fig_parameters$overlap_cent)

# Telomeres
num_overlap_telo_fig_parameters <- get_num_regions(n_bams_list_fig_parameters$overlap_telo)

# Short Arms
num_overlap_s_arms_fig_parameters <- get_num_regions(n_bams_list_fig_parameters$overlap_s_arms)

# Bio
#num_overlap_bio_fig_parameters <- get_num_regions(n_bams_list_fig_parameters$overlap_bio)

# Genes
#num_overlap_genes_fig_parameters <- get_num_regions(n_bams_list_fig_parameters$overlap_genes)

# Exons
#num_overlap_exons_fig_parameters <- get_num_regions(n_bams_list_fig_parameters$overlap_exons)

# CDS
#num_overlap_cds_fig_parameters <- get_num_regions(n_bams_list_fig_parameters$overlap_cds)

# Centromeres
cent_tw <- centromeres |> width() |> sum()
t_w_overlap_cent_fig_parameters <- get_t_w(
  n_bams_list_fig_parameters$all,
  centromeres
)

# Telomeres
telo_tw <- telomeres |> width() |> sum()
t_w_overlap_telo_fig_parameters <- get_t_w(
  n_bams_list_fig_parameters$all,
  telomeres
)

# Short Arms
s_arms_tw <- short_arms |> width() |> sum()
t_w_overlap_s_arms_fig_parameters <- get_t_w(
  n_bams_list_fig_parameters$all,
  short_arms
)

# Bio
#bio_tw <- bio_gr |> width() |> sum()
#t_w_overlap_bio_fig_parameters <- get_t_w(
#  n_bams_list_fig_parameters$all,
#  bio_gr
#)

# Genes
#genes_tw <- genes_gr |> width() |> sum()
#t_w_overlap_genes_fig_parameters <- get_t_w(
#  n_bams_list_fig_parameters$all,
#  genes_gr
#)

# Exons
#exons_tw <- exons_gr |> width() |> sum()
#t_w_overlap_exons_fig_parameters <- get_t_w(
#  n_bams_list_fig_parameters$all,
#  exons_gr
#)

# CDS
#cds_tw <- cds_gr |> width() |> sum()
#t_w_overlap_cds_fig_parameters <- get_t_w(
#  n_bams_list_fig_parameters$all,
#  cds_gr
#)


# clustering fig 2 pt1
bed_list_fig_2 <- names(n_bams_list_fig_2$all)
names_fig_2 <- unname(names_dict[
  sub(paste0(".bed", "$"), "", bed_list_fig_2)
])
gr_list_fig_2 <- n_bams_list_fig_2$all[bed_list_fig_2]

# Jaccard width
# Make Jaccard width matrix
overlap_jaccard_count_fig_2 <- overlap_matrix(
  LIST = gr_list_fig_2,
  Names = names_fig_2,
  fun = jaccard_count,
  Vend = 1
)
# Jaccard width
# Make Jaccard width matrix
overlap_jaccard_width_fig_2 <- overlap_matrix(
  LIST = gr_list_fig_2,
  Names = names_fig_2,
  fun = jaccard_width,
  Vend = 1
)
# Forbes coefficient
# Forbes overlap matrix
overlap_forbes_width_fig_2 <- overlap_matrix(
  LIST = gr_list_fig_2,
  Names = names_fig_2,
  fun = forbes_coeff_width,
  Vend = 1,
  genome_size = genome_size,
  only_autosomal = FALSE
)
#r clustering fig 2 pt2}
# Perform hierarchical clustering using Ward's method
j_c_distance_matrix_fig_2 <- 1 - overlap_jaccard_count_fig_2 # Convert Jaccard coefficients to distances
diag(j_c_distance_matrix_fig_2) <- 0 # Set diagonal to 0 since distance from an element to itself is 0
j_c_distance_matrix_fig_2 <- as.dist(j_c_distance_matrix_fig_2)
overlap_jaccard_count_hc <- hclust(
  j_c_distance_matrix_fig_2,
  method = "ward.D2"
)

j_w_distance_matrix_fig_2 <- 1 - overlap_jaccard_width_fig_2 # Convert Jaccard coefficients to distances
diag(j_w_distance_matrix_fig_2) <- 0 # Set diagonal to 0 since distance from an element to itself is 0
j_w_distance_matrix_fig_2 <- as.dist(j_w_distance_matrix_fig_2)
overlap_jaccard_width_hc <- hclust(
  j_w_distance_matrix_fig_2,
  method = "ward.D2"
)


epsilon <- 1e-10 # Define a small constant epsilon to avoid division by zero
f_w_distance_matrix_fig_2 <- 1 / (overlap_forbes_width_fig_2 + epsilon) # Convert Forbes coefficients to distances
diag(f_w_distance_matrix_fig_2) <- 0 # Set diagonal to 0 since distance from an element to itself is 0
f_w_distance_matrix_fig_2 <- as.dist(f_w_distance_matrix_fig_2)
overlap_forbes_width_hc <- hclust(
  f_w_distance_matrix_fig_2,
  method = "ward.D2"
)

# clustering fig 2 pt1
bed_list_fig_2_exp <- names(n_bams_list_fig_2_exp$all)
names_fig_2_exp <- unname(names_dict[
  sub(paste0(".bed", "$"), "", bed_list_fig_2_exp)
])
gr_list_fig_2_exp <- n_bams_list_fig_2_exp$all[bed_list_fig_2_exp]

# Jaccard width
# Make Jaccard width matrix
overlap_jaccard_count_fig_2_exp <- overlap_matrix(
  LIST = gr_list_fig_2_exp,
  Names = names_fig_2_exp,
  fun = jaccard_count,
  Vend = 1
)
# Jaccard width
# Make Jaccard width matrix
overlap_jaccard_width_fig_2_exp <- overlap_matrix(
  LIST = gr_list_fig_2_exp,
  Names = names_fig_2_exp,
  fun = jaccard_width,
  Vend = 1
)
# Forbes coefficient
# Forbes overlap matrix
overlap_forbes_width_fig_2_exp <- overlap_matrix(
  LIST = gr_list_fig_2_exp,
  Names = names_fig_2_exp,
  fun = forbes_coeff_width,
  Vend = 1,
  genome_size = genome_size,
  only_autosomal = FALSE
)
#r clustering fig 2 pt2}
# Perform hierarchical clustering using Ward's method
j_c_distance_matrix_fig_2_exp <- 1 - overlap_jaccard_count_fig_2_exp # Convert Jaccard coefficients to distances
diag(j_c_distance_matrix_fig_2_exp) <- 0 # Set diagonal to 0 since distance from an element to itself is 0
j_c_distance_matrix_fig_2_exp <- as.dist(j_c_distance_matrix_fig_2_exp)
overlap_jaccard_count_fig_2_exp_hc <- hclust(
  j_c_distance_matrix_fig_2_exp,
  method = "ward.D2"
)

j_w_distance_matrix_fig_2_exp <- 1 - overlap_jaccard_width_fig_2_exp # Convert Jaccard coefficients to distances
diag(j_w_distance_matrix_fig_2_exp) <- 0 # Set diagonal to 0 since distance from an element to itself is 0
j_w_distance_matrix_fig_2_exp <- as.dist(j_w_distance_matrix_fig_2_exp)
overlap_jaccard_width_fig_2_exp_hc <- hclust(
  j_w_distance_matrix_fig_2_exp,
  method = "ward.D2"
)


epsilon <- 1e-10 # Define a small constant epsilon to avoid division by zero
f_w_distance_matrix_fig_2_exp <- 1 / (overlap_forbes_width_fig_2_exp + epsilon) # Convert Forbes coefficients to distances
diag(f_w_distance_matrix_fig_2_exp) <- 0 # Set diagonal to 0 since distance from an element to itself is 0
f_w_distance_matrix_fig_2_exp <- as.dist(f_w_distance_matrix_fig_2_exp)
overlap_forbes_width_fig_2_exp_hc <- hclust(
  f_w_distance_matrix_fig_2_exp,
  method = "ward.D2"
)

#r clustering fig aln pt1}

bed_list_fig_aligners <- names(n_bams_list_fig_aligners$all)
names_fig_aligners <- unname(names_dict[
  sub(paste0(".bed", "$"), "", bed_list_fig_aligners)
])
gr_list_fig_aligners <- n_bams_list_fig_aligners$all[bed_list_fig_aligners]

# Jaccard width
# Make Jaccard width matrix
overlap_jaccard_count_aln <- overlap_matrix(
  LIST = gr_list_fig_aligners,
  Names = names_fig_aligners,
  fun = jaccard_count,
  Vend = 1
)
# Jaccard width
# Make Jaccard width matrix
overlap_jaccard_width_aln <- overlap_matrix(
  LIST = gr_list_fig_aligners,
  Names = names_fig_aligners,
  fun = jaccard_width,
  Vend = 1
)
# Forbes coefficient
# Forbes overlap matrix
overlap_forbes_width_aln <- overlap_matrix(
  LIST = gr_list_fig_aligners,
  Names = names_fig_aligners,
  fun = forbes_coeff_width,
  Vend = 1,
  genome_size = genome_size,
  only_autosomal = FALSE
)
#r clustering fig aln pt2}
# Perform hierarchical clustering using Ward's method
j_c_distance_matrix_aln <- 1 - overlap_jaccard_count_aln # Convert Jaccard coefficients to distances
diag(j_c_distance_matrix_aln) <- 0 # Set diagonal to 0 since distance from an element to itself is 0
j_c_distance_matrix_aln <- as.dist(j_c_distance_matrix_aln)
overlap_jaccard_count_hc_aln <- hclust(
  j_c_distance_matrix_aln,
  method = "ward.D2"
)

j_w_distance_matrix_aln <- 1 - overlap_jaccard_width_aln # Convert Jaccard coefficients to distances
diag(j_w_distance_matrix_aln) <- 0 # Set diagonal to 0 since distance from an element to itself is 0
j_w_distance_matrix_aln <- as.dist(j_w_distance_matrix_aln)
overlap_jaccard_width_hc_aln <- hclust(
  j_w_distance_matrix_aln,
  method = "ward.D2"
)


epsilon <- 1e-10 # Define a small constant epsilon to avoid division by zero
f_w_distance_matrix_aln <- 1 / (overlap_forbes_width_aln + epsilon) # Convert Forbes coefficients to distances
diag(f_w_distance_matrix_aln) <- 0 # Set diagonal to 0 since distance from an element to itself is 0
f_w_distance_matrix_aln <- as.dist(f_w_distance_matrix_aln)
overlap_forbes_width_hc_aln <- hclust(
  f_w_distance_matrix_aln,
  method = "ward.D2"
)
#r clustering fig aln 2 pt1}

bed_list_fig_aligners_2 <- names(n_bams_list_fig_aligners_2$all)
names_fig_aligners_2 <- unname(names_dict[
  sub(paste0(".bed", "$"), "", bed_list_fig_aligners_2)
])
gr_list_fig_aligners_2 <- n_bams_list_fig_aligners_2$all[bed_list_fig_aligners_2]

# Jaccard width
# Make Jaccard width matrix
overlap_jaccard_count_aln_2 <- overlap_matrix(
  LIST = gr_list_fig_aligners_2,
  Names = names_fig_aligners_2,
  fun = jaccard_count,
  Vend = 1
)
# Jaccard width
# Make Jaccard width matrix
overlap_jaccard_width_aln_2 <- overlap_matrix(
  LIST = gr_list_fig_aligners_2,
  Names = names_fig_aligners_2,
  fun = jaccard_width,
  Vend = 1
)
# Forbes coefficient
# Forbes overlap matrix
overlap_forbes_width_aln_2 <- overlap_matrix(
  LIST = gr_list_fig_aligners_2,
  Names = names_fig_aligners_2,
  fun = forbes_coeff_width,
  Vend = 1,
  genome_size = genome_size,
  only_autosomal = FALSE
)
#r clustering fig aln 2 pt2}
# Perform hierarchical clustering using Ward's method
j_c_distance_matrix_aln_2 <- 1 - overlap_jaccard_count_aln_2 # Convert Jaccard coefficients to distances
diag(j_c_distance_matrix_aln_2) <- 0 # Set diagonal to 0 since distance from an element to itself is 0
j_c_distance_matrix_aln_2 <- as.dist(j_c_distance_matrix_aln_2)
overlap_jaccard_count_hc_aln_2 <- hclust(
  j_c_distance_matrix_aln_2,
  method = "ward.D2"
)

j_w_distance_matrix_aln_2 <- 1 - overlap_jaccard_width_aln_2 # Convert Jaccard coefficients to distances
diag(j_w_distance_matrix_aln_2) <- 0 # Set diagonal to 0 since distance from an element to itself is 0
j_w_distance_matrix_aln_2 <- as.dist(j_w_distance_matrix_aln_2)
overlap_jaccard_width_hc_aln_2 <- hclust(
  j_w_distance_matrix_aln_2,
  method = "ward.D2"
)

epsilon <- 1e-10 # Define a small constant epsilon to avoid division by zero
f_w_distance_matrix_aln_2 <- 1 / (overlap_forbes_width_aln_2 + epsilon) # Convert Forbes coefficients to distances
diag(f_w_distance_matrix_aln_2) <- 0 # Set diagonal to 0 since distance from an element to itself is 0
f_w_distance_matrix_aln_2 <- as.dist(f_w_distance_matrix_aln_2)
overlap_forbes_width_hc_aln_2 <- hclust(
  f_w_distance_matrix_aln_2,
  method = "ward.D2"
)
#r clustering fig parameters pt1}

bed_list_fig_parameters <- names(n_bams_list_fig_parameters$all)
names_fig_parameters <- unname(names_dict[
  sub(paste0(".bed", "$"), "", bed_list_fig_parameters)
])
gr_list_fig_parameters <- n_bams_list_fig_parameters$all[bed_list_fig_parameters]

# Jaccard width
# Make Jaccard width matrix
overlap_jaccard_count_para <- overlap_matrix(
  LIST = gr_list_fig_parameters,
  Names = names_fig_parameters,
  fun = jaccard_count,
  Vend = 1
)
# Jaccard width
# Make Jaccard width matrix
overlap_jaccard_width_para <- overlap_matrix(
  LIST = gr_list_fig_parameters,
  Names = names_fig_parameters,
  fun = jaccard_width,
  Vend = 1
)
# Forbes coefficient
# Forbes overlap matrix
overlap_forbes_width_para <- overlap_matrix(
  LIST = gr_list_fig_parameters,
  Names = names_fig_parameters,
  fun = forbes_coeff_width,
  Vend = 1,
  genome_size = genome_size,
  only_autosomal = FALSE
)
overlap_forbes_width_para_gs <- overlap_matrix(
  LIST = c(gr_list_fig_parameters, gr_list_fig_2),
  Names = c(names_fig_parameters, names_fig_2),
  fun = forbes_coeff_width,
  Vend = 1,
  genome_size = genome_size,
  only_autosomal = FALSE
)
#r clustering fig parameters pt2}
# Perform hierarchical clustering using Ward's method
j_c_distance_matrix_para <- 1 - overlap_jaccard_count_para # Convert Jaccard coefficients to distances
diag(j_c_distance_matrix_para) <- 0 # Set diagonal to 0 since distance from an element to itself is 0
j_c_distance_matrix_para <- as.dist(j_c_distance_matrix_para)
overlap_jaccard_count_para_hc <- hclust(
  j_c_distance_matrix_para,
  method = "ward.D2"
)

j_w_distance_matrix_para <- 1 - overlap_jaccard_width_para # Convert Jaccard coefficients to distances
diag(j_w_distance_matrix_para) <- 0 # Set diagonal to 0 since distance from an element to itself is 0
j_w_distance_matrix_para <- as.dist(j_w_distance_matrix_para)
overlap_jaccard_width_para_hc <- hclust(
  j_w_distance_matrix_para,
  method = "ward.D2"
)


epsilon <- 1e-10 # Define a small constant epsilon to avoid division by zero
f_w_distance_matrix_para <- 1 / (overlap_forbes_width_para + epsilon) # Convert Forbes coefficients to distances
diag(f_w_distance_matrix_para) <- 0 # Set diagonal to 0 since distance from an element to itself is 0
f_w_distance_matrix_para <- as.dist(f_w_distance_matrix_para)
overlap_forbes_width_para_hc <- hclust(
  f_w_distance_matrix_para,
  method = "ward.D2"
)

epsilon <- 1e-10 # Define a small constant epsilon to avoid division by zero
f_w_distance_matrix_para_gs <- 1 / (overlap_forbes_width_para_gs + epsilon) # Convert Forbes coefficients to distances
diag(f_w_distance_matrix_para_gs) <- 0 # Set diagonal to 0 since distance from an element to itself is 0
f_w_distance_matrix_para_gs <- as.dist(f_w_distance_matrix_para_gs)
overlap_forbes_width_para_gs_hc <- hclust(
  f_w_distance_matrix_para_gs,
  method = "ward.D2"
)

# figure 1 color setup
list_names_fig_1 <- unname(names_dict[
  sub(paste0(".bed", "$"), "", names(n_bams_all_fig_1))
])
list_colors_fig_1 <- ggsci::pal_lancet(
  palette = "lanonc",
  alpha = 1
)(
  length(list_names_fig_1)
)
# Define colors
list_colors_fig_1 <- c(list_colors_fig_1, "darkgrey", "black")
color_order_fig_1 <- c(list_names_fig_1, "Gaps", "Shared")
names(list_colors_fig_1) <- color_order_fig_1
list_colors_fig_1 <- list_colors_fig_1[c(names_dict, "Gaps", "Shared")]
list_colors_fig_1 <- na.omit(list_colors_fig_1)
#r figure 2 color setup}
list_names_fig_2 <- unname(names_dict[
  sub(paste0(".bed", "$"), "", names(n_bams_all_fig_2))
])
list_colors_fig_2 <- ggsci::pal_lancet(
  palette = "lanonc",
  alpha = 1
)(
  length(unique(
    c(list_names_fig_1, list_names_fig_2)
  ))
)
new_colors_fig_2 <- list_colors_fig_2[
  ! list_colors_fig_2 %in% list_colors_fig_1
]
new_names_fig_2 <- list_names_fig_2[
  ! list_names_fig_2 %in% names(list_colors_fig_1)
]
names(new_colors_fig_2) <- new_names_fig_2

# Define colors
list_colors_fig_2 <- c(list_colors_fig_1, new_colors_fig_2)
list_colors_fig_2 <- list_colors_fig_2[c(names_dict, "Gaps", "Shared")]
list_colors_fig_2 <- na.omit(list_colors_fig_2)

#r figure 2 exp color setup}
list_names_fig_2_exp <- unname(names_dict[
  sub(paste0(".bed", "$"), "", names(n_bams_all_fig_2_exp))
])
list_colors_fig_2_exp <- ggsci::pal_lancet(
  palette = "lanonc",
  alpha = 1
)(
  length(unique(
    c(list_names_fig_1, list_names_fig_2_exp)
  ))
)
new_colors_fig_2_exp <- list_colors_fig_2_exp[
  ! list_colors_fig_2_exp %in% list_colors_fig_2
]
new_names_fig_2_exp <- list_names_fig_2_exp[
  ! list_names_fig_2_exp %in% names(list_colors_fig_2)
]
names(new_colors_fig_2_exp) <- new_names_fig_2_exp

# Define colors
list_colors_fig_2_exp <- c(list_colors_fig_2, new_colors_fig_2_exp)
list_colors_fig_2_exp <- list_colors_fig_2_exp[c(names_dict, "Gaps", "Shared")]
list_colors_fig_2_exp <- na.omit(list_colors_fig_2_exp)

#r figure mm10 color setup}
list_names_fig_mm10 <- unname(names_dict[
  sub(paste0(".bed", "$"), "", names(n_bams_all_fig_mm10))
])
list_colors_fig_mm10 <- ggsci::pal_lancet(
  palette = "lanonc",
  alpha = 1
)(
  length(list_names_fig_mm10)
)
# Define colors
list_colors_fig_mm10 <- c(list_colors_fig_mm10, "darkgrey", "black")
color_order_fig_mm10 <- c(list_names_fig_mm10, "Gaps", "Shared")
names(list_colors_fig_mm10) <- color_order_fig_mm10
list_colors_fig_mm10 <- list_colors_fig_mm10[c(names_dict, "Gaps", "Shared")]
list_colors_fig_mm10 <- na.omit(list_colors_fig_mm10)
#r figure aligners color setup}
list_names_fig_aligners <- unname(names_dict[
  sub(paste0(".bed", "$"), "", names(n_bams_all_fig_aligners))
])
list_colors_fig_aligners <- ggsci::pal_lancet(
  palette = "lanonc",
  alpha = 1
)(
  length(list_names_fig_aligners)
)
# Define colors
list_colors_fig_aligners <- c(list_colors_fig_aligners, "darkgrey", "black")
color_order_fig_aligners <- c(list_names_fig_aligners, "Gaps", "Shared")
names(list_colors_fig_aligners) <- color_order_fig_aligners
list_colors_fig_aligners <- list_colors_fig_aligners[c(names_dict, "Gaps", "Shared")]
list_colors_fig_aligners <- na.omit(list_colors_fig_aligners)
#r figure aligners 2 color setup}
list_names_fig_aligners_2 <- unname(names_dict[
  sub(paste0(".bed", "$"), "", names(n_bams_all_fig_aligners_2))
])
list_colors_fig_aligners_2 <- ggsci::pal_lancet(
  palette = "lanonc",
  alpha = 1
)(
  length(list_names_fig_aligners_2)
)
# Define colors
list_colors_fig_aligners_2 <- c(list_colors_fig_aligners_2, "darkgrey", "black")
color_order_fig_aligners_2 <- c(list_names_fig_aligners_2, "Gaps", "Shared")
names(list_colors_fig_aligners_2) <- color_order_fig_aligners_2
list_colors_fig_aligners_2 <- list_colors_fig_aligners_2[c(names_dict, "Gaps", "Shared")]
list_colors_fig_aligners_2 <- na.omit(list_colors_fig_aligners_2)
#r figure merged color setup}
list_names_fig_merged <- unname(names_dict[
  sub(paste0(".bed", "$"), "", names(n_bams_all_fig_merged))
])
list_colors_fig_merged <- ggsci::pal_lancet(
  palette = "lanonc",
  alpha = 1
)(
  length(list_names_fig_merged)
)
# Define colors
list_colors_fig_merged <- c(list_colors_fig_merged, "darkgrey", "black")
color_order_fig_merged <- c(list_names_fig_merged, "Gaps", "Shared")
names(list_colors_fig_merged) <- color_order_fig_merged
list_colors_fig_merged <- list_colors_fig_merged[c(names_dict, "Gaps", "Shared")]
list_colors_fig_merged <- na.omit(list_colors_fig_merged)

