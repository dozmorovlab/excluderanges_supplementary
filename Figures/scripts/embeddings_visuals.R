library(uwot)
library(ggplot2)
library(ggrepel)
library(patchwork)

# Load embeddings
# Define the directory path
directory <- file.path("data", "embeddings")

# List all files in the directory (assuming they are .csv files)
file_list <- list.files(directory, pattern = "\\.csv$", full.names = TRUE)

# Create a named list of data frames
file_path <- "data/embeddings/figure_S4D.csv"

if (!exists("point_size")) {
  point_size <- 2
}
if (!exists("text_size")) {
  text_size <- 3
}

get_plot <- function(
  embeddings,
  out_mode = "MDS",
  preprocessed = FALSE,
  point_size_ = point_size,
  text_size_ = text_size
) {

  if (!preprocessed) {
    # Set the first column as row names and remove it from the data frame
    rownames(embeddings) <- gsub(
      ".bed|hg38 | STAR 1k| \\+ Centromere",
      "",
      embeddings[, 1]
    )
    embeddings <- embeddings[, -1]

    rownames(embeddings) <- gsub(
      " 36bp List",
      " 36bp",
      rownames(embeddings)
    )
    rownames(embeddings) <- gsub(
      " 101bp List",
      " 101bp",
      rownames(embeddings)
    )
    rownames(embeddings) <- gsub(
      " List",
      " Blacklist",
      rownames(embeddings)
    )
    rownames(embeddings) <- gsub(
      "HS \\+ LM \\+ CM",
      "HS \\+ LM",
      rownames(embeddings)
    )
  }

  if (out_mode == "PCA") {
    # Get PCA
    pca_data <- prcomp(embeddings)

    # Plot PCA
    pca_plot <- ggplot(
      data = data.frame(
        pca_data$x,
        embeddings,
        samples = rownames(embeddings),
        stringsAsFactors = FALSE
      ),
      aes(
        x = as.numeric(PC1),
        y = as.numeric(PC2),
        label = samples,
        color = samples
      )
    ) +
      geom_point(size = point_size_, alpha = 1) +
      geom_text_repel(
        color = "black",
        size = text_size_,
        max.overlaps = getOption(
          "ggrepel.max.overlaps",
          default = 100
        )
      ) +
      geom_hline(yintercept = 0, color = "gray65") +
      geom_vline(xintercept = 0, color = "gray65") +
      labs(
        #title = "Principal Component Analysis",
        title = NULL,
        color = NULL
      ) +
      scale_x_continuous(
        name = paste0(
          "PC1, ",
          round(summary(pca_data)$importance[2, 1] * 100, digits = 2),
          "% variability"
        )
      ) +
      scale_y_continuous(
        name = paste0(
          "PC2, ",
          round(summary(pca_data)$importance[2, 2] * 100, digits = 2),
          "% variability"
        )
      ) +
      theme_minimal() +
      theme(legend.position = "none")
    return(pca_plot)
  } else {
    # Get MDS
    mds_data <- cmdscale(dist(embeddings))

    # Plot MDS
    mds_plot <- ggplot(
      data = data.frame(
        mds_data,
        embeddings,
        samples = rownames(embeddings),
        stringsAsFactors = FALSE
      ),
      aes(
        x = as.numeric(mds_data[, 1]),
        y = as.numeric(mds_data[, 2]),
        label = samples,
        color = samples
      )
    ) +
      geom_point(size = point_size_) +
      geom_text_repel(color = "black", size = text_size_) +
      labs(
        title = basename(file_path),
        #title = "Principal Coordinate Analysis",
        color = NULL
      ) +
      scale_x_continuous(
        name = "Principal Coordinate 1"
      ) +
      scale_y_continuous(
        name = "Principal Coordinate 2"
      ) +
      theme_minimal() +
      theme(
        legend.position = "none"
      )
    return(mds_plot)
  }
}

embeddings_plot_list <- sapply(file_list, function(file_path) {
  # Read the file into a data frame
  embeddings <- read.table(
    file_path,
    header = TRUE,
    sep = ",",
    check.names = FALSE
  )
  return(get_plot(embeddings))
}, simplify = FALSE)

# Read the file into a data frame
embeddings_all <- read.table(
  "data/embeddings/all_unique.csv",
  header = TRUE,
  sep = ",",
  check.names = FALSE
)
embeddings_para <- read.table(
  "data/embeddings/figure_4H.csv",
  header = TRUE,
  sep = ",",
  check.names = FALSE
)
embeddings_gs <- read.table(
  "data/embeddings/figure_S4D.csv",
  header = TRUE,
  sep = ",",
  check.names = FALSE
)

# Set the first column as row names and remove it from the data frame
rownames(embeddings_all) <- gsub(
  ".bed|hg38 | STAR 1k| \\+ Centromere",
  "",
  embeddings_all[, 1]
)
embeddings_all <- embeddings_all[, -1]
rownames(embeddings_all) <- gsub(
  " 36bp List",
  " 36bp",
  rownames(embeddings_all)
)
rownames(embeddings_all) <- gsub(
  " 101bp List",
  " 101bp",
  rownames(embeddings_all)
)
rownames(embeddings_all) <- gsub(
  " List",
  " Blacklist",
  rownames(embeddings_all)
)
rownames(embeddings_all) <- gsub(
  "HS \\+ LM \\+ CM",
  "HS \\+ LM",
  rownames(embeddings_all)
)

rownames(embeddings_para) <- gsub(
  ".bed|hg38 | STAR 1k| \\+ Centromere",
  "",
  embeddings_para[, 1]
)
embeddings_para <- embeddings_para[, -1]
rownames(embeddings_para) <- gsub(
  " 36bp List",
  " 36bp",
  rownames(embeddings_para)
)
rownames(embeddings_para) <- gsub(
  " 101bp List",
  " 101bp",
  rownames(embeddings_para)
)
rownames(embeddings_para) <- gsub(
  " List",
  " Blacklist",
  rownames(embeddings_para)
)
rownames(embeddings_para) <- gsub(
  "HS \\+ LM \\+ CM",
  "HS \\+ LM",
  rownames(embeddings_para)
)

rownames(embeddings_gs) <- gsub(
  ".bed|hg38 | STAR 1k| \\+ Centromere",
  "",
  embeddings_gs[, 1]
)
embeddings_gs <- embeddings_gs[, -1]
rownames(embeddings_gs) <- gsub(
  " 36bp List",
  " 36bp",
  rownames(embeddings_gs)
)
rownames(embeddings_gs) <- gsub(
  " 101bp List",
  " 101bp",
  rownames(embeddings_gs)
)
rownames(embeddings_gs) <- gsub(
  " List",
  " Blacklist",
  rownames(embeddings_gs)
)
rownames(embeddings_gs) <- gsub(
  "HS \\+ LM \\+ CM",
  "HS \\+ LM",
  rownames(embeddings_gs)
)


embeddings_before_names <- rownames(rbind(embeddings_gs, embeddings_para))
embeddings_before <- embeddings_all[embeddings_before_names, ]

# Get PCA
#pca_data_before <- prcomp(embeddings_before)
#pca_data_before_df <- pca_data_before$x
#importance_before <- summary(pca_data_before)$importance

# Get MDS
mds_data_before <- cmdscale(dist(embeddings_before))

# Plot MDS
embeddings_plot_before <- ggplot(
  data = data.frame(
    mds_data_before,
    embeddings_before,
    samples = ifelse(
      rownames(embeddings_before) %in% rownames(embeddings_gs),
      rownames(embeddings_before),
      ""
    ),
    stringsAsFactors = FALSE
  ),
  aes(
    x = as.numeric(mds_data_before[, 1]),
    y = as.numeric(mds_data_before[, 2]),
    label = samples,
    color = samples
  )
) +
  geom_point(size = point_size, alpha = 1) +
  geom_text_repel(
    color = "black",
    size = text_size,
    max.overlaps = getOption(
      "ggrepel.max.overlaps",
      default = 100
    ),
  ) +
  labs(
    title = basename(file_path),
    color = NULL
  ) +
  scale_x_continuous(
    name = "Principal Coordinate 1"
  ) +
  scale_y_continuous(
    name = "Principal Coordinate 2"
  ) +
  theme_minimal() +
  theme(
    legend.position = "none"
  )

names(embeddings_plot_list) <- basename(names(embeddings_plot_list))

fig_2_expanded_lists <- c(
  "GitHub Blacklist",
  "Generated Blacklist",
  "High Signal",
  "Low Mappability",
  "HS + LM",
  "Kundaje Unified",
  "Nordin CUT&RUN",
  "GreyListChIP"
)
fig_2_expanded_embeddings <- embeddings_all[fig_2_expanded_lists, ]

fig_2_expanded_plot <- get_plot(
  fig_2_expanded_embeddings,
  preprocessed = TRUE
) +
  labs(title = NULL)
