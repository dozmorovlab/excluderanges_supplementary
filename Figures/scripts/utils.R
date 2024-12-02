#!/usr/local/bin/Rscript

# =============================================================================
# LIBARIES 
# =============================================================================
shhh <- suppressPackageStartupMessages # It's a library, so shhh!
shhh(library(GenomicRanges))
shhh(library(ggplot2))
shhh(library(rtracklayer))
shhh(library(dplyr))
shhh(library(ComplexHeatmap))
shhh(library(gridExtra))
shhh(library(gdata))
shhh(library(AnnotationHub))
shhh(library(TxDb.Hsapiens.UCSC.hg38.knownGene))

# =============================================================================
# TEST GRANGES OBJECTS
# =============================================================================
## create unit tests
test1 <- data.frame(
    chrom=c('chr1', 'chr1', 'chr1', 'chr1'),
    start=c(1, 7, 7, 12),
    end = c(5, 10,10,15)
    ) |> GenomicRanges::makeGRangesFromDataFrame() |> reduce()
test2 <- data.frame(
    chrom=c('chr1', 'chr1'),
    start=c(3, 20),
    end = c(14,25)
    ) |> GenomicRanges::makeGRangesFromDataFrame() |> reduce()
test_list <- list(test1, test2)

# =============================================================================
# DATA PREPARATION FUNCTIONS
# =============================================================================
# Read in bed file as GRanges object
# Output: GRanges object
#' @param bed_file path to bed file
#' @param merge_gaps boolean, fill in gaps or not. Used primarily within
#   centromeres, telomeres, and short arms
#' @param keepColumn boolean, keep metadata columns
#' @param bedGraph if TRUE, then the 4th column in bedGraph will be used as score
bed_to_GRanges <- function(bed_file, merge_gaps = FALSE, keepColumns = FALSE, bedGraph = FALSE){
    # read in file as dataframe
    bed <- readr::read_tsv(bed_file, col_names = FALSE)#[1:3]
    # change column names
    if (bedGraph) { # a 4-column bedGraph file
      colnames(bed) <- c('chr', 'start', 'end', 'score')
    } else { # First 3 columns of a BED file
      colnames(bed) <- c('chr', 'start', 'end')
    }
    # make GRanges object, remove overlap ranges with reduce()
    GRange <- bed |> makeGRangesFromDataFrame(keep.extra.columns = keepColumns, starts.in.df.are.0based = TRUE)
    GRange <- GRange |> sort() |> keepStandardChromosomes(pruning.mode="tidy") # |> reduce()
    if (merge_gaps == TRUE){
        GRange <- GRange |> reduce(drop.empty.ranges = TRUE, 10000000)
    }
    return(GRange)
}

# Make a list of GRanges objects from input bed files
# Output: list of GRanges objects
#' @param path_to_sets path to bed files to made into list of GRanges objects
#' @param remove_gaps boolean; whether or not to remove gaps from each GRanges object
#' @param path_to_gaps path to gap bed files; needed if remove_gaps = TRUE
#' @param keepColumn boolean, keep metadata columns
#' @param bedGraph if TRUE, then the 4th column in bedGraph will be used as score
make_excludeList <- function(
    path_to_sets, remove_gaps=TRUE, path_to_gaps, keepColumns = FALSE, bedGraph = FALSE
    ){
    # Read in excludable set files
    filenames <- list.files(path_to_sets, full.names=TRUE) |>
        stringr::str_sort(numeric=TRUE)
    # Make analysis list
    all_excludeGR_list <- lapply(
        filenames, 
        function(FILE){bed_to_GRanges(FILE, keepColumns = keepColumns, bedGraph = bedGraph)}
    )
    # Add names
    names(all_excludeGR_list) <- basename(filenames)
    # Remove intersection with gaps if remove_gaps = True
    if (remove_gaps == T){
        # get gap file names
        gap_files <- list.files(path_to_gaps, full.names=T)
        # create list of gap files as GRanges objects
        gap_list <- lapply(
            gap_files, 
            function(SET){bed_to_GRanges(SET)}
        )
        # fix gaps in gaps
        gap_list <- lapply(
            gap_list, 
            function(set){reduce(set, drop.empty.ranges = TRUE, 10000000)}
        )
        # combine all gaps into one GRanges object
        all_gaps <- do.call("c", gap_list) |> sort()
        # # Get regions overlapping gaps
        # removed_gaps <- lapply(
        #     all_excludeGR_list,
        #     function(SET){subsetByOverlaps(SET, all_gaps)}
        # )
        # # Remove regions overlapping gaps
        # all_excludeGR_list <- mapply(
        #     function(SET, GAPS){GenomicRanges::setdiff(SET, GAPS)},
        #     all_excludeGR_list,
        #     removed_gaps
        # )
        # Remove regions overlapping gaps
        # https://www.biostars.org/p/263214/
        for (i in seq_along(all_excludeGR_list)) {
          # Check if there are any overlaps with gaps, do nothing if none
          if (length(queryHits(findOverlaps(all_excludeGR_list[[i]], all_gaps, type="any"))) > 0) {
            all_excludeGR_list[[i]] <- all_excludeGR_list[[i]][-queryHits(findOverlaps(all_excludeGR_list[[i]], all_gaps, type="any")),] 
          }
        }
    }
    return(all_excludeGR_list)
}

# Make vector of names from input bed files, used for plotting
# Output: vector of names
#' @param path_to_sets path to bed files
#' @param remove_prefix boolean, whether or not to remove genome prefix
#' @param remove_version boolean, whether or not to remove version prefix
make_names <- function(
    path_to_sets, remove_prefix = TRUE, remove_version = FALSE
    ){
    # Read in excludable set files
    Names <- list.files(path_to_sets, full.names=F) |>
        stringr::str_sort(numeric=TRUE)
    # Make short names from filenames, useful for plotting
    Names <- lapply(Names, function(X) sub('.bed', '', X))
    # remove prefix, e.g. 'hg38.'
    if (remove_prefix == TRUE){
        Names <- lapply(
            Names, function(X){
                sub(paste0(genome, '.'), '', X)
            }
        ) |> unlist()
    }
    # remove version number, must exist at beginning of filename, e.g.
    #   1.something.bed, 2.something.bed, etc
    if (remove_version == TRUE){
        Names <- lapply(
            Names, function(X){
                sub('.*\\.', '', X)
            }
        ) |> unlist()
    }
    return(Names |> as.character())
}

# Download GRanges objects from Google Drive
# Output: GRanges object
#' @param gID 33 character google drive ID, e.g. 1EyUiDqVLjVmofRl3cGGI2Mj1swRBEfQK
#' @param output save file downloaded file with this name
get_gdrive_data <- function(gID, output){
    first_part <- "https://drive.google.com/uc?export=download&id="
    # download file 
    if (!file.exists(file.path(output))){
        download.file(
            url = paste0(first_part, gID),
            destfile = output, method = 'auto'
        )
    }
    # read in bedfile
    SET <- bed_to_GRanges(output)
    return(SET)
}

# Download GRanges objects annotation hub
# Output: GRanges object
#' @param ID annotation hub ID
#' @param output save file downloaded file with this name
get_data <- function(ID, output){
    # download file if it doesn't exist
    if (!file.exists(output)){
        Set <- AnnotationHub()[[ID]] %>% 
            sort() %>% 
            keepStandardChromosomes(pruning.mode="tidy")
        export.bed(Set, con=output)
        }
}

# =============================================================================
# CALCULATION FUNCTIONS
# =============================================================================
# function to calculate percent of grA overlapping grB by width
#' @param grA GRanges object
#' @param grB GRanges object
percent_covers_width <- function(grA, grB){
    # get size of set a
    grA_size <- grA %>% width() %>% sum()
    # get intersects of set a and b
    intersects <- GenomicRanges::intersect(grA, grB, ignore.strand=TRUE) %>% reduce %>% width %>% sum
    # get percent covers value
    percent <- 100 * intersects / grA_size
    # if 0 < value < 1 then write "< 1.0"
    if ((percent > 0) & (percent < 1)){percent <- as.character('< 1.0')}
    else{percent <- as.character(round(100*intersects/grA_size, 2))}
    return(percent)
}

# function to calculate percent of grA overlapping grB by count
#' @param grA GRanges object
#' @param grB GRanges object
percent_covers_count <- function(grA, grB){
    # get number of regions in a
    grA_size <- grA %>% length 
    # get number of intersects between a and b
    intersects <- subsetByOverlaps(grA, grB) %>% length 
    # calculate percent covers count
    percent <- as.character(round(100*intersects/grA_size, 2))
    # write output as: X% (N-intersects/N-regions-in-a)
    report <- paste0(percent, '% (', intersects, '/',grA_size,')')
    return(report)
}

# Get width of intersection between grA and grB
#' @param grA GRanges object
#' @param grB GRanges object
intersection_width <- function(grA, grB){
    Intersect <- GenomicRanges::intersect(grA, grB) |> width() |> sum()
    return(Intersect)
}

# Get number of grA regions intersecting grB regions
#' @param grA GRanges object
#' @param grB GRanges object
intersection_count <- function(grA, grB){
    Intersect <- subsetByOverlaps(grA, grB) |> length()
    return(Intersect)
}

# Get overlap coefficients (in terms of width) between two sets, i.e. what 
#   proportion of the smaller set is contained within the larger?
#' @param grA GRanges object
#' @param grB GRanges object
overlap_coeff_width <- function(grA, grB) {
    # get smaller set by count
    if (length(grA) <= length(grB)){
        smaller_set <- grA
        larger_set <- grB
    }else{
        smaller_set <- grB
        larger_set <- grA
    }
    # get intersection of width
    IW <- GenomicRanges::intersect(grA, grB, ignore.strand=TRUE) %>% reduce %>% width %>% sum
    # get width of smaller set
    MW <- smaller_set %>% width %>% sum
    # calculate coefficient
    overlap_coeff <- IW/MW
    return(overlap_coeff)
}

# Get overlap coefficients (in terms of number) between two sets, i.e. what
#   proportion of regions of the smaller set intersect with the larger?
#' @param grA GRanges object
#' @param grB GRanges object
overlap_coeff_count <- function(grA, grB) {
    # get smaller set by count
    if (length(grA) <= length(grB)){
        smaller_set <- grA
        larger_set <- grB
    }else{
        smaller_set <- grB
        larger_set <- grA
    }
    # get intersection of counts
    IC <- subsetByOverlaps(smaller_set, larger_set) %>% length()
    # get smaller count
    MC <- smaller_set %>% length
    # calculate coefficient
    coeff <- IC/MC
    return(coeff)
}

# Get jaccard width between two sets.
#' @param grA GRanges object
#' @param grB GRanges object
jaccard_width <- function(grA, grB) {
    # intersection of a and b
    IW_1 <- GenomicRanges::intersect(grA, grB, ignore.strand=TRUE) %>% reduce %>% width %>% sum
    IW_2 <- GenomicRanges::intersect(grB, grA, ignore.strand=TRUE) %>% reduce %>% width %>% sum
    IW <- mean(IW_1, IW_2)

    # union of a and b
    UW <- GenomicRanges::union(grA, grB, ignore.strand=TRUE) %>% GenomicRanges::reduce() %>% GenomicRanges::width() %>% sum
    # calculate coefficient
    coeff <- IW / UW
    return(coeff)
}

# Get jaccard count between two sets.
#' @param grA GRanges object
#' @param grB GRanges object
jaccard_count <- function(grA, grB) {

    # intersection of a and b
    IC_1 <- subsetByOverlaps(grA, grB) %>% length()
    IC_2 <- subsetByOverlaps(grB, grA) %>% length()
    IC <- mean(IC_1, IC_2)

    # union of a and b
    UC <- GenomicRanges::union(grA, grB, ignore.strand = TRUE) %>% length()

    # calculate coefficient
    coeff <- IC / UC
    return(coeff)
}

# Get size of reference genome
#' @param genome_id reference genome, one of: "hg38"
#' @param only_autosomal boolean, whether to remove X,Y, and M chromosomes or not
#' @param include_chrs list of chromosome names to include
get_genome_size <- function(genome_id="hg38", only_autosomal=TRUE, include_chr=NULL){

    # get genome info
    chrom_data <- GenomeInfoDb::getChromInfoFromUCSC(
                        genome = genome_id,
                        assembled.molecules.only = TRUE
    )
    if (only_autosomal == TRUE){
        # remove chromosomes M, X, Y, as not all excludable sets have them
        chrom_data <- chrom_data[chrom_data$chrom != 'chrM', ]
        chrom_data <- chrom_data[chrom_data$chrom != 'chrX', ]
        chrom_data <- chrom_data[chrom_data$chrom != 'chrY', ]
    } else if (! is.null(include_chr)) {
        chrom_data <- chrom_data[chrom_data$chrom %in% include_chr, ]
    }
    # sum length of chromosomes
    genome_length <- sum(chrom_data$size)
    return(genome_length)
}

# Get size of reference genome chromosomes
#' @param genome_id reference genome, one of: "hg38"
#' @param only_autosomal boolean, whether to remove X,Y, and M chromosomes or not
get_chrom_sizes <- function(genome_id, only_autosomal=TRUE){
    if (genome_id == 'hg38'){
        # get genome info
        chrom_data <- GenomeInfoDb::getChromInfoFromUCSC(genome = genome_id, assembled.molecules.only = TRUE)
    }
    if (only_autosomal==TRUE){
        # remove chromosomes M, X, Y, as not all excludable sets have them
        chrom_data <- chrom_data[chrom_data$chrom != 'chrM', ]
        chrom_data <- chrom_data[chrom_data$chrom != 'chrX', ]
        chrom_data <- chrom_data[chrom_data$chrom != 'chrY', ]
    }
    chrom_data <- chrom_data %>% select(chrom, size)
    chrom_data$start <- 1
    colnames(chrom_data) <- c('chrom', 'end', 'start')
    chrom_data <- chrom_data %>% GenomicRanges::makeGRangesFromDataFrame()
    return(chrom_data)
}


# Calculate forbes coefficient
#' @param grA GRanges object
#' @param grBGRanges object
#' @param genome_size size of genome corresponding to grA and grB; default is hg38 size
#' @param only_autosomal boolean, whether to remove X,Y, and M chromosomes or not; default TRUE
forbes_coeff_width <- function(grA, grB, genome_size=2875001522, only_autosomal=TRUE){
    # remove X,Y, and M chromosomes, as not all excludable sets have them
    if (only_autosomal==TRUE){
        grA <- dropSeqlevels(grA, c('chrX', 'chrY', 'chrM'), pruning.mode='tidy')
        grB <- dropSeqlevels(grB, c('chrX', 'chrY', 'chrM'), pruning.mode='tidy')
    }
    # get total width of intersection
    Intersection_width <- GenomicRanges::intersect(grA, grB, ignore.strand=TRUE) %>% reduce %>% width %>% sum
    # get total width of A
    A_cardinality <- grA %>% reduce %>% width %>% sum
    # get total width of B
    B_cardinality <- grB %>% reduce %>% width %>% sum
    forbes_coeff = (genome_size*Intersection_width)/(as.numeric(A_cardinality)*as.numeric(B_cardinality))
    return(forbes_coeff)
}

# Calculate forbes coefficient for each chromosome, then take mean
#' @param grA GRanges object
#' @param grB GRanges object
#' @param genome_name name of genome corresponding to grA and grB
#' @param only_autosomal boolean, whether to remove X,Y, and M chromosomes or not; default TRUE
forbes_coeff_width_by_chrom <- function(grA, grB, genome_name='hg38', only_autosomal=TRUE){
    # get intersection lengths
    intersections <- GenomicRanges::intersect(grA, grB, ignore.strand=TRUE)
    int_df <- intersections %>% as.data.frame
    intersections <- intersections %>% split(seqnames(intersections)) %>% width %>% sum
    intersections <- intersections[which(intersections != 0)]
    # get lengths of gra and grb
    grA_df <- grA %>% as.data.frame
    grB_df <- grB %>% as.data.frame
    grA_df <- grA_df[which(grA_df$seqnames %in% int_df$seqnames),]
    grB_df <- grB_df[which(grB_df$seqnames %in% int_df$seqnames),]
    grA_chroms <- grA_df %>% makeGRangesFromDataFrame
    grB_chroms <- grB_df %>% makeGRangesFromDataFrame
    grA_chroms <- grA_chroms %>% split(seqnames(grA_chroms)) %>% width %>% sum
    grB_chroms <- grB_chroms %>% split(seqnames(grB_chroms)) %>% width %>% sum
    grA_chroms <- grA_chroms[which(grA_chroms != 0)]
    grB_chroms <- grB_chroms[which(grB_chroms != 0)]
    # Get reference genome chromosome lengths
    Chroms <- get_chrom_sizes(genome_name, only_autosomal)
    chroms_df <- Chroms %>% as.data.frame
    chroms_df <- chroms_df[which(chroms_df$seqnames %in% int_df$seqnames),]
    Chroms <- chroms_df %>% makeGRangesFromDataFrame
    Chroms <- Chroms %>% split(seqnames(Chroms)) %>% width %>% sum
    Chroms <- Chroms[which(Chroms != 0)]
    # Numbers are too large for R, work with log10 numbers instead
    #   original eq: forbes_coeff <- (genome_size * Intersection_width) / (A_cardinality * B_cardinality)
    log10_forbes <- (log10(Chroms) + log10(intersections)) - ((log10(grA_chroms)) + log10(grB_chroms))
    # convert back and take mean
    forbes_coeff = 10^log10_forbes %>% mean
    return(forbes_coeff)
}

# Cluster regions based on a specified width threshold
#' @author Brydon Wall
#' @param gr GRanges object
#' @param merge_dist merge regions with a neighbor distance <= than this
cluster_gr <- function(gr, merge_dist) {
    # # Invert gr to get all neighbors as regions
    # # The last regions get erased, but we add back later
    # inverted_gr <- GenomicRanges::gaps(gr)
    # # Keep only neighbors that are > merge_dist
    # subset_gr <- inverted_gr[GenomicRanges::width(inverted_gr) > merge_dist]
    # # Invert again to get original regions, now merged
    # # Take the union with the original gr to add the last regions back
    # out_gr <- GenomicRanges::union(GenomicRanges::gaps(subset_gr), gr) |> reduce()
    # return(out_gr)
  # Ensure the GRanges object is sorted
  gr <- sort(gr)
  # Use the reduce function to merge ranges
  reduced_gr <- reduce(gr, min.gapwidth = merge_dist)
  return(reduced_gr)
}

# ===============================================================================
# PLOTTING FUNCTIONS
# ===============================================================================
# Create matrix of overlap coefficient for each pairwise combination
#   in input list.
# Output: matrix of overlap measured by @param fun
#' @param LIST list of GRanges objects
#' @param names vector of strings describing LIST
#' @param fun function to measure overlap, e.g. jaccard_count
#' @param Vend number of version, e.g. if you have 10 versions of each file,
#   use Vend = 10
overlap_matrix <- function(LIST, Names, fun, Vend = 1, ...){
    if (Vend == 1){
        mtx_to_plot <- matrix(data = 0, nrow = length(LIST), ncol = length(LIST))
        for (i in 1:length(LIST)){
            for (j in 1:length(LIST)){
                if (i == j){ # remove diagonal values using NA instead of 1
                    #mtx_to_plot[i, j] <- 1.00
                    mtx_to_plot[i, j] <- NA
                }
                else {
                    mtx_to_plot[i, j] <- suppressWarnings(do.call(fun, list(LIST[[i]], LIST[[j]], ...)))
                    if (is.nan(mtx_to_plot[i,j])){ mtx_to_plot[i,j] = 0 }
                    mtx_to_plot[j, i] <- mtx_to_plot[i,j]
                }
            }
        }
        rownames(mtx_to_plot) <- Names
        colnames(mtx_to_plot) <- Names
        return(mtx_to_plot)
    }
    # Version numbers
    versions <- seq(1, Vend)
    subLists <- lapply(
        versions,
        function(v){
            # Subset list to specific version
            subList <- LIST[which(grepl(paste0('^',v,'\\.'), Names))]
            # Get list names
            subNames <- Names[which(grepl(paste0('^',v,'\\.'), Names))]
            # Remove version number, for column/row labelling
            subNames <- sapply(
                subNames, function(X){sub(paste0('^',v,'\\.'), '', X)}
            )
            # Correlation matrix, empty
            mtx_to_plot <- matrix(
                data = 0, nrow = length(subList), ncol = length(subList)
            )
            for (i in 2:length(subList)){
                for (j in 1:length(subList)){
                    if (i == j){ # remove diagonal values using NA instead of 1
                        #mtx_to_plot[i, j] <- 1.00
                        mtx_to_plot[i, j] <- NA
                    }
                    else {
                        mtx_to_plot[i, j] <- suppressWarnings(fun(subList[[i]], subList[[j]]))
                        if (is.nan(mtx_to_plot[i,j])){ mtx_to_plot[i,j] = 0 }
                        mtx_to_plot[j, i] <- mtx_to_plot[i,j]
                    }
                }
            }
            colnames(mtx_to_plot) <- rownames(mtx_to_plot) <- subNames
            return(mtx_to_plot)
        }
    )
    # Get mean of matrices in list
    mtx_to_plot <- Reduce("+", subLists) / length(subLists)
    return(mtx_to_plot)
}

# Visualize input matrix as a heatmap.
#' @param mtx matrix 
#' @param cell_size int, size of square cells
#' @param Title string, title of plot; default "heatmap"
#' @param Fontsize double, size of font proportional to cell size; default 0.8
Heatmap <- function(
    mtx, cell_size, Title="heatmap", Fontsize=0.8,
    show_col = FALSE, legend_side = "bottom", cluster_rows = FALSE,
    legend_direction = "horizontal", show_center = TRUE
){
    # order matrix
    mtx <- mtx[order(rowSums(mtx, na.rm=T), decreasing=TRUE),
                     order(colSums(mtx, na.rm=T), decreasing=TRUE)]
    # make quantiles from matrix values
    quants <- mtx %>% quantile(probs = seq(0, 1, 0.1), na.rm = T) %>% unique() %>% sort()
    # make color vector from quantiles
    Colors <- colorRampPalette(c('floralwhite',
                                 'mistyrose',
                                 'red'))(length(quants))
    # make color vector for values
    col_val <- circlize::colorRamp2(quants, Colors, transparency=0)
    # Lower half of heatmap
    H1 <- ComplexHeatmap::Heatmap(mtx, rect_gp = gpar(type = 'none'), col = col_val,
                            width = ncol(mtx)*unit(cell_size, "mm"),
                            height = nrow(mtx)*unit(cell_size, "mm"),
                            # no clustering, has been done already
                            cluster_rows = cluster_rows, cluster_columns = FALSE,
                            # no row or column names
                            show_row_names = FALSE, show_column_names = FALSE,
                            # legend parameters
                            heatmap_legend_param = list(title='Overlap',
                                                        legend_direction = legend_direction,
                                                        #legend_width=unit(100, 'mm')
                                                        legend_width=ncol(mtx)*unit(cell_size, "mm")
                                                        ),
                            # cell formatting
                            cell_fun = function(j,i,x,y,w,h,fill) {
                                # if diagonal, set rownames
                                if (i == j) {
                                    grid.rect(x,y,w,h, gp = gpar(fill='gray71', col='black'))
                                    if (show_center){
                                        grid.text(paste0(colnames(mtx)[i], '\n', round(sum(mtx[i,], na.rm = TRUE), 2)),
                                                gp = gpar(fontsize = Fontsize*cell_size, col = 'black', fontface = 'bold'), x, y)
                                    } else {
                                        grid.text(paste0(round(sum(mtx[i,], na.rm = TRUE), 2)),
                                                gp = gpar(fontsize = Fontsize*cell_size, col = 'black', fontface = 'bold'), x, y)
                                    }

                                }
                                # otherwise fill with appropriate color
                                if (i > j){grid.rect(x,y,w,h, gp = gpar(fill=col_val(mtx[i,j]), col='black'))}
                            }
    )
    # Upper half of heatmap
    H2 <- ComplexHeatmap::Heatmap(mtx, rect_gp = gpar(type = 'none'),
                            column_title = Title,
                            width = ncol(mtx)*unit(cell_size, "mm"),
                            height = nrow(mtx)*unit(cell_size, "mm"),
                            # no clustering
                            cluster_rows = FALSE, cluster_columns = FALSE,
                            # show only row names
                            show_row_names = TRUE, show_column_names = show_col,
                            # no legend
                            show_heatmap_legend = FALSE,
                            # cell formatting
                            cell_fun = function(j,i,x,y,w,h,fill) {
                                # if diagonal, set blank
                                if (i == j) {grid.text(NULL)}
                                else if ((i<j) & (mtx[i,j] > 0.001)){
                                    # make square grey
                                    grid.rect(x,y,w,h,gp = gpar(fill='white', col='gray80'))
                                    # show overlap value
                                    grid.text(sprintf("%.2f", mtx[i,j]), x,y, gp = gpar(fontsize=Fontsize*cell_size, col='black'))
                                }else if ((i < j) & (mtx[i,j] < 0.01) & (mtx[i,j] > 0)){
                                    # make square grey
                                    grid.rect(x,y,w,h,gp = gpar(fill='white', col='gray80'))
                                    # show overlap value
                                    grid.text('<0.01', x,y, gp = gpar(fontsize=Fontsize*cell_size, col='black'))
                                }else if ((i < j) & (mtx[i,j] == 0)){
                                    # make square grey
                                    grid.rect(x,y,w,h,gp = gpar(fill='white', col='gray80'))
                                    # show overlap value
                                    grid.text('0.000', x,y, gp = gpar(fontsize=Fontsize*cell_size, col='black'))
                                }
                            }
    )
    # Visualize heatmap
    H <-draw(H1 + H2,
         heatmap_legend_side = legend_side,
         ht_gap = unit(-(nrow(mtx)*cell_size), 'mm')
    )
    return(H)
}

# get mean value of input mtx
#' @param mtx matrix
#' @param mtx_desc string describing mtx
get_mean_overlap <- function(mtx, mtx_desc="input matrix"){
    # Find mean of matrix, ignore diagonal
    Mean <- mean(mtx, na.rm=TRUE) %>% round(4)
    # Find standard deviate of matrix, ignore diagonal
    Std <- sd(mtx, na.rm=TRUE) %>% round(4)
    # verbose report
    Report <- paste('The mean value in ', mtx_desc, ' is:', Mean, '+/-', Std)
    return(Report)
}

# Make bar plot of X and Y vectors
#' @param Dframe input dataframe
#' @param X character vector of labels
#' @param Y numeric vector of values
#' @param Y_title title for Y-axis
#' @param Title title for plot
#' @param short boolean, whether to use short label names
#' @param Sort boolean, whether or not to sort by values
#' @param remove_standard boolean, whether to remove standard set from plot
#' @param Legend position of legend
bar_plot <- function(
    Dframe, X, Y, FILL = NULL, X_title = NULL,
    Y_title = NULL, Title = NULL, Legend='right'
){
    PLOT <- ggplot(
        Dframe, aes(x = X, y = Y, fill = FILL)
    ) +
    geom_col(color = 'black') +
    labs(
        y = X_title,
        x = Y_title,
        title = Title
    ) + 
    theme_classic() +
    theme(
        axis.text.x = element_text(size=Size, angle = 0, hjust = .5), 
        legend.position = Legend,
        axis.title.x = element_text(size=Size), 
        axis.text.y = element_text(size=Size, angle = 0, vjust = 0, hjust = 0), 
        axis.title.y = element_text(size=Size), 
        plot.title  = element_text(size=Size+2, hjust=0.5)
    ) 
    return(PLOT)
}

# Box plot function
#' @param Dframe input dataframe
#' @param X Dframe column to plot on X axis
#' @param Y Dframe column to plot on Y axis
#' @param FILL group by this vector
#' @param X_title x axis title
#' @param Y_title y axis title
#' @param Title title for plot
box_plot <- function(
        Dframe = NULL, X, Y, FILL = NULL,
        X_title = NULL, Y_title = NULL, Title = NULL
    ){
    PLOT <- ggplot(Dframe, aes(x = X, y = Y, fill = FILL)) +
        geom_boxplot(outlier.color = NA) +
        geom_jitter(height=0) +
        labs(
            x = X_title,
            y = Y_title, 
            title = Title
        ) + 
        theme_classic() +
        theme(
            legend.position = 'none',
            axis.text.x = element_text(
                size=Size, angle = 45, hjust = 1
            ), 
            axis.title.x = element_text(
                size=Size
            ), 
            axis.text.y = element_text(
                size=Size, angle = 0, vjust = -.5, hjust = 1
            ), 
            plot.title  = element_text(
                size=Size, hjust=0.5
            )
        ) 
    return(PLOT)
}

# Show width distributions of regions widths in each GRanges objects in @param excludeGR_list
# Output: width distribution plot
#' @param excludeGR_list list of GRanges objects
#' @param Names vector assigning names to excludeGR_list
width_distribution <- function(excludeGR_list, Names){
    # find the width distribution of each excludable set
    Width = lapply(excludeGR_list, function(x) as.integer(width(x)))
    Width <- Width |> unlist() |> as.integer()
    # each distribution needs a name assigned to it
    Source <- mapply(
        function(name, grange_object){
            rep(name, length(grange_object))
        },
        Names, excludeGR_list
    )
    Source <- Source |> unlist() |> as.character()
    # make width distribution dataframe, all files
    width_dist_data <- data.frame(Width, Source)
    # Remove versions from names
    Unique_names <- sapply(
        Names, 
        function(NAME){sub('^.*\\.','', NAME)}
        ) |> unique()
    PLOTS <- lapply(
        Unique_names,
        function(NAME){
            # Get filter by nBams
            Plot_dataframe <- width_dist_data |> 
                filter(
                    grepl(paste0('\\.', NAME, '$'), Source)
                )

            # plot distribution of region widths
            WD <- ggplot(
                    Plot_dataframe, 
                    aes(x = log10(Width), fill = Source)
                ) +
                geom_density(color = 'black', alpha=0.4) +
                labs(
                    x = 'log10 region widths (bases)',
                    y='',
                    title = NAME
                ) +
                theme_classic() +
                theme(
                    legend.position = 'none',
                    axis.text.x = element_text(
                        size=Size, angle = 0, hjust = 0.5
                    ),
                    axis.title.x = element_text(
                        size=Size
                    ),
                    axis.text.y = element_text(
                        size=Size, angle = 0, vjust = .5, hjust = 1
                    ),
                    plot.title  = element_text(
                        size=Size+2, hjust=0.5
                    )
                )
            return(WD)
        }
    )
    # Put plots in a list
    plots <- lapply(
        seq(1, length(PLOTS)),
        function(X){PLOTS[[X]]}
    )
    # Call list into a grid
    Width_distributions <- do.call("grid.arrange", c(plots, ncol = 1))
    return(Width_distributions)
}

# Show width distributions of regions widths in each GRanges objects in @param excludeGR_list
# Output: width distribution plot
#' @param excludeGR_list list of GRanges objects
#' @param names vector assigning names to excludeGR_list
#' @param short boolean, should set names be shortened?
width_distribution1 <- function(
  excludeGR_list, Names, short = TRUE, 
  Legend = 'right', Title = 'Width distribution'
  ){
    # make vector of abbreviated names
    if (short == TRUE){
      name_short = unlist(lapply(Names, function(name) name_dict[name]))
    } else {
      name_short = Names
    }
    # find the width distribution of each excludable set
    Width = lapply(excludeGR_list, function(x) as.integer(width(x)))
    Width <- Width %>% unlist() %>% as.integer()
    # each distribution needs a name assigned to it
    Source <- mapply(function(name, grange_object) rep(name, length(grange_object)), name_short, excludeGR_list)
    Source <- Source %>% unlist() %>% as.character()
    # Find the median of each width distribution
    Median <- lapply(excludeGR_list, function(set) rep(median(width(set)), length(set)))
    Median <- Median %>% unlist() %>% as.integer()
    # make width distribution dataframe
    width_dist_data <- data.frame(Width, Source, Median)
    width_dist_data <- width_dist_data[order(width_dist_data$Median),]
    # plot distribution of region widths
    WD <- ggplot(
      width_dist_data,
      aes(x = log10(Width), y = reorder(Source, Median), fill = Source)
      ) +
      ggridges::geom_density_ridges() +
      scale_fill_manual(
        values = Set_palette,
        labels = rev_dict[names(Set_palette)],
        breaks = rev(width_dist_data$Source %>% unique())
      ) +
      labs(
        x = 'log10 region widths (bases)',
        y='',
        title = Title
      ) +
      theme_classic() +
      theme(
        legend.position = Legend,
        axis.text.x = element_text(size=Size, angle = 0, hjust = 0.5),
        axis.title.x = element_text(size=Size),
        axis.text.y = element_text(size=Size, angle = 0, vjust = .5, hjust = 1),
        plot.title  = element_text(size=Size+2, hjust=0.5)
      )
    return(WD)
}

# Merge all cetromere ranges per chr
# Output: path_to_cent merged
#' @author Brydon Wall
#' @param path_to_cent a path to the centromere BED file
process_centromeres <- function(path_to_cent) {
  
  # Extract the directory and filename
  directory <- dirname(path_to_cent)
  filename <- basename(path_to_cent)

  # Add "_merged" to the filename
  new_filename <- paste0(
    tools::file_path_sans_ext(filename),
    "_merged.",
    tools::file_ext(filename)
  )

  # Create the new file path
  new_filepath <- file.path(directory, new_filename)

  gr <- bed_to_GRanges(path_to_cent)

  # Find the minimum start and maximum end positions within each chromosome
  min_start <- tapply(start(gr), seqnames(gr), min)
  max_end <- tapply(end(gr), seqnames(gr), max)

  # Create a new GRanges object with the combined ranges
  combined_gr <- GRanges(
  seqnames = names(min_start),
  ranges = IRanges(start = min_start, end = max_end)
  )

  # Save the new BED
  save_bed(combined_gr, new_filepath)

  file.remove(path_to_cent)
}

# Downloads gap regions for T2T, hg38, hg19, mm39, and mm10.
download_data <- function(){
# Create data dir
  dir.create(file.path('.', 'data'), recursive = TRUE, showWarnings = FALSE)
# Get T2T gaps --------------------------------------------------------------
  dir_gaps <- file.path('.','data','T2T_gaps')
  if (!dir.exists(dir_gaps)) {
    dir.create(dir_gaps, recursive=T, showWarnings=FALSE)
    T2T.telomere <- get_data('AH107350', file.path(dir_gaps, 'T2T.telomere.gap.bed'))
    T2T.cen_mask <- get_data('AH107349', file.path(dir_gaps, 'T2T.cen_mask.gap.bed'))
# T2T cenmask is combination of short arms and telomeres. Separate.
# Process short arms
    T2T.cen_mask <- bed_to_GRanges(file.path(dir_gaps, 'T2T.cen_mask.gap.bed'))
    T2T.short_arms <- as.data.frame(T2T.cen_mask)[-c(1:12, 16:20, 23), ] # keep only short arm chroms
    T2T.short_arms[1,]$end = 15500000
    T2T.short_arms[2,]$end = 10000000
    T2T.short_arms[3,]$end = 16500000
    T2T.short_arms[4,]$end = 11000000
    T2T.short_arms[5,]$end = 12700000
# Save short arms
    export.bed(T2T.short_arms, con = file.path(dir_gaps, 'T2T.short_arms.gap.bed'))
# Process centromeres
    T2T.centromeres <- T2T.cen_mask %>% as.data.frame
    T2T.centromeres[13,]$start = 15500000
    T2T.centromeres[14,]$start = 10000000
    T2T.centromeres[15,]$start = 16500000
    T2T.centromeres[21,]$start = 11000000
    T2T.centromeres[22,]$start = 12700000
# Save centromeres
    export.bed(T2T.centromeres, con = file.path(dir_gaps, 'T2T.centromeres.gap.bed'))

process_centromeres(file.path(dir_gaps, 'T2T.centromeres.gap.bed'))

# Remove cenmask file
    file.remove(file.path(dir_gaps, 'T2T.cen_mask.gap.bed'))
  }
  
# Get hg38 gaps ------------------------------------------------------------
  dir_gaps <- file.path('.','data','hg38_gaps')
  if (!dir.exists(dir_gaps)) {
    dir.create(dir_gaps, recursive=T, showWarnings=FALSE)
    hg38.centromeres <- get_data('AH107354', file.path(dir_gaps, 'hg38.UCSC.centromere.gap.bed'))

    process_centromeres(file.path(dir_gaps, 'hg38.UCSC.centromere.gap.bed'))

    hg38.telomeres   <- get_data('AH107355', file.path(dir_gaps, 'hg38.UCSC.telomere.gap.bed'))
    hg38.short_arms  <- get_data('AH107356', file.path(dir_gaps, 'hg38.UCSC.short_arms.gap.bed'))
  }
  
# Get hg19 gaps ------------------------------------------------------------
  dir_gaps <- file.path('.','data','hg19_gaps')
  if (!dir.exists(dir_gaps)) {
    dir.create(dir_gaps, recursive=T, showWarnings=FALSE)
    hg19.centromere <- get_data('AH107360', file.path(dir_gaps, 'hg19.UCSC.centromere.gap.bed'))

    process_centromeres(file.path(dir_gaps, 'hg19.UCSC.centromere.gap.bed'))

    hg19.short_arm  <- get_data('AH107362', file.path(dir_gaps, 'hg19.UCSC.short_arm.gap.bed'))
    hg19.telomere   <- get_data('AH107361', file.path(dir_gaps, 'hg19.UCSC.telomere.gap.bed'))
  }
  
# Get mm39 gaps ------------------------------------------------------------
  dir_gaps <- file.path('.','data','mm39_gaps')
  if (!dir.exists(dir_gaps)) {
    dir.create(dir_gaps, recursive=T, showWarnings=FALSE)
    mm39.centromere <- get_data('AH107367', file.path(dir_gaps,'mm39.UCSC.centromere.gap.bed'))

    process_centromeres(file.path(dir_gaps,'mm39.UCSC.centromere.gap.bed'))

    mm39.short_arm <- get_data('AH107369', file.path(dir_gaps, 'mm39.UCSC.short_arm.gap.bed'))
    mm39.telomere <- get_data('AH107368', file.path(dir_gaps, 'mm39.UCSC.telomere.gap.bed'))
  }
  
# Get mm10 gaps ------------------------------------------------------------
  dir_gaps <- file.path('.','data','mm10_gaps')
  if (!dir.exists(dir_gaps)) {
    dir.create(dir_gaps, recursive=T, showWarnings=FALSE)
    mm10.centromere <- get_data('AH107372', file.path(dir_gaps, 'mm10.UCSC.centromere.gap.bed'))

    process_centromeres(file.path(dir_gaps, 'mm10.UCSC.centromere.gap.bed'))

    mm10.short_arm  <- get_data('AH107374', file.path(dir_gaps, 'mm10.UCSC.short_arm.gap.bed'))
    mm10.telomere   <- get_data('AH107373', file.path(dir_gaps, 'mm10.UCSC.telomere.gap.bed'))
  }
}

# Fetch genome knownGene data
# Output: genome_exons.bed, genome_genes.bed, genome_cds.bed
#' @author Brydon Wall
#' @param exclude_chr a vector with chrs to exclude
#' @param exclude_alts whether to exclude alts (e.g. _alt or _random or _fix or chrUn)
get_genome_knownGenes <- function(
    exclude_chr = c("chrX", "chrY", "chrM"),
    exclude_alts = TRUE,
    genome = "hg38",
    gene_data = TxDb.Hsapiens.UCSC.hg38.knownGene
) {

  filter_chromosomes <- function(dataset, exclude_chr) {
    chrs <- unique(as.vector(seqnames(dataset)))
    if (exclude_alts) {
      alts <- chrs[grepl("_alt|_random|_fix|chrUn", chrs)]
    } else {
      alts <- c()
    }
    return(dropSeqlevels(dataset, c(exclude_chr, alts), pruning.mode = "tidy"))
  }

  # Check for files
  exons_file <- file.path(".", "data", "exons", paste0(genome, "_exons.bed"))
  genes_file <- file.path(".", "data", "genes", paste0(genome, "_genes.bed"))
  cds_file <- file.path(".", "data", "cds", paste0(genome, "_cds.bed"))

  # Exons
  if (!file.exists(exons_file)) {
    dir.create(dirname(exons_file))
    genome_exons <- exons(gene_data)
    genome_exons <- filter_chromosomes(genome_exons, exclude_chr)
    save_bed(genome_exons, exons_file)
  }

  # Genes
  if (!file.exists(genes_file)) {
    dir.create(dirname(genes_file))
    genome_genes <- genes(gene_data)
    seqnames(genome_genes)
    genome_genes <- filter_chromosomes(genome_genes, exclude_chr)
    save_bed(sort(genome_genes), genes_file)
  }

  # CDS
  if (!file.exists(cds_file)) {
    dir.create(dirname(cds_file))
    genome_cds <- cds(gene_data)
    genome_cds <- filter_chromosomes(genome_cds, exclude_chr)
    save_bed(genome_cds, cds_file)
  }

}
# From http://www.sthda.com/english/wiki/ggplot2-error-bars-quick-start-guide-r-software-and-data-visualization
#+++++++++++++++++++++++++
# Function to calculate the mean and the standard deviation
# for each group
#+++++++++++++++++++++++++
# data : a data frame
# varname : the name of a column containing the variable
#to be summariezed
# groupnames : vector of column names to be used as
# grouping variables
data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
}


#' Save GRanges as BED file
#' @param gr GRanges object
#' @param filename output file name
save_bed <- function(gr, filename, ignore_strand = TRUE) {
    export(
        gr,
        con = filename,
        format = "bed",
        ignore.strand = ignore_strand
    )
}
