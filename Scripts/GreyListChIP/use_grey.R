library(GreyListChIP)
library(BiocParallel)
library(Rsamtools)
library(GenomicAlignments)
library(SummarizedExperiment)

#https://github.com/ccagc/QDNAseq/issues/100

# Parse
# Function to check if an option is present in the command line arguments
get_option_presence <- function(option_name) {
  any(grepl(paste0("--", option_name), commandArgs()))
}

# Check if the "--new_func" option is present
new_func <- get_option_presence("new_func")

# Get the values of "--bed" and "--bam" options
bed_index <- grep("--bed", commandArgs())
bam_index <- grep("--bam", commandArgs())
kary_index <- grep("--kary", commandArgs())
n_cores_index <- grep("--n_cores", commandArgs())
seed_index <- grep("--seed", commandArgs())
gap_index <- grep("--max_gap", commandArgs())
yield_size_index <- grep("--yield_size", commandArgs())

# Check if "--bed" and "--bam" options are provided
if (length(bed_index) == 0 || length(bam_index) == 0) {
  stop("Error: Both --bed and --bam options must be provided.")
}

# Extract values of "--bed" and "--bam" options
bed <- commandArgs()[bed_index + 1]
bam <- commandArgs()[bam_index + 1]
kary <- commandArgs()[kary_index + 1]
n_cores <- commandArgs()[n_cores_index + 1]
seed <- commandArgs()[seed_index + 1]
max_gap <- as.numeric(commandArgs()[gap_index + 1])
yield_size <- commandArgs()[yield_size_index + 1]

set.seed(seed)

bam_path <- file.path(bam)
bed_path <- file.path(bed)

kary_path <- file.path(
  here::here(),
  kary
)

if (new_func) {
  countReads <- function(obj, bamFile) {
    fd <- BamFile(bamFile, yieldSize = yield_size)
    counts <- summarizeOverlaps(
      obj@tiles,
      fd,
      inter.feature = FALSE,
      ignore.strand = TRUE,
      mode = "IntersectionStrict"
    )
    obj@counts <- assays(counts)[[1]][, 1]
    obj@files <- c(obj@files, bamFile)
    return(obj)
  }
}


gl <- new("GreyList", karyoFile = kary_path)
gl <- countReads(gl, bam_path)
gl <- calcThreshold(gl, cores = n_cores)
gl <- makeGreyList(gl, maxGap = max_gap)

export(gl, con = bed_path)