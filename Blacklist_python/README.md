# README

## Overview

This repository contains two Python scripts for generating excludable regions (blacklists) from BAM files and mappability data, inspired by [Boyle Lab's blacklist.cpp](https://github.com/Boyle-Lab/Blacklist/blob/master/blacklist.cpp). These tools enable robust computation of genomic regions unsuitable for analysis due to high noise or poor mappability.

### Scripts

- **[blacklist_literal.py](blacklist_literal.py)**:  
  A direct translation of the original `blacklist.cpp` into Python, line by line, designed for compatibility with its inputs and outputs.
  - **Inputs**:
    - `BAM file(s)`: Aligned reads for analysis.
    - `mappability file(s)`: Binary uint8 files from [umap](https://github.com/hoffmangroup/umap) or [Hoffman Lab](https://bismap.hoffmanlab.org/).
  - **Outputs**:
    - `BED file`: A terminal output containing excludable regions in BED format with annotations.

- **[blacklist2.py](blacklist2.py)**:  
  An enhanced version of `blacklist.cpp`, introducing new features and configurability.  
  - **Inputs**:
    - `BAM file(s)` and `mappability file(s)` as described above.
  - **Outputs**:
    - `BED file`: A detailed BED file with excludable regions and annotations, output to a file or the terminal.
  - **Parameters**:
    (Refer to the "Parameters" section below for detailed options.)

---

## Features

### `blacklist_literal.py`
- Faithful implementation of `blacklist.cpp`.
- Identical input and output requirements for seamless integration.

### `blacklist2.py`
- Highly configurable with additional parameters for:
  - Custom bin sizes and overlaps.
  - Fine-tuned region-specific analysis.
  - Output flexibility (file or terminal).
  - Advanced threshold control for weak and strong percentiles.
  - Compatibility with original and improved algorithms.

---

## Usage

### Example: `blacklist_literal.py`
```bash
python blacklist_literal.py \
    --bam ./input/sample.bam \
    --mappability ./mappability/sample.mappability
```

### Example: `blacklist2.py`
```bash
python blacklist2.py \
    -b ./input/sample1.bam,./input/sample2.bam \
    -m ./mappability/ \
    -r chr1,chr2 \
    -i 500 \
    -p 50 \
    -g 100 \
    -u 50 \
    -w 0.98 \
    -s 0.999 \
    -o ./output/blacklist.bed \
    -v
```

### Command-line Options for `blacklist2.py`:

| Parameter          | Description                                                                                          | Default             |
|--------------------|------------------------------------------------------------------------------------------------------|---------------------|
| `-b, --bams`       | Comma-separated list of BAM files or directories to process.                                         | `./input`           |
| `-m, --mappability`| Directory containing mappability files.                                                             | `./mappability`     |
| `-r, --regions`    | Restrict analysis to specific regions (comma-separated, e.g., `chr1,chr2`).                          | `all`               |
| `-i, --bin`        | Size of bins to use for analysis.                                                                   | `1000`              |
| `-p, --overlap`    | Overlap between bins.                                                                               | `100`               |
| `-g, --bridge`     | Size of bridges between regions.                                                                    | `200`               |
| `-u, --unique`     | Minimum length for a read to be considered unique.                                                  | `36`                |
| `-w, --weak`       | Percentile threshold for weak regions.                                                              | `0.99`              |
| `-s, --strong`     | Percentile threshold for strong regions.                                                            | `0.999`             |
| `-o, --output`     | Path to save the output `.bed` file.                                                                | `./b2output.bed`    |
| `-v, --view`       | Print results to terminal instead of saving to a file.                                              | Disabled            |
| `-n, --noMerge`    | Prevent merging of overlapping regions.                                                             | Disabled            |

---

## Installation and Prerequisites

The scripts require Python (>=3.8) and the following Python packages:
- `numpy`
- `pandas`
- `pysam`

You can install the required dependencies using the provided Conda environment file.

### Conda Environment Setup

1. Create the environment:
   ```bash
   conda env create -f blacklist-env.yml
   ```

2. Activate the environment:
   ```bash
   conda activate blacklist-env
   ```
   