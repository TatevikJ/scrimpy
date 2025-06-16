# SCRIMPy: Single Cell Replication Inference from Multiome data using Python

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A523.04.0-23aa62.svg)](https://www.nextflow.io/)

## Table of Contents

- [Overview](#overview)
- [Quick Start](#quick-start)
  - [Prerequisites](#prerequisites)
  - [Installation](#installation)
  - [Basic Usage](#basic-usage)
- [Input Requirements](#input-requirements)
  - [Sample Sheet Format](#sample-sheet-format)
  - [Required Parameters](#required-parameters)
  - [Optional Parameters](#optional-parameters)
- [Output Structure](#output-structure)
  - [Key Output Files](#key-output-files)


## Overview

**SCRIMPy** is a Nextflow pipeline that processes single-cell ATAC-seq fragment data to infer DNA replication states of individual cells. The main steps are:

1. Filter ATAC fragments using provided list of cell barcodes and retain only autosomal fragments for each cell.
2. Infer replication timing from fragment coverage if replication timing data is not provided.
3. Assign replication timing rank to each fragment.
4. Identify overlapping fragments.
5. Use obtained replication and overlap information of fragments to compute per-cell metrics file.

## Quick Start

### Prerequisites

- [Nextflow](https://www.nextflow.io/docs/latest/getstarted.html#installation) (>=23.04.2)
- Java 11 or later

### Installation

```bash
git clone https://github.com/TatevikJ/scrimpy.git
cd scrimpy
```

### Basic Usage

```bash
nextflow run main.nf -c nextflow.config 
```

## Input Requirements

### Required Parameters

| Parameter | Description |
|-----------|-------------|
| `input_paths_csv` | Path to sample sheet (CSV format) |
| `outdir` | Output directory |

### Optional Parameters

| Parameter           | Description                                                                 | Default |
|---------------------|-----------------------------------------------------------------------------|---------|
| `rt_tsv`            | Path to replication timing data                                             | `None` (required if `rt_chrom_sizes` and `rt_fragments` are not provided)|
| `rt_chrom_sizes`    | Path to `chrom.sizes` file for generating genomic bins for RT inference    | `None` (required if `rt_tsv` is not provided) |
| `rt_fragments`      | Path to a fragment file (`atac_fragments.tsv.gz`) used for RT inference     | `None` (required if `rt_tsv` is not provided) |
| `mem_per_cpu`       | Memory allocated per CPU                                                    | `40.GB` |
| `cpus`              | Number of CPUs to use per task                                              | `6` |

### Sample Sheet Format

The pipeline requires a CSV sample sheet with the following columns:

```csv
sample,fragments,barcodes
sample1,/path/to/sample1/atac_fragments.tsv.gz,/path/to/sample1/barcode_list.tsv
sample2,/path/to/sample2/atac_fragments.tsv.gz,/path/to/sample2/barcode_list.tsv
```

## Output Structure

```
outdir/
├── reports/                 # Pipeline execution reports
├── sample1/                 # Sample 1 analysis results
│   └── metrics_per_cell.csv # Per-cell metrics for sample 1
└── sample2/                 # Sample 2 analysis results
    └── metrics_per_cell.csv # Per-cell metrics for sample 2
```

