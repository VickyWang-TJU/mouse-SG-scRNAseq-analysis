# Single-cell RNA-seq analysis of mouse thoracolumbar sympathetic ganglia

This repository contains an R script used for the analysis and visualization of 10x Genomics single-cell RNA-seq data from mouse thoracolumbar sympathetic ganglia (SG) collected at the T12/T13/L1 levels.

The analysis includes quality control, preprocessing, clustering, manual cell-type annotation, sympathetic neuron subclustering, and marker gene visualization.

## Data

Processed single-cell RNA-seq count matrices will be deposited in GEO as:

- `barcodes.tsv.gz`
- `features.tsv.gz`
- `matrix.mtx.gz`

for each SG sample.

The GEO accession number will be added upon release.

## Repository structure

```text
.
├── README.md
├── scripts/
│   └── SG_scRNAseq_analysis.R
├── data/
│   ├── data1/
│   │   ├── barcodes.tsv.gz
│   │   ├── features.tsv.gz
│   │   └── matrix.mtx.gz
│   └── data2/
│       ├── barcodes.tsv.gz
│       ├── features.tsv.gz
│       └── matrix.mtx.gz
└── output/
```

## Requirements

The analysis was performed in R using the following packages:

- Seurat
- tidyverse
- ggplot2
- patchwork
- Matrix
- harmony
- SingleR
- celldex
- UCell
- stringr
- plyr

Package versions may vary depending on the R environment. Users are advised to install recent compatible versions of these packages.

## Input files

The script expects two 10x Genomics-formatted input directories:

```text
data/data1/
data/data2/
```

Each input directory should contain:

```text
barcodes.tsv.gz
features.tsv.gz
matrix.mtx.gz
```

The script automatically detects the project root when run from the `scripts/` directory and reads input files from:

```text
./data/data1
./data/data2
```

## Running the analysis

Run the main analysis script from the repository root:

```bash
Rscript scripts/SG_scRNAseq_analysis.R
```

All generated plots will be saved to:

```text
./output/
```

## Analysis workflow

The main analysis script performs the following steps:

1. Load 10x Genomics single-cell count matrices using `Read10X`.
2. Create Seurat objects for `data1` and `data2`.
3. Merge the two Seurat objects.
4. Calculate mitochondrial gene percentage.
5. Generate quality-control violin plots.
6. Filter cells based on detected genes, UMI counts, and mitochondrial gene percentage.
7. Normalize data and identify highly variable genes.
8. Scale data and perform PCA.
9. Generate a PCA elbow plot.
10. Perform neighbor finding, clustering, and UMAP dimensionality reduction.
11. Annotate major SG cell populations manually based on known marker genes.
12. Visualize annotated cell populations and marker gene expression.
13. Subset sympathetic neurons.
14. Perform quality control, normalization, clustering, and UMAP analysis for sympathetic neurons.
15. Generate selected visualization plots for sympathetic neuron subclusters and marker genes.

## Cell-type annotation

Major cell populations are annotated manually in the script, including:

- Sympathetic neurons
- Schwann cells
- Vascular cells
- Satellite glia
- Fibroblasts
- Immune cells

## Output

The script generates representative analysis outputs, including:

- quality-control violin plots
- PCA elbow plot
- UMAP clustering plots
- manual cell-type annotation plots
- marker gene dot plots
- `Gusb` violin and feature plots
- sympathetic neuron subcluster plots
- selected `Gusb` and `Th` expression plots

All output files are saved directly in:

```text
./output/
```

## Notes

This repository provides the analysis code used for data processing and figure generation.

Raw sequencing files and processed count matrices will be deposited in GEO.
