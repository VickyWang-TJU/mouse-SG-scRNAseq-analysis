# Single-cell RNA-seq analysis of mouse thoracolumbar sympathetic ganglia

This repository contains R scripts used for the analysis and visualization of 10x Genomics single-cell RNA-seq data from mouse thoracolumbar sympathetic ganglia (SG) collected at the T12/T13/L1 levels.

The analysis includes quality control, preprocessing, clustering, cell-type annotation, sympathetic neuron subclustering, marker gene visualization, and differential gene expression analysis.

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
└── output/
    ├── figures/
    └── tables/
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

The script expects 10x Genomics-formatted matrix files for each sample:

```text
barcodes.tsv.gz
features.tsv.gz
matrix.mtx.gz
```

The input data directory should be modified in the script before running:

```r
data_path <- "path/to/data/"
```

## Analysis workflow

The main analysis script performs the following steps:

1. Load 10x Genomics single-cell count matrices using `Read10X`.
2. Create Seurat objects and merge samples.
3. Perform quality control based on:
   - number of detected genes
   - total UMI counts
   - mitochondrial gene percentage
4. Normalize data and identify highly variable genes.
5. Perform PCA, UMAP dimensionality reduction and clustering.
6. Annotate major SG cell populations using canonical marker genes.
7. Visualize marker gene expression using UMAP feature plots, violin plots and dot plots.
8. Subset sympathetic neurons and perform subclustering.
9. Identify cluster-enriched genes and visualize selected marker genes.
10. Analyze expression of genes related to GAG metabolism, including `Gusb`.

## Cell-type annotation

Major cell populations were annotated based on canonical marker genes, including:

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
- UMAP plots
- marker gene dot plots
- feature plots
- sympathetic neuron subcluster plots
- differential expression results
- selected gene expression visualizations

## Notes

Before running the script, users should update local file paths and sample names according to their own directory structure.

This repository provides the analysis code used for data processing and figure generation. Raw sequencing files and processed count matrices are deposited separately in GEO.
