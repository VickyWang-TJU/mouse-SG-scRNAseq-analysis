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
