library(patchwork)  
library(harmony)
library(BiocManager)
library(Seurat)
library(tidyverse)
library(celldex)
library(Matrix)
library(SingleR)
library(ggplot2)
library(plyr)
library(UCell)
library(stringr)

get_script_dir <- function() {
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_arg) > 0) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg[1]))))
  }
  getwd()
}

script_dir <- get_script_dir()
project_root <- if (basename(script_dir) == "scripts") {
  normalizePath(file.path(script_dir, ".."))
} else {
  normalizePath(getwd())
}

output_dir <- file.path(project_root, "output")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

save_plot <- function(filename, ...) {
  ggsave(filename = file.path(output_dir, filename), ...)
}

data_path <- file.path(project_root, "data")
data1 <- Read10X(file.path(data_path, "data1"))
data2 <- Read10X(file.path(data_path, "data2"))
data1_seurat <- CreateSeuratObject(data1, project = "data1", 
                                  min.cells = 3, min.features = 200)
data2_seurat <- CreateSeuratObject(data2, project = "data2", 
                                  min.cells = 3, min.features = 200)
seurat <- merge(data1_seurat, y = data2_seurat, 
                add.cell.ids = c("data1", "data2"))
seurat[["percent.mt"]] <- PercentageFeatureSet(seurat, pattern = "^mt-")
p_qc <- VlnPlot(seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
                pt.size = 0.1, ncol = 3) +
  ggtitle("Quality Control Metrics")
print(p_qc)
save_plot("QC_metrics.pdf", plot = p_qc, width = 12, height = 5)
library(patchwork)

vln_plots <- VlnPlot(seurat, 
                     features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
                     pt.size = 0.1, 
                     group.by = "single_sample",
                     combine = FALSE)

for(i in 1:3) {
  vln_plots[[i]] <- vln_plots[[i]] + 
    theme(axis.text.x = element_blank(),
          axis.title.x = element_blank())
}

p_qc <- wrap_plots(vln_plots, nrow = 1) +
  plot_annotation(title = "Quality Control Metrics", 
                  theme = theme(plot.title = element_text(hjust = 0.5, size = 16)))

print(p_qc)
save_plot("QC_metrics.pdf", plot = p_qc, width = 12, height = 5)

seurat <- subset(seurat, 
                 subset = nFeature_RNA > 500 & nFeature_RNA < 6000 &
                   percent.mt < 10 &
                   nCount_RNA > 1000 & nCount_RNA < 30000)

seurat <- NormalizeData(seurat)
seurat <- FindVariableFeatures(seurat, nfeatures = 2000)
seurat <- ScaleData(seurat)
seurat <- RunPCA(seurat, npcs = 30)
p_elbow <- ElbowPlot(seurat, ndims = 30) + 
  geom_vline(xintercept = 20, linetype = "dashed", color = "red") +
  ggtitle("PCA Elbow Plot")
print(p_elbow)
save_plot("PCA_elbow_plot.pdf", plot = p_elbow, width = 5, height = 5)
seurat <- FindNeighbors(seurat, dims = 1:10)
seurat <- FindClusters(seurat, resolution = 0.2)
seurat <- RunUMAP(seurat, dims = 1:10)
p_cluster <- DimPlot(seurat, reduction = "umap", 
                     label = TRUE, label.size = 4, repel = TRUE) +
  ggtitle("UMAP Clustering")
print(p_cluster)
celltype_annotation <- list(
  Sympathetic_Neurons = c("1", "2", "8"),
  Schwann_Cells = c("0"),
  Vascular_Cells = c("3"),
  Satellite_Glia = c("7"),
  Fibroblasts = c("4", "6"),
  Immune_Cells = c("5")
)

seurat$celltype_manual <- as.character(Idents(seurat))
for(celltype in names(celltype_annotation)) {
  clusters <- celltype_annotation[[celltype]]
  seurat$celltype_manual[seurat$celltype_manual %in% clusters] <- celltype
}

Idents(seurat) <- factor(seurat$celltype_manual, 
                         levels = names(celltype_annotation))

p_annotation <- DimPlot(seurat, label = FALSE, repel = FALSE, 
                        label.size = 4, pt.size = 0.5) +
  ggtitle("Manual Cell Type Annotation")
print(p_annotation)
marker_genes <- c("Th","Sox10", "Pecam1", "Pdgfrb", "Acta2", "Slc1a3" , "Col1a1","Cd3d", "Cd68")

p <- DotPlot(seurat, features = marker_genes) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("Marker Genes Validation") +
  labs(x = "Cell Types", y = "Genes")

print(p)
save_plot("marker_validation_dotplot.pdf", plot = p, width = 8, height =5)
VlnPlot(seurat, features = "Gusb")
save_plot("Gusb_vlnplot.pdf", width = 8, height = 5)
FeaturePlot(seurat, features = "Gusb")
save_plot("Gusb_FeaturePlot.pdf", width =4, height = 5)
symp_cells <- WhichCells(seurat, idents = "Sympathetic_Neurons")
seurat_symp <- subset(seurat, cells = symp_cells)
seurat_symp <- subset(seurat_symp,
                      subset = nFeature_RNA > 400 & 
                        nFeature_RNA < 5000 & 
                        percent.mt < 15 & 
                        nCount_RNA > 600 & 
                        nCount_RNA < 30000)

print(paste("质控后细胞数:", ncol(seurat_symp)))
seurat_symp <- NormalizeData(seurat_symp)
seurat_symp <- FindVariableFeatures(seurat_symp)
seurat_symp <- ScaleData(seurat_symp)
seurat_symp <- RunPCA(seurat_symp)
seurat_symp <- FindNeighbors(seurat_symp, dims = 1:15)
seurat_symp <- FindClusters(seurat_symp, resolution = 0.2)
seurat_symp <- RunUMAP(seurat_symp, dims = 1:15,
                       n.neighbors = 10, min.dist = 0.4)
p_cluster <- DimPlot(seurat_symp, label = TRUE) + 
  ggtitle("Sympathetic Neurons Subclusters")
print(p_cluster)
save_plot("seurat_symp.pdf", p_cluster, width = 5, height = 5)
seurat_symp_clean <- seurat_symp
my_colors <- c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf")
p1 <- DimPlot(seurat_symp_clean, 
              cols = my_colors[1:length(unique(seurat_symp_clean$seurat_clusters))],
              pt.size = 1) + 
  ggtitle("Sympathetic Neurons")
print(p1)
save_plot("DimPlot_Sympathetic_Neurons.pdf", p1, width = 5, height = 5)
p2 <- FeaturePlot(seurat_symp_clean, features = 'Gusb', label = TRUE, pt.size = 1,
                  cols = c("lightgrey", "#d62728")) +
  ggtitle("Gusb Expression")
p2$data <- p2$data[order(p2$data[, "Gusb"], decreasing = FALSE), ]
print(p2)
save_plot("FeaturePlot_Gusb.pdf", plot = p2, width = 5, height = 5)
p3 <- DotPlot(seurat_symp_clean, features = 'Gusb',
              cols = c("lightgrey", "#d62728")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(p3)
save_plot("DotPlot_Gusb_red.pdf", p3, width = 3, height = 5)
p4 <- FeaturePlot(seurat_symp_clean, features = 'Th', label = TRUE, pt.size = 1,
                  cols = c("lightgrey", "#d62728")) +
  ggtitle("Th Expression")
print(p4)
save_plot("FeaturePlot_Th.pdf", plot = p4, width = 5, height = 5)
