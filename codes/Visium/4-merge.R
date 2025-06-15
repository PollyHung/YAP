library(Seurat)
library(dplyr)
library(magrittr)
library(tidyr)
library(ggplot2)
library(stringr)
library(scales)

setwd("~/Desktop/YAP1/")

## Load in Data and Cleave other stuff -----------------------------------------
visium <- readRDS("data/spatial/collapsed/spatial.rds")
xenium <- readRDS("~/Desktop/Collection/spatial-seq/10xXenium/Seurat/ST_xenium1_so.rds")

genes <- intersect(rownames(visium), rownames(xenium))

counts <- visium@assays$RNA$counts[genes, ]
visium <- CreateSeuratObject(counts, meta.data = visium@meta.data)
visium <- JoinLayers(visium)

counts <- xenium@assays$RNA$counts[genes, ]
xenium <- CreateSeuratObject(counts, meta.data = xenium@meta.data)
xenium <- JoinLayers(xenium)

rm(counts)


## Prepare Metadata ------------------------------------------------------------
xenium@meta.data$cells <- rownames(xenium@meta.data)
xenium@meta.data$sample <- xenium@meta.data$samples
visium_meta <- read.csv("data/spatial/collapsed/labelled_metadata.csv", row.names = 1)
visium_meta$cell.types <- visium_meta$cell_type
visium_meta$nCount_RNA <- visium_meta$nCount_Spatial
visium_meta$nFeature_RNA <- visium_meta$nFeature_Spatial
visium_meta$YAP1 <- NULL
visium_meta$patients <- visium_meta$sample

meta_cols <- intersect(colnames(xenium@meta.data), colnames(visium_meta))
visium_meta <- visium_meta[, meta_cols]
xenium_meta <- xenium@meta.data[, meta_cols]

metadata <- rbind(visium_meta, xenium_meta)

## Merge Two Objects ---------------------------------------------------
merged <- merge(visium, xenium, project = "merged")
cells_id <- intersect(rownames(metadata), colnames(merged))
merged <- subset(merged, cells = cells_id)
metadata <- metadata[colnames(merged), ]
merged@meta.data <- metadata
merged@meta.data$cell.types <- gsub("CARD.", "", merged@meta.data$cell.types)
merged@meta.data$cell.types <- ifelse(merged@meta.data$cell.types == "T.cell", "TNK.cell", merged@meta.data$cell.types)
merged@meta.data$cell.types <- ifelse(merged@meta.data$cell.types == "Myeloid.cell", "Monocyte", merged@meta.data$cell.types)
merged@meta.data$cell.types <- ifelse(merged@meta.data$cell.types == "Malignant", "Ovarian.cancer.cell", merged@meta.data$cell.types)
merged@meta.data$cell.types <- ifelse(merged@meta.data$cell.types == "Endothelial.cell", "Endothelial", merged@meta.data$cell.types)

merged <- JoinLayers(merged)

saveRDS(merged, "data/merged_visium_xenium.rds")


## Preprocessing ---------------------------------------------------------------
merged[["percent.mt"]] <- PercentageFeatureSet(merged, pattern = "^MT-")
merged <- NormalizeData(merged, normalization.method = "LogNormalize", scale.factor = 10000)
merged <- FindVariableFeatures(merged, selection.method = "vst", nfeatures = 3000)
regev_lab_cell_cycle_genes <- read.delim("data/regev_lab_cell_cycle_genes.txt", header=F) %>% unlist
s.genes <- regev_lab_cell_cycle_genes[1:43]
g2m.genes <- regev_lab_cell_cycle_genes[44:97]
merged <- CellCycleScoring(object = merged, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
merged <- ScaleData(merged, vars.to.regress = c("Phase", "orig.ident"), features = VariableFeatures(merged))
merged <- RunPCA(merged, features = VariableFeatures(merged))
merged <- harmony::RunHarmony(merged, 
                              group.by.vars = "sample", 
                              reduction = "pca", 
                              assay.use = "RNA", 
                              reduction.save = "harmony", 
                              lambda = 1)
merged <- FindNeighbors(merged, dims = 1:10, reduction = "harmony")
merged <- FindClusters(merged, resolution = 0.5, verbose = TRUE)
merged <- RunUMAP(merged, dims = 1:10, reduction="harmony")

DimPlot(merged) + NoLegend()
saveRDS(merged, "data/merged_visium_xenium.rds")


## Add Fibroblast Stiffness Scores --------------------------------------
fibro <- as.data.frame.matrix(table(merged$sample, merged$cell.types))
fibro$All.Cells <- rowSums(fibro)
fibro$Fibroblast.Pct <- 100*(fibro$Fibroblast / fibro$All.Cells)
cutoff <- median(fibro$Fibroblast.Pct)
fibro$Stiff <- ifelse(fibro$Fibroblast.Pct > cutoff, "high_stiff", "low_stiff")
merged$Stiffness <- ifelse(merged$sample %in% rownames(fibro)[fibro$Stiff == "high_stiff"], "high_stiff", "low_stiff")


## Add YAP1 scores -----------------------------------------------------
# load in the gene features 
YAP1 <- readxl::read_xlsx("~/Desktop/YAP1/data/YAP1.xlsx", col_names = F)
YAP <- YAP1[1, 3:202] %>% unlist %>% unname %>% intersect(rownames(merged))
YAP_S94A <- YAP1[2, 3:202] %>% unlist %>% unname %>% intersect(rownames(merged))
YAP_5SA <- YAP1[3, 3:202] %>% unlist %>% unname %>% intersect(rownames(merged))
# add original module scores 
merged <- AddModuleScore(merged, features = list(YAP), name = "YAP1", ctrl = 40)
merged <- AddModuleScore(merged, features = list(YAP_5SA), name = "YAP1_5SA", ctrl = 40)
merged <- AddModuleScore(merged, features = list(YAP_S94A), name = "YAP1_S94A", ctrl = 40)
# rename the columns 
merged$YAP1 <- merged$YAP11; merged$YAP11 <- NULL
merged$YAP1_5SA <- merged$YAP1_5SA1; merged$YAP1_5SA1 <- NULL
merged$YAP1_S94A <- merged$YAP1_S94A1; merged$YAP1_S94A1 <- NULL
# scale the scores 
merged$YAP1.norm <- rescale(merged$YAP1, to = c(0, 1))
merged$YAP1_5SA.norm <- rescale(merged$YAP1_5SA, to = c(0, 1))
merged$YAP1_S94A.norm <- rescale(merged$YAP1_S94A, to = c(0, 1))


## Add Cancer Heterogeneity ----------------------------------------------------
meta <- merged@meta.data
meta$cells <- rownames(meta)
df <- data.frame(row.names = unique(meta$sample), 
                 YAP1.hetero = rep(NA, length(unique(meta$sample))), 
                 YAP1_5SA.hetero = rep(NA, length(unique(meta$sample))), 
                 YAP1_S94A.hetero = rep(NA, length(unique(meta$sample))))
for(i in unique(meta$sample)){
  meta_sub <- meta %>% dplyr::filter(sample == i) %>% dplyr::filter(cell.types == "Ovarian.cancer.cell")
  df[i, "YAP1.hetero"] <- DescTools::Gini(meta_sub$YAP1.norm)
  df[i, "YAP1_5SA.hetero"] <- DescTools::Gini(meta_sub$YAP1.norm)
  df[i, "YAP1_S94A.hetero"] <- DescTools::Gini(meta_sub$YAP1.norm)
}
meta <- merge(meta, df, by.x = "sample", by.y = "row.names") %>% 
  tibble::column_to_rownames(var = "cells")
merged@meta.data <- meta[colnames(merged), ]

saveRDS(merged, "data/merged_visium_xenium.rds")


## Add back the metadata to each sub sample ------------------------------------
directories <- list.dirs("~/Desktop/YAP1/data/spatial/GSE211956", recursive = F)
meta <- merged@meta.data
for(i in directories){
  
  setwd(i)
  
  spatial <- readRDS("output/seurat.rds")
  spatial@meta.data$YAP1 <- NULL
  spatial@meta.data$YAP1_S94A <- NULL
  spatial@meta.data$YAP1_5SA <- NULL
  sample_id <- gsub("/Users/polly_hung/Desktop/YAP1/data/spatial/GSE211956/", "", i)
  sub_meta <- meta %>% dplyr::filter(sample == sample_id) %>% 
    dplyr::select("cell.types", "Stiffness", "YAP1", "YAP1_5SA", "YAP1_S94A", "YAP1.norm", "YAP1_5SA.norm", "YAP1_S94A.norm", "YAP1.hetero")
  
  rownames(sub_meta) <- gsub(sample_id, "", rownames(sub_meta))
  rownames(sub_meta) <- gsub("_", "", rownames(sub_meta))
  sub_meta2 <- merge(spatial@meta.data, sub_meta, by = "row.names") %>% 
    tibble::column_to_rownames(var = "Row.names")
  sub_meta2 <- sub_meta2[rownames(spatial@meta.data), ]
  
  spatial@meta.data <- sub_meta2
  SpatialFeaturePlot(spatial, features = "YAP1_5SA.norm", pt.size.factor = 3, image.alpha = 0) 
  
  saveRDS(spatial, "output/labelled2.rds")

}

