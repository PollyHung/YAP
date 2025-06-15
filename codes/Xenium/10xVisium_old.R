library(Seurat)
library(dplyr)
library(magrittr)
library(tidyr)
library(ggplot2)


## Load in Data and Cleave other stuff ---------------------------------
obj <- readRDS("data/collapsed_spatial/collapsed_spatial.rds")
test <- readRDS("~/Desktop/Collection/spatial-seq/10xXenium/Seurat/ST_Test1_so.rds")

genes <- intersect(rownames(test), rownames(obj))

counts <- obj@assays$RNA$counts[genes, ]
obj <- CreateSeuratObject(counts, meta.data = obj@meta.data)
obj <- JoinLayers(obj)

counts <- test@assays$RNA$counts[genes, ]
test <- CreateSeuratObject(counts, meta.data = test@meta.data)
test <- JoinLayers(test)

rm(counts)


## Prepare Metadata 
test@meta.data$cells <- rownames(test@meta.data)
test@meta.data$sample <- test@meta.data$samples
obj@meta.data$patients <- obj@meta.data$sample


## Merge Two Objects ---------------------------------------------------
merged <- merge(obj, test, project = "merged")
obj <- obj@meta.data
test <- test@meta.data

merged@meta.data <- merged@meta.data %>% dplyr::select(sample, patients, cells, nCount_RNA, nFeature_RNA)
merged <- JoinLayers(merged)
merged$orig.ident <- ifelse(grepl("GSM", merged$sample), "set1", "set2")


## Preprocessing 
merged[["percent.mt"]] <- PercentageFeatureSet(merged, pattern = "^MT-")
merged <- NormalizeData(merged, normalization.method = "LogNormalize", scale.factor = 10000)
merged <- FindVariableFeatures(merged, selection.method = "vst", nfeatures = 5000)
regev_lab_cell_cycle_genes <- read.delim("data/regev_lab_cell_cycle_genes.txt", header=F) %>% unlist
s.genes <- regev_lab_cell_cycle_genes[1:43]
g2m.genes <- regev_lab_cell_cycle_genes[44:97]
merged <- CellCycleScoring(object = merged, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
merged <- ScaleData(merged, vars.to.regress = c("Phase", "orig.ident"), features = VariableFeatures(object = merged))
merged <- RunPCA(merged, features = VariableFeatures(object = merged))
DimPlot(merged, group.by = "Phase", reduction = "pca") ## scale for Phase 
DimPlot(merged, split.by = "orig.ident", group.by = "orig.ident", reduction = "pca") ## scale for datasets 

merged <- ScaleData(merged, vars.to.regress = c("orig.ident", "sample"))
merged <- RunPCA(merged, features = VariableFeatures(object = merged))

merged <- harmony::RunHarmony(merged, 
                              group.by.vars = c("sample", "orig.ident"), 
                              reduction = "pca", 
                              assay.use = "RNA", 
                              reduction.save = "harmony", 
                              lambda = 1)
merged <- FindNeighbors(merged, dims = 1:20, reduction = "harmony")
merged <- FindClusters(object = merged, resolution = c(0.1, 0.3, 0.5, 0.7), verbose = TRUE)
merged <- RunUMAP(merged, dims = 1:10, reduction="harmony")

DimPlot(merged) + NoLegend()


## Add Back the Cell Labels 
test <- test %>% dplyr::select(cells, cell.types)
obj$cell.types <- "unknown"
obj <- obj %>% dplyr::select(cells, cell.types)
meta_combined <- rbind(test, obj)
meta_combined <- meta_combined[colnames(merged), ]
meta_combined <- merge(merged@meta.data, meta_combined, by = "cells") %>% 
  tibble::column_to_rownames(var = "cells")
meta_combined$cells <- rownames(meta_combined)
meta_combined <- meta_combined[colnames(merged), ]
merged@meta.data <- meta_combined


## Color by cell types 
custom_colors <- c("grey", "red", "green", "blue", "purple", "pink", "orange")
names(custom_colors) <- unique(merged$cell.types)
DimPlot(merged, group.by = "cell.types", cols = custom_colors)


## Add Fibroblast Stiffness Scores --------------------------------------
fibro <- as.data.frame.matrix(table(obj$samples, obj$cell.types))
fibro$All.Cells <- rowSums(fibro)
fibro$Fibroblast.Pct <- 100*(fibro$Fibroblast / fibro$All.Cells)
fibro$Stiff <- ifelse(fibro$Fibroblast.Pct > 50, "high_stiff", "low_stiff")
obj$Stiffness <- ifelse(obj$samples %in% rownames(fibro)[fibro$Stiff == "high_stiff"], "high_stiff", "low_stiff")


## Add YAP1 scores -----------------------------------------------------
# load in the gene features 
YAP1 <- readxl::read_xlsx("~/Desktop/YAP1/data/YAP1.xlsx", col_names = F)
YAP <- YAP1[1, 3:202] %>% unlist %>% unname %>% intersect(rownames(obj))
YAP_S94A <- YAP1[2, 3:202] %>% unlist %>% unname %>% intersect(rownames(obj))
YAP_5SA <- YAP1[3, 3:202] %>% unlist %>% unname %>% intersect(rownames(obj))
# add original module scores 
obj <- AddModuleScore(obj, features = list(YAP), name = "YAP1", slot = "scale.data")
obj <- AddModuleScore(obj, features = list(YAP_5SA), name = "YAP1_5SA", slot = "scale.data")
obj <- AddModuleScore(obj, features = list(YAP_S94A), name = "YAP1_S94A", slot = "scale.data")
# rename the columns 
obj$YAP1 <- obj$YAP11; obj$YAP11 <- NULL
obj$YAP1_5SA <- obj$YAP1_5SA1; obj$YAP1_5SA1 <- NULL
obj$YAP1_S94A <- obj$YAP1_S94A1; obj$YAP1_S94A1 <- NULL
# scale the scores 
obj$YAP1.norm <- rescale(obj$YAP1, to = c(0, 1))
obj$YAP1_5SA.norm <- rescale(obj$YAP1_5SA, to = c(0, 1))
obj$YAP1_S94A.norm <- rescale(obj$YAP1_S94A, to = c(0, 1))



## Add YAP1 colours ----------------------------------------------------
# Make Malignant cell Colours by Bin
add_colour <- function(obj, var, bin_name, color_name){
  # create color 
  cancer <- obj@meta.data %>% dplyr::filter(cell.types == "Malignant")
  cancer[[bin_name]] <- cut(cancer[[var]], 
                            breaks = quantile(cancer[[var]], probs = seq(0, 1, length.out = 11), na.rm = TRUE), 
                            include.lowest = TRUE)
  bin_colors <- colorRampPalette(c("#344CB7", "white", "#E72929"))(10)
  names(bin_colors) <- levels(unique(cancer[[bin_name]])) ## ordered names 
  cancer[[color_name]] <- bin_colors[cancer[[bin_name]]]
  cancer <- cancer[ , c(bin_name, color_name)]
  # add color 
  temp <- merge(obj@meta.data, cancer, by = "row.names", all.x = TRUE) %>% 
    tibble::column_to_rownames(var = "Row.names")
  obj@meta.data <- temp[colnames(obj), ]
  return(obj)
}
obj <- add_colour(obj = obj, var = "YAP1", bin_name = "YAP1.bins", color_name = "YAP1.color")
obj <- add_colour(obj = obj, var = "YAP1_5SA", bin_name = "YAP1_5SA.bins", color_name = "YAP1_5SA.color")
obj <- add_colour(obj = obj, var = "YAP1_S94A", bin_name = "YAP1_S94A.bins", color_name = "YAP1_S94A.color")

# Add other Colors 
obj$YAP1.color <- ifelse(is.na(obj$YAP1.color), ifelse(obj$cell.types == "Fibroblast", "#73EC8B", "#3C3D37"), obj$YAP1.color)
obj$YAP1_5SA.color <- ifelse(is.na(obj$YAP1_5SA.color), ifelse(obj$cell.types == "Fibroblast", "#73EC8B", "#3C3D37"), obj$YAP1_5SA.color)
obj$YAP1_S94A.color <- ifelse(is.na(obj$YAP1_S94A.color), ifelse(obj$cell.types == "Fibroblast", "#73EC8B", "#3C3D37"), obj$YAP1_S94A.color)


## Add Cancer Heterogeneity ----------------------------------------------------
meta <- obj@meta.data
meta$cells <- rownames(meta)
df <- data.frame(row.names = unique(meta$samples), 
                 YAP1.hetero = rep(NA, length(unique(meta$samples))), 
                 YAP1_5SA.hetero = rep(NA, length(unique(meta$samples))), 
                 YAP1_S94A.hetero = rep(NA, length(unique(meta$samples))))
for(i in unique(meta$samples)){
  meta_sub <- meta %>% dplyr::filter(samples == i) %>% dplyr::filter(cell.types == "Malignant")
  df[i, "YAP1.hetero"] <- Gini(meta_sub$YAP1.norm)
  df[i, "YAP1_5SA.hetero"] <- Gini(meta_sub$YAP1.norm)
  df[i, "YAP1_S94A.hetero"] <- Gini(meta_sub$YAP1.norm)
}
meta <- merge(meta, df, by.x = "samples", by.y = "row.names") %>% 
  tibble::column_to_rownames(var = "cells")
obj@meta.data <- meta[colnames(obj), ]

