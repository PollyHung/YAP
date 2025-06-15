library(Seurat)
library(dplyr)
library(readxl)
library(magrittr)
library(ggplot2)


## Metadata 
meta <- read.csv("metadata/GSE211956/metadata.csv", row.names = 1) %>% 
  dplyr::select("cells", "sample", "percent.mt", "mitoRatio", "log10GenesPerUMI", 
                "S.Score", "G2M.Score", "Phase", "YAP1", "YAP1_S94A", "YAP1_5SA", 
                "stiffness", "not_stiffness", "fibroblast_stiffness")
meta_list <- split(meta, as.factor(meta$sample))


# ## YAP1 signatures 
# YAP1 <- readxl::read_xlsx("YAP1.xlsx", col_names = F)
# YAP <- YAP1[1, 3:202] %>% unlist %>% unname
# YAP_S94A <- YAP1[2, 3:202] %>% unlist %>% unname
# YAP_5SA <- YAP1[3, 3:202] %>% unlist %>% unname
# 
# ## Fibriboasts signatures 
# fibroblasts <- c("COL11A1", "COMP", "FN1", "VCAN", "CTSB", "COL1A1")
# not_fibro <- c("LAMB1", "LAMC1", "HSPG2", "CTSG", "COL15A1", "LAMA4", "VWF", "COL6A6", "ABI3BP", "TNXB")
# t_cells <- c("CD4", "CD3E", "CD3D", "CD8A", "CD8B", "GNLY", "GZMB", "NCR1")
# myeloid <- c("CD14", "ITGAX", "HLA-DRA", "CD68", "CD163", "CSF1R", "ITGAM", "CD1C", "FCGR3A", "ADGRE1")
# b_cells <- c("CD19", "MS4A1", "BLNK", "BANK1", "IGHD", "IGHM", "IGKC")
# endothelial <- c("PECAM1", "VWF", "ICAM1", "FLT1", "KDR")

## for each sample.. 
samples <- list.dirs("data/spatial/GSE211956/", recursive = F, full.names = F)

for (sample in samples) {
  
  # read in the pre-processed spatial object 
  path <- paste0("data/spatial/GSE211956/", sample, "/output/seurat.rds")
  plot <- paste0("data/spatial/GSE211956/", sample, "/plots/")
  seurat <- readRDS(path)
  
  # # add module score 
  # seurat <- AddModuleScore(seurat, features = list(YAP), name = "YAP1")
  # seurat <- AddModuleScore(seurat, features = list(YAP_S94A), name = "YAP1_S94A")
  # seurat <- AddModuleScore(seurat, features = list(YAP_5SA), name = "YAP1_5SA")
  # seurat <- AddModuleScore(seurat, features = list(fibroblasts), name = "fibroblasts")
  # seurat <- AddModuleScore(seurat, features = list(not_fibro), name = "not_fibroblasts")
  # seurat <- AddModuleScore(seurat, features = list(t_cells), name = "T_cells")
  # seurat <- AddModuleScore(seurat, features = list(myeloid), name = "myeloid")
  # seurat <- AddModuleScore(seurat, features = list(b_cells), name = "B_cells")
  # seurat <- AddModuleScore(seurat, features = list(endothelial), name = "endothelial")
  # 
  # # edit names 
  # seurat@meta.data$YAP1 <- seurat@meta.data$YAP11; seurat@meta.data$YAP11 <- NULL
  # seurat@meta.data$YAP1_S94A <- seurat@meta.data$YAP1_S94A1; seurat@meta.data$YAP1_S94A1 <- NULL
  # seurat@meta.data$YAP1_5SA <- seurat@meta.data$YAP1_5SA1; seurat@meta.data$YAP1_5SA1 <- NULL
  # seurat@meta.data$fibroblast <- seurat@meta.data$fibroblasts1; seurat@meta.data$fibroblasts1 <- NULL
  # seurat@meta.data$not_fibroblast <- seurat@meta.data$not_fibroblasts1; seurat@meta.data$not_fibroblasts1 <- NULL
  # seurat@meta.data$T_cells <- seurat@meta.data$T_cells1; seurat@meta.data$T_cells1 <- NULL
  # seurat@meta.data$myeloid <- seurat@meta.data$myeloid1; seurat@meta.data$myeloid1 <- NULL
  # seurat@meta.data$B_cells <- seurat@meta.data$B_cells1; seurat@meta.data$B_cells1 <- NULL
  # seurat@meta.data$endothelial <- seurat@meta.data$endothelial1; seurat@meta.data$endothelial1 <- NULL
  
  ## Add the new metadata 
  meta_old <- seurat@meta.data %>% 
    dplyr::select("orig.ident", "nCount_Spatial", "nFeature_Spatial", "nCount_SCT", "nFeature_SCT")
  meta_new <- meta_list[[sample]]
  rownames(meta_new) <- gsub(sample, "", rownames(meta_new))
  rownames(meta_new) <- gsub("_", "", rownames(meta_new))
  meta_new$cells <- rownames(meta_new)
  meta_update <- merge(meta_old, meta_new, by = "row.names") %>% 
    tibble::column_to_rownames(var = "Row.names")
  meta_update <- meta_update[rownames(seurat@meta.data), ]
  seurat@meta.data <- meta_update
  
  ## Plot 
  p <- SpatialFeaturePlot(seurat, features = c("YAP1", "YAP1_S94A", "YAP1_5SA"), pt.size.factor = 3, image.alpha = 0) &  
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, limits = c(-0.2, 0.20))
  ggsave(file.path(plot, "YAP.png"), p, width = 7, height = 4, units = "in", dpi = 600)
  
  p2 <- SpatialFeaturePlot(seurat, features = c("stiffness"), pt.size.factor = 2, image.alpha = 0) &  
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, limits = c(-1, 2))
  ggsave(file.path(plot, "stiffness.png"), p2, width = 7, height = 4, units = "in", dpi = 600) 
  # save 
  # saveRDS(object = seurat, file = path)
}













# ## Clustering 
# p = DimPlot(spatial, group.by = "RNA_snn_res.0.5")
# 
# 
# ## Add YAP1 module scores 
# YAP1 = readxl::read_xlsx("data/YAP1.xlsx", col_names = F)
# YAP = YAP1[1, 3:202] %>% unlist %>% unname
# YAP_S94A = YAP1[2, 3:202] %>% unlist %>% unname
# YAP_5SA = YAP1[3, 3:202] %>% unlist %>% unname
# 
# spatial = AddModuleScore(spatial, features = list(YAP), name = "YAP1")
# spatial = AddModuleScore(spatial, features = list(YAP_S94A), name = "YAP1_S94A")
# spatial = AddModuleScore(spatial, features = list(YAP_5SA), name = "YAP1_5SA")
# 
# spatial@meta.data$YAP1 = spatial@meta.data$YAP11; spatial@meta.data$YAP11 = NULL
# spatial@meta.data$YAP1_S94A = spatial@meta.data$YAP1_S94A1; spatial@meta.data$YAP1_S94A1 = NULL
# spatial@meta.data$YAP1_5SA = spatial@meta.data$YAP1_5SA1; spatial@meta.data$YAP1_5SA1 = NULL
# 
# meta = spatial@meta.data
# # write.csv(meta, "metadata/GSE211956/metadata.csv", quote = F, row.names = T, col.names = T)
# 
# ## Stiffness 
# stiffness = c("COL11A1", "COMP", "FN1", "VCAN", "CTSB")
# down_reg = c("LAMB1", "LAMC1", "HSPG2", "CTSG", "COL15A1", "LAMA4", "VWF", "COL6A6", "ABI3BP", "TNXB")
# spatial = AddModuleScore(spatial, features = list(stiffness), name = "stiffness")
# spatial = AddModuleScore(spatial, features = list(down_reg), name = "not_stiffness")
# spatial@meta.data$stiffness = spatial@meta.data$stiffness1; spatial@meta.data$stiffness1 = NULL
# spatial@meta.data$not_stiffness = spatial@meta.data$not_stiffness1; spatial@meta.data$not_stiffness1 = NULL
# spatial$fibroblast_stiffness = spatial$stiffness - spatial$not_stiffness
# 
# meta = spatial@meta.data
# 
# write.csv()
# 
# ## Auto Labelling 
# ref = fetchReference("hpca", "2024-02-26")
# sce = as.SingleCellExperiment(spatial)
# singleR_results = SingleR(test = sce, ref = ref, labels = ref$label.main)
# spatial[['SingleR.labels']] = singleR_results$labels
# 























