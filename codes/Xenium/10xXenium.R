library(dplyr)
library(magrittr)
library(Seurat)
library(DropletUtils)
library(Matrix)
library(ggplot2)
library(ggpubr)
library(lpSolveAPI)
library(EBImage)
library(scales)
library(umap)
set.seed(123)

setwd("~/Desktop/Collection/spatial-seq/10xXenium/Seurat/")

## Processing ------------------------------------------------------------------
seurat <- readRDS("ST_Discovery_so.rds")

assay5 <- as(spatial[["RNA"]], Class = "Assay5")
seurat <- CreateSeuratObject(assay5, meta.data = spatial@meta.data)
seurat <- JoinLayers(seurat)
seurat$cell.types <- gsub("_LC", "", seurat$cell.types)

rm(spatial)
rm(assay5)



# Add Stiffness Scores 
stiffness <- c("COL11A1", "COMP", "FN1", "VCAN", "CTSB") %>% intersect(rownames(seurat))
seurat <- AddModuleScore(seurat, features = list(stiffness), name = "stiffness", ctrl = 40, slot = "scale.data")
seurat$stiffness <- seurat$stiffness1; seurat$stiffness1 <- NULL

# ## Use K-means clustering to determine which sample has high-stiffness Tumour 
# fibroblast <- seurat@meta.data %>%
#   dplyr::filter(cell.types == "Fibroblast") %>%
#   dplyr::select(samples, cell.types,cell.subtypes, stiffness)
# kmeans_result <- kmeans(fibroblast$stiffness, centers = 2)  # 2 clusters for high and low stiffness
# fibroblast$cluster <- as.factor(kmeans_result$cluster)
# k_means <- table(fibroblast$samples, fibroblast$cluster) %>% as.data.frame.matrix()
# colnames(k_means) <- c("low stiffness", "high stiffness")
# k_means$pct_high <- k_means$`high stiffness`/(k_means$`low stiffness` + k_means$`high stiffness`) * 100
# 
# # Remove samples with less than 10 fibroblast cells 
# k_means <- k_means %>% dplyr::filter(`low stiffness` + `high stiffness` > 10)
# 
# # Samples with more than 40% of Fibroblast Cells being classified as High-stiffness 
# # are considered as high-stiff samples 
# n <- 50
# seurat$sample_group <- ifelse(seurat$samples %in% rownames(k_means)[k_means$pct_high > n],
#                               "high_stiffness", "low_stiffness")

## Use Mean Stiffness By Sample, > 0 is highly stiffed, < 0 is not stiffed 
mean_stiffness <- seurat@meta.data %>%
  group_by(samples) %>%        
  summarize(mean_stiffness = mean(stiffness, na.rm = TRUE)) 

seurat$sample_group <- ifelse(seurat$samples %in% mean_stiffness$samples[mean_stiffness$mean_stiffness > 0],
                              "high_stiffness", "low_stiffness")
ggplot(seurat@meta.data, aes(x = sample_group, y = stiffness)) + geom_violin()


# YAP signatures 
YAP1 <- readxl::read_xlsx("~/Desktop/YAP1/data/YAP1.xlsx", col_names = F)
YAP <- YAP1[1, 3:202] %>% unlist %>% unname %>% intersect(rownames(seurat))
YAP_S94A <- YAP1[2, 3:202] %>% unlist %>% unname %>% intersect(rownames(seurat))
YAP_5SA <- YAP1[3, 3:202] %>% unlist %>% unname %>% intersect(rownames(seurat))

seurat <- AddModuleScore(seurat, features = list(YAP), name = "YAP", ctrl = 40, slot = "scale.data")
seurat <- AddModuleScore(seurat, features = list(YAP_5SA), name = "YAP_5SA", ctrl = 40, slot = "scale.data")
seurat <- AddModuleScore(seurat, features = list(YAP_S94A), name = "YAP_S94A", ctrl = 40, slot = "scale.data")

seurat$YAP1_norm <- rescale(seurat$YAP1, to = c(0, 1))
seurat$YAP_5SA1_norm <- rescale(seurat$YAP_5SA1, to = c(0, 1))
seurat$YAP_S94A1_norm <- rescale(seurat$YAP_S94A1, to = c(0, 1))


## Make Colour ------------------------------------------------------------------
p <- ggplot(meta, aes(x = x, y = y)) + 
  geom_point(aes(color = cell.types), size = 1, alpha = 1) +
  scale_color_manual(values = c("Fibroblast" = "#00FF9C", "Monocyte" = "black", "TNK.cell" = "black",
                                "B.cell" = "black", "Endothelial" = "black", "Mast.cell" = "black", "Malignant" = "white")) 
pd1 <- ggplot_build(p)$data[[1]]

meta <- meta %>% mutate(color = ifelse(cell.types == "Malignant", YAP1, NA))
p2 <- ggplot(meta, aes(x = x, y = y)) + geom_point(aes(color = color), size = 1) + 
  scale_color_gradient2(midpoint = 0, low = "#006BFF", mid = "#FFF100", high = "#E52020", na.value = "#3C3D37", 
                        limits = c(range(meta$YAP1)[[1]], range(meta$YAP1)[[2]])) 
pd2 <- ggplot_build(p2)$data[[1]]

meta <- meta %>% mutate(color = ifelse(cell.types == "Malignant", YAP_S94A1, NA))
p3 <- ggplot(meta, aes(x = x, y = y)) + geom_point(aes(color = color), size = 1) + 
  scale_color_gradient2(midpoint = 0, low = "#006BFF", mid = "#FFF100", high = "#E52020", na.value = "#3C3D37", 
                        limits = c(range(meta$YAP_S94A1)[[1]], range(meta$YAP_S94A1)[[2]])) 
pd3 <- ggplot_build(p3)$data[[1]]

meta <- meta %>% mutate(color = ifelse(cell.types == "Malignant", YAP_5SA1, NA))
p4 <- ggplot(meta, aes(x = x, y = y)) + geom_point(aes(color = color), size = 1) + 
  scale_color_gradient2(midpoint = 0, low = "#006BFF", mid = "#FFF100", high = "#E52020", na.value = "#3C3D37", 
                        limits = c(range(meta$YAP_5SA1)[[1]], range(meta$YAP_5SA1)[[2]])) 
pd4 <- ggplot_build(p4)$data[[1]]
meta$color <- NULL

# Append Colour 
meta$pd1 <- pd1$colour
meta$YAP1_colour <- pd2$colour
meta$YAP1_S94A_colour <- pd3$colour
meta$YAP1_5SA_colour <- pd4$colour
meta <- meta %>% mutate(YAP1_colour = ifelse(pd1 == "white", YAP1_colour, pd1))
meta <- meta %>% mutate(YAP1_S94A_colour = ifelse(pd1 == "white", YAP1_S94A_colour, pd1))
meta <- meta %>% mutate(YAP1_5SA_colour = ifelse(pd1 == "white", YAP1_5SA_colour, pd1))
meta <- meta[colnames(seurat), ]
seurat@meta.data <- meta

saveRDS(seurat, "ST_Discovery_so.rds")

## Plot ------------------------------------------------------------------------
id <- unique(meta$samples)
folder <- "~/Desktop/YAP1/plot/10xXenium/quickView/Discovery/"

for(sample in id){
  meta_sub <- meta %>% dplyr::filter(meta$samples %in% sample)
  
  # YAP1
  p <- ggplot(meta_sub, aes(x = x, y = y, color = YAP1_colour)) + 
    geom_point(size = 2, show.legend = FALSE) + scale_color_identity() + theme_minimal() + 
    labs(title = sample, subtitle = unique(meta_sub$sample_group), color = "Score", caption = "YAP1") 
  ggsave(paste0(folder, "YAP1/", sample, ".png"), p, width = 6, height = 4, dpi = 600, units = "in")
  
  # YAP1 S5A
  p <- ggplot(meta_sub, aes(x = x, y = y, color = YAP1_5SA_colour)) + 
    geom_point(size = 2, show.legend = FALSE) + scale_color_identity() + theme_minimal() + 
    labs(title = sample, subtitle = unique(meta_sub$sample_group), color = "Score", caption = "YAP1 S8A") 
  ggsave(paste0(folder, "YAP1 5SA/", sample, ".png"), p, width = 6, height = 4, dpi = 600, units = "in")
  
  p <- ggplot(meta_sub, aes(x = x, y = y, color = YAP1_S94A_colour)) + 
    geom_point(size = 2, show.legend = FALSE) + scale_color_identity() + theme_minimal() + 
    labs(title = sample, subtitle = unique(meta_sub$sample_group), color = "Score", caption = "YAP1 S94A") 
  ggsave(paste0(folder, "YAP1 S94A/", sample, ".png"), p, width = 6, height = 4, dpi = 600, units = "in")
}










