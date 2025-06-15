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
library(future.apply)
library(progressr)
set.seed(123)

setwd("~/Desktop/Collection/spatial-seq/10xXenium/Seurat/")

## Processing ------------------------------------------------------------------
seurat <- readRDS("ST_Discovery_so.rds")

# Add Stiffness Scores 
stiffness <- c("COL11A1", "COMP", "FN1", "VCAN", "CTSB") %>% intersect(rownames(seurat))
seurat <- AddModuleScore(seurat, features = list(stiffness), name = "stiffness", ctrl = 40, slot = "scale.data")
seurat$stiffness <- seurat$stiffness1; seurat$stiffness1 <- NULL

## Use Mean Stiffness By Sample, > 0 is highly stiffed, < 0 is not stiffed 
mean_stiffness <- seurat@meta.data %>%
  group_by(samples) %>%        
  summarize(mean_stiffness = mean(stiffness, na.rm = TRUE)) 

seurat$sample_group <- ifelse(seurat$samples %in% mean_stiffness$samples[mean_stiffness$mean_stiffness > 0],
                              "high_stiffness", "low_stiffness")
ggplot(seurat@meta.data, aes(x = sample_group, y = stiffness)) +
  geom_violin(fill = "lightblue", alpha = 0.5) + 
  geom_boxplot(width = 0.1, color = "black", fill = "white", alpha = 1, outlier.size = 0.1) + 
  labs(title = "Stiffness", x = "Sample Group", y = "Stiffness") +
  theme_minimal() + 
  stat_compare_means(method = "t.test", label = "p.format")
ggsave("~/Desktop/YAP1/data/10xXenium/stiffness_discovery.png", width = 3, height = 4, units = "in", dpi = 600)

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


## Overview --------------------------------------------------------------------
meta <- seurat@meta.data
meta$cells <- rownames(meta)
cancer <- meta %>% dplyr::filter(cell.types == "Malignant")
wilcox.test(cancer$YAP1[cancer$sample_group == "high_stiffness"], cancer$YAP1[cancer$sample_group == "low_stiffness"])
wilcox.test(cancer$YAP_5SA1[cancer$sample_group == "high_stiffness"], cancer$YAP_5SA1[cancer$sample_group == "low_stiffness"])
wilcox.test(cancer$YAP_S94A1[cancer$sample_group == "high_stiffness"], cancer$YAP_S94A1[cancer$sample_group == "low_stiffness"])


## Plot ------------------------------------------------------------------------
# Make Colour 
p <- ggplot(meta, aes(x = x, y = y)) + 
  geom_point(aes(color = cell.types), size = 1, alpha = 1) +
  scale_color_manual(values = c("Fibroblast" = "#00FF9C", "Monocyte" = "#A27B5C", "TNK.cell" = "#3F4F44",
                                "B.cell" = "#2C3930", "Endothelial" = "#B6CBBD", "Mast.cell" = "#DCD7C9", "Malignant" = "white")) 
pd1 <- ggplot_build(p)$data[[1]]

meta <- meta %>% mutate(color = ifelse(cell.types == "Malignant", YAP1, NA))
p2 <- ggplot(meta, aes(x = x, y = y)) + geom_point(aes(color = color), size = 1) + 
  scale_color_gradient2(midpoint = 0, low = "#006BFF", mid = "#FFF100", high = "#E52020", na.value = "#3C3D37", 
                        limits = c(range(meta$YAP_S94A1)[[1]], range(meta$YAP_S94A1)[[2]])) 
pd2 <- ggplot_build(p2)$data[[1]]

meta <- meta %>% mutate(color = ifelse(cell.types == "Malignant", YAP_S94A1, NA))
p3 <- ggplot(meta, aes(x = x, y = y)) + geom_point(aes(color = color), size = 1) + 
  scale_color_gradient2(midpoint = 0, low = "#006BFF", mid = "#FFF100", high = "#E52020", na.value = "#3C3D37", 
                        limits = c(range(meta$YAP_S94A1)[[1]], range(meta$YAP_S94A1)[[2]])) 
pd3 <- ggplot_build(p3)$data[[1]]

meta <- meta %>% mutate(color = ifelse(cell.types == "Malignant", YAP_5SA1, NA))
p4 <- ggplot(meta, aes(x = x, y = y)) + geom_point(aes(color = color), size = 1) + 
  scale_color_gradient2(midpoint = 0, low = "#006BFF", mid = "#FFF100", high = "#E52020", na.value = "#3C3D37", 
                        limits = c(range(meta$YAP_S94A1)[[1]], range(meta$YAP_S94A1)[[2]])) 
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

## Use EBI to plot the real plot 
seg_path <- paste0("~/Desktop/Collection/spatial-seq/10xXenium/Discovery/", sample, "/", sample, "_whole-cell_03.csv")
cellseg = read.csv(seg_path)
colnames(cellseg) <- unlist(lapply(colnames(cellseg), function(x) {
  strsplit(x, split = "X0.")[[1]][2]
}))
colnames(cellseg)[1] <- "0"
cellmask <- cellseg
cellmask[cellmask != 0] <- 0
if("black" == "white") {
  cellmask[cellmask == 0] <- 1
} else {
  cellmask[cellmask == 0] <- 0
}
cellmask <- EBImage::Image(as.matrix(cellmask))
cellmask <- EBImage::channel(cellmask, "rgb")

# meta_sub <- meta %>% dplyr::filter(samples == sample)
# color_mapping <- setNames(meta_sub$YAP1_colour, rownames(meta_sub))

plan(multisession) 

# 1. Precompute all cell indices and colors upfront
cell_index_map <- setNames(
  as.integer(gsub(".*_c", "", rownames(meta_sub))),
  rownames(meta_sub)
)

# 2. Vectorized color conversion for ALL cells
rgb_matrix <- t(col2rgb(meta_sub$YAP1_colour)) / 255
rownames(rgb_matrix) <- rownames(meta_sub)

# 3. Parallel processing with progress bar
with_progress({
  p <- progressor(along = rownames(meta_sub))
  
  future_lapply(rownames(meta_sub), function(cell_id) {
    p() # Update progress bar
    
    # Get precomputed values
    cell_num <- cell_index_map[[cell_id]]
    rgb_vals <- rgb_matrix[cell_id, ]
    
    # Create mask once
    cell_mask <- cellseg == cell_num
    
    # Modify RGB channels in single operation
    cellmask@.Data[,,1][cell_mask] <- rgb_vals[1]
    cellmask@.Data[,,2][cell_mask] <- rgb_vals[2]
    cellmask@.Data[,,3][cell_mask] <- rgb_vals[3]
    
    NULL # Return nothing to minimize memory
  }, future.seed = TRUE)
})

EBImage::writeImage(cellmask, "~/Desktop/custom_colored_mask.png")


