# ------------------------------------------------------------------------
# Title: Preprocessing of the 10x Xenium Dataset
# Description: This script performs preprocessing of the 10x Xenium dataset, 
#              specifically calculating YAP1 and associated mutant scores, 
#              fibroblast stiffness scores, and the copy number of YAP1. 
#              Additionally, it adds color coding to the datasets based on 
#              binned YAP1 and associated mutant scores for enhanced 
#              visualization.
# Author: Polly
# Date: 18 March 2025
# ------------------------------------------------------------------------


## Library Packages, Define Functions ------------------------------------
library(Seurat)
library(dplyr)
library(magrittr)
library(DescTools)
library(scales)
library(ggplot2)
background = "black"
low_qc_color = 0

## Processing (only need to be done once) --------------------------------
obj <- readRDS("~/Desktop/Collection/spatial-seq/10xXenium/Seurat/ST_Test1_so.rds")
assay5 <- as(obj[["RNA"]], Class = "Assay5")
obj <- CreateSeuratObject(assay5, meta.data = obj@meta.data)
obj <- JoinLayers(obj)
obj$cell.types <- gsub("_LC", "", obj$cell.types)


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


## Rename Bins to set up for colouring                
obj@meta.data <- obj@meta.data %>%
  mutate(YAP1.bins.named = case_when(
    YAP1.bins == "[-0.942,-0.367]" ~ "cancer1",
    YAP1.bins == "(-0.367,-0.276]" ~ "cancer2",
    YAP1.bins == "(-0.276,-0.214]" ~ "cancer3",
    YAP1.bins == "(-0.214,-0.163]" ~ "cancer4",
    YAP1.bins == "(-0.163,-0.115]" ~ "cancer5",
    YAP1.bins == "(-0.115,-0.0651]" ~ "cancer6",
    YAP1.bins == "(-0.0651,-0.00987]" ~ "cancer7",
    YAP1.bins == "(-0.00987,0.0619]" ~ "cancer8",
    YAP1.bins == "(0.0619,0.166]" ~ "cancer9",
    YAP1.bins == "(0.166,1.25]" ~ "cancer10",
    TRUE ~ ifelse(cell.types == "Fibroblast", "fibroblast", "others")
  ))
obj@meta.data <- obj@meta.data %>%
  mutate(YAP1_5SA.bins.named = case_when(
    YAP1_5SA.bins == "[-0.9,-0.351]" ~ "cancer1",
    YAP1_5SA.bins == "(-0.351,-0.266]" ~ "cancer2",
    YAP1_5SA.bins == "(-0.266,-0.205]" ~ "cancer3",
    YAP1_5SA.bins == "(-0.205,-0.155]" ~ "cancer4",
    YAP1_5SA.bins == "(-0.155,-0.108]" ~ "cancer5",
    YAP1_5SA.bins == "(-0.108,-0.0618]" ~ "cancer6",
    YAP1_5SA.bins == "(-0.0618,-0.0104]" ~ "cancer7",
    YAP1_5SA.bins == "(-0.0104,0.0518]" ~ "cancer8",
    YAP1_5SA.bins == "(0.0518,0.145]" ~ "cancer9",
    YAP1_5SA.bins == "(0.145,0.821]" ~ "cancer10",
    TRUE ~ ifelse(cell.types == "Fibroblast", "fibroblast", "others")
  ))
obj@meta.data <- obj@meta.data %>%
  mutate(YAP1_S94A.bins.named = case_when(
    YAP1_S94A.bins == "[-0.985,-0.347]" ~ "cancer1",
    YAP1_S94A.bins == "(-0.347,-0.255]" ~ "cancer2",
    YAP1_S94A.bins == "(-0.255,-0.191]" ~ "cancer3",
    YAP1_S94A.bins == "(-0.191,-0.136]" ~ "cancer4",
    YAP1_S94A.bins == "(-0.136,-0.0842]" ~ "cancer5",
    YAP1_S94A.bins == "(-0.0842,-0.0277]" ~ "cancer6",
    YAP1_S94A.bins == "(-0.0277,0.0342]" ~ "cancer7",
    YAP1_S94A.bins == "(0.0342,0.117]" ~ "cancer8",
    YAP1_S94A.bins == "(0.117,0.253]" ~ "cancer9",
    YAP1_S94A.bins == "(0.253,1.39]" ~ "cancer10",
    TRUE ~ ifelse(cell.types == "Fibroblast", "fibroblast", "others")
  ))


## Plot Each Point By Label 
meta <- obj@meta.data
meta_sub <- meta %>% dplyr::filter(samples == "SMI6K_2_F00033")

cellseg <- read.csv("~/Desktop/Collection/spatial-seq/10xXenium/Test/segments/20240215_193538_S1_C902_P99_N99_F00033.TIF.csv")
colnames(cellseg) <- unlist(lapply(colnames(cellseg), function(x) {
  strsplit(x, split = "X0.")[[1]][2]
}))
colnames(cellseg)[1] <- "0"
cellmask <- cellseg
cellmask[cellmask != 0] <- low_qc_color
if (background == "white") {
  cellmask[cellmask == 0] <- 1
} else {
  cellmask[cellmask == 0] <- 0
}
cellmask <- EBImage::Image(as.matrix(cellmask))
cellmask <- EBImage::channel(cellmask, "rgb")

celltypes <- meta_sub$YAP1_5SA.bins.named
names(celltypes) <- rownames(meta_sub)
levels <- setdiff(unique(celltypes), cont_field)

cell_2_rgb <- list()
cell_2_rgb[["others"]] = c(unname(col2rgb("#3C3D37")[,1]))
cell_2_rgb[["fibroblast"]] = c(unname(col2rgb("darkgrey")[,1]))
cell_2_rgb[["cancer1"]] = c(unname(col2rgb("#438381")[,1]))
cell_2_rgb[["cancer2"]] = c(unname(col2rgb("#74B86D")[,1]))
cell_2_rgb[["cancer3"]] = c(unname(col2rgb("#8DBF7A")[,1]))
cell_2_rgb[["cancer4"]] = c(unname(col2rgb("#A8C68B")[,1]))
cell_2_rgb[["cancer5"]] = c(unname(col2rgb("#C4D8A4")[,1]))
cell_2_rgb[["cancer6"]] = c(unname(col2rgb("#EBE4BA")[,1]))
cell_2_rgb[["cancer7"]] = c(unname(col2rgb("#C8B79A")[,1]))
cell_2_rgb[["cancer8"]] = c(unname(col2rgb("#B78B72")[,1]))
cell_2_rgb[["cancer9"]] = c(unname(col2rgb("#A85A4D")[,1]))
cell_2_rgb[["cancer10"]] = c(unname(col2rgb("#953126")[,1]))


for (type in levels) {
  cellids = names(celltypes[celltypes == type])
  cellidx = unlist(lapply(cellids, function(x) as.integer(strsplit(x,
                                                                   split = "c")[[1]][2])))
  celltype_mask <- t(apply(cellseg, 1, function(x) {
    x %in% cellidx
  }))
  celltype_rgb_1 <- cellmask@.Data[, , 1]
  celltype_rgb_1[celltype_mask] <- as.numeric(cell_2_rgb[[type]][1])/255
  cellmask[, , 1] <- celltype_rgb_1
  celltype_rgb_2 <- cellmask@.Data[, , 2]
  celltype_rgb_2[celltype_mask] <- as.numeric(cell_2_rgb[[type]][2])/255
  cellmask[, , 2] <- celltype_rgb_2
  celltype_rgb_3 <- cellmask@.Data[, , 3]
  celltype_rgb_3[celltype_mask] <- as.numeric(cell_2_rgb[[type]][3])/255
  cellmask[, , 3] <- celltype_rgb_3
}

EBImage::writeImage(cellmask, files = "~/Desktop/SMI6K_2_F00033.png")



test2 <- readRDS("~/Desktop/Collection/spatial-seq/10xXenium/Seurat/ST_Test2_HGSC1_so.rds")




