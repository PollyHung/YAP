library(Seurat)
library(dplyr)
library(magrittr)
library(ggplot2)

## Read in RDS
# fibroblast <- readRDS("data/scRNA/fibroblasts.rds")

## YAP1 signatures  
YAP_sigs <- readxl::read_xlsx("data/YAP1.xlsx", col_names = F)
YAP <- YAP_sigs[1, 3:202] %>% unlist %>% unname %>% intersect(rownames(fibroblast))
# YAP_S94A <- YAP_sigs[2, 3:202] %>% unlist %>% unname %>% intersect(rownames(fibroblast))
# YAP_5SA <- YAP_sigs[3, 3:202] %>% unlist %>% unname %>% intersect(rownames(fibroblast))
YAP_S94A <- read.delim("data/scRNA/YAP_S94A.txt", header = F) %>% unlist %>% unname
YAP_5SA <- read.delim("data/scRNA/YAP_5SA.txt", header = F) %>% unlist %>% unname

## Stiffness 
stiffness <- c("COL11A1", "COMP", "FN1", "VCAN", "CTSB", "COL1A1")

## add module scores 
fibroblast <- AddModuleScore(fibroblast, features = list(YAP), name = "YAP1")
fibroblast <- AddModuleScore(fibroblast, features = list(YAP_S94A), name = "YAP1S94A")
fibroblast <- AddModuleScore(fibroblast, features = list(YAP_5SA), name = "YAP15SA")
fibroblast <- AddModuleScore(fibroblast, features = list(stiffness), name = "stiffness")
fibroblast$YAP1 <- fibroblast$YAP11; fibroblast$YAP11 <- NULL
fibroblast$YAP1S94A <- fibroblast$YAP1S94A1; fibroblast$YAP1S94A1 <- NULL
fibroblast$YAP15SA <- fibroblast$YAP15SA1; fibroblast$YAP15SA1 <- NULL
fibroblast$stiffness <- fibroblast$stiffness1; fibroblast$stiffness1 <- NULL

## YAP1+ cells
fibroblast <- AddModuleScore(fibroblast, features = "YAP1", name = "YAP1_positive")
fibroblast$YAP1_positive <- fibroblast$YAP1_positive1; fibroblast$YAP1_positive1 <- NULL
fibroblast$YAP1_status <- ifelse(fibroblast$YAP1_positive > 0 & fibroblast$YAP1 > 0, "YAP1+", "YAP1-")

## the two mutants 
YAP1 <- subset(fibroblast, subset = YAP1_status == "YAP1+")
YAP1 <- RunPCA(YAP1, features = c(YAP_5SA, YAP_S94A))
YAP1 <- FindNeighbors(YAP1, dims = 1:15)
YAP1 <- FindClusters(YAP1, resolution = 0.2)
YAP1 <- RunUMAP(YAP1, dims = 1:10)

diff_exp <- FindAllMarkers(YAP1)
YAP_5SA_diff_exp <- diff_exp %>% dplyr::filter(gene %in% YAP_5SA)
YAP_5SA_diff_exp$diff_pct <- YAP_5SA_diff_exp$pct.1 - YAP_5SA_diff_exp$pct.2
YAP_5SA_diff_exp <- YAP_5SA_diff_exp %>% dplyr::filter(diff_pct > 0.1)

FeaturePlot(fibroblast, features = "rna_YAP1", split.by = "YAP1_status")
FeaturePlot(YAP1, features = c("YAP1S94A", "YAP15SA"), pt.size = 0.1)
DimPlot(YAP1, pt.size = 0.1)

YAP1$YAP1S94A_status <- ifelse(YAP1$RNA_snn_res.0.2 %in% c("2"), "YAP1-S94A+", "YAP1-S94A-")
DimPlot(YAP1, group.by = "YAP1S94A_status", pt.size = 0.1)

YAP1$YAP15SA_status <- ifelse(YAP1$RNA_snn_res.0.2 %in% c("0", "3"), "YAP1-5SA+", "YAP1-5SA-")
DimPlot(YAP1, group.by = "YAP15SA_status", pt.size = 0.1)

YAP1$YAP1_status <- ifelse(YAP1$RNA_snn_res.0.2 %in% c("2"), "YAP1-S94A+", "YAP1 wt")
YAP1$YAP1_status <- ifelse(YAP1$RNA_snn_res.0.2 %in% c("0", "3"), "YAP1-5SA+", YAP1$YAP1_status)

DimPlot(YAP1, group.by = "YAP1_status", pt.size = 0.1)




YAP1S94A_diff_exp <- FindMarkers(YAP1, ident.1 = "YAP1-S94A+", "YAP1-S94A-", group.by = "YAP1S94A_status")
YAP1S94A_diff_exp <- YAP1S94A_diff_exp %>% dplyr::arrange(-avg_log2FC)
YAP1S94A_diff_exp$genes <- rownames(YAP1S94A_diff_exp)
write.table(YAP1S94A_diff_exp[, c("genes", "avg_log2FC")], "~/Desktop/YAP1S94A.rnk", sep = "\t", quote = F, row.names = F, col.names = F)

