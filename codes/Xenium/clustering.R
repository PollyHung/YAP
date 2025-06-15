library(magrittr)
library(dplyr)
library(Seurat)
library(umap)
library(ggplot2)


setwd("~/Desktop/Collection/spatial-seq/10xXenium/")

## Read in 
seurat <- readRDS("Seurat/ST_Test1_so.rds")
YAP1 <- readxl::read_xlsx("~/Desktop/YAP1/data/YAP1.xlsx", col_names = F)
YAP <- YAP1[1, 3:202] %>% unlist %>% unname %>% intersect(rownames(seurat))
YAP_S94A <- YAP1[2, 3:202] %>% unlist %>% unname %>% intersect(rownames(seurat))
YAP_5SA <- YAP1[3, 3:202] %>% unlist %>% unname %>% intersect(rownames(seurat))

## Fetch Data 
YAP1 <- FetchData(seurat, vars = c("samples", "sample_group", "cell.types", "stiffness", YAP), layer = "scale.data")
YAP1_S94A <- FetchData(seurat, vars = c("samples", "sample_group", "cell.types", "stiffness", YAP_S94A), layer = "scale.data")
YAP1_5SA <- FetchData(seurat, vars = c("samples", "sample_group", "cell.types", "stiffness", YAP_5SA), layer = "scale.data")

## Summarise 
YAP1_cancer <- YAP1 %>% 
  dplyr::filter(cell.types == "Malignant") %>% 
  dplyr::group_by(samples) %>% 
  summarise(across(starts_with("SV2B"):last_col(), mean, na.rm = TRUE)) %>% 
  tibble::column_to_rownames(var = "samples")
YAP1_group <- distinct(YAP1[, c("samples", "sample_group")])
rownames(YAP1_group) <- YAP1_group$samples

YAP1_5SA_cancer <- YAP1_5SA %>% 
  dplyr::filter(cell.types == "Malignant") %>% 
  dplyr::group_by(samples) %>% 
  summarise(across(starts_with("SV2B"):last_col(), mean, na.rm = TRUE)) %>% 
  tibble::column_to_rownames(var = "samples")
YAP1_5SA_group <- distinct(YAP1_5SA[, c("samples", "sample_group")])
rownames(YAP1_5SA_group) <- YAP1_5SA_group$samples

YAP1_S94A_cancer <- YAP1_S94A %>% 
  dplyr::filter(cell.types == "Malignant") %>% 
  dplyr::group_by(samples) %>% 
  summarise(across(starts_with("SPARCL1"):last_col(), mean, na.rm = TRUE)) %>% 
  tibble::column_to_rownames(var = "samples")
YAP1_S94A_group <- distinct(YAP1_S94A[, c("samples", "sample_group")])
rownames(YAP1_S94A_group) <- YAP1_S94A_group$samples
  

## Run UAMP
YAP1_umap <- umap(YAP1_cancer)
YAP1_umap <- data.frame(YAP1_umap$layout)
YAP1_umap <- merge(YAP1_group, YAP1_umap, by = "row.names") %>% tibble::column_to_rownames(var = "Row.names")
ggplot(YAP1_umap, aes(x = X1, y = X2, label = samples, colour = sample_group)) + 
  geom_point() + labs(title = "Cluster Samples by Cancer Cell YAP1 Expression", 
                      subtitle = "Coloured by stiffness; UMAP reduction", 
                      caption = "Fibroblast stiffness does not differnetiate 
                      YAP1 expression in cancer cells; each dot is a sample") + theme_bw()
ggsave("~/Desktop/YAP1/plot/summary/umap.test.yap1.png", width = 6, height = 4, units = "in", dpi = 600)

YAP1_5SA_umap <- umap(YAP1_5SA_cancer)
YAP1_5SA_umap <- data.frame(YAP1_5SA_umap$layout)
YAP1_5SA_umap <- merge(YAP1_5SA_group, YAP1_5SA_umap, by = "row.names") %>% tibble::column_to_rownames(var = "Row.names")
ggplot(YAP1_5SA_umap, aes(x = X1, y = X2, label = samples, colour = sample_group)) + 
  geom_point() + labs(title = "Cluster Samples by Cancer Cell YAP1 5SA Expression", 
                      subtitle = "Coloured by stiffness; UMAP reduction", 
                      caption = "Fibroblast stiffness does not differnetiate 
                      YAP1 5SA expression in cancer cells; each dot is a sample") + theme_bw()
ggsave("~/Desktop/YAP1/plot/summary/umap.test.yap1_5sa.png", width = 6, height = 4, units = "in", dpi = 600)

YAP1_S94A_umap <- umap(YAP1_S94A_cancer)
YAP1_S94A_umap <- data.frame(YAP1_S94A_umap$layout)
YAP1_S94A_umap <- merge(YAP1_S94A_group, YAP1_S94A_umap, by = "row.names") %>% tibble::column_to_rownames(var = "Row.names")
ggplot(YAP1_S94A_umap, aes(x = X1, y = X2, label = samples, colour = sample_group)) + 
  geom_point() + labs(title = "Cluster Samples by Cancer Cell YAP1 S94A Expression", 
                      subtitle = "Coloured by stiffness; UMAP reduction", 
                      caption = "Fibroblast stiffness does not differnetiate 
                      YAP1 S94A expression in cancer cells; each dot is a sample") + theme_bw()
ggsave("~/Desktop/YAP1/plot/summary/umap.test.yap1_s94a.png", width = 6, height = 4, units = "in", dpi = 600)









