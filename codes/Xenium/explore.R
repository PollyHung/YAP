library(Seurat)
library(dplyr)
library(magrittr)
library(ggplot2)
library(ggpubr)
library(DescTools)


## Round 1
test <- readRDS("~/Desktop/Collection/spatial-seq/10xXenium/Seurat/ST_Test1_so.rds")
meta <- test@meta.data


## Fibroblast Stiffness new Definition 
cell_types <- table(meta$samples, meta$cell.types) %>% as.data.frame.matrix()
cell_types$all_cells <- rowSums(cell_types)
cell_types$Fibroblast_pct <- (cell_types$Fibroblast/cell_types$all_cells)*100
cell_types$cluster <- ifelse(cell_types$Fibroblast_pct > 50, "high_stiff", "low_stiff")
meta$clusterByFibroCount <- ifelse(test$samples %in% rownames(cell_types)[cell_types$cluster == "high_stiff"], "high_stiff", "low_stiff")


## Heterogeneity per sample by Gini Coefficient 
df <- data.frame(row.names = unique(meta$samples), 
                 YAP1 = rep(NA, length(unique(meta$samples))), 
                 YAP1_5SA = rep(NA, length(unique(meta$samples))), 
                 YAP1_S94A = rep(NA, length(unique(meta$samples))))
for(i in unique(meta$samples)){
  meta_sub <- meta %>% dplyr::filter(samples == i) %>% dplyr::filter(cell.types == "Malignant")
  df[i, "YAP1"] <- Gini(meta_sub$YAP1_norm)
  df[i, "YAP1_5SA"] <- Gini(meta_sub$YAP_5SA1_norm)
  df[i, "YAP1_S94A"] <- Gini(meta_sub$YAP_S94A1_norm)
}
df2 <- meta %>% dplyr::select(samples, clusterByFibroCount) %>% distinct()
df <- merge(df2, df, by.x = "samples", by.y = "row.names")
colnames(df) <- c("samples", "clusterByFibroCount", "YAP1.hetero", "YAP1_5SA.hetero", "YAP1_S94A.hetero")
ggplot(df, aes(x = clusterByFibroCount, y = YAP1_5SA.hetero)) + geom_boxplot() + 
  stat_compare_means(method = 't.test') + 
  labs(title = "Cancer Cell Heterogeneity", 
       subtitle = "Heterogeneity by YAP1 5SA 
by Gini Coefficient") +
  theme_bw()
ggsave("plot/Cancer Heterogeneity by YAP1 5SA between High stiffness and low stiffness.png", 
       width = 3, height = 4, units = "in", dpi = 600)



## Cancer YAP1 5SA score by Stiffness 
df2 <- meta %>% 
  dplyr::filter(cell.types == "Malignant") %>% 
  dplyr::group_by(samples) %>% 
  dplyr::summarise(YAP_5SA1_avg = mean(YAP_5SA1))
df3 <- meta %>% dplyr::select(samples, clusterByFibroCount) %>% distinct
df2 <- merge(df2, df3, by = "samples")
ggplot(df2, aes(x = clusterByFibroCount, y = YAP_5SA1_avg)) + geom_boxplot() + 
  stat_compare_means(method = 't.test') + 
  labs(title = "YAP1 5SA score by Stiffness") +
  theme_bw()
ggsave("plot/YAP1 5SA score by stiffness score.png", 
       width = 3, height = 3.5, units = "in", dpi = 600)



## Copy Number 
tcga <- read.delim("~/Desktop/table.tsv")
tcga$chromosome <- unname(unlist(lapply(strsplit(tcga$Cytoband, split = "q"), function(x){x[[1]]})))
tcga <- tcga %>% dplyr::filter(chromosome == "11") %>% dplyr::filter(Spearman.s.Correlation > 0.5)
test <- AddModuleScore(test, features = list(tcga$Correlated.Gene), name = "chr11q22")
tcga <- test@meta.data %>% 
  dplyr::group_by(samples) %>% 
  summarise(copyNumber = mean(chr11q221))
tcga$status <- ifelse(tcga$copyNumber > 0, "amp", "wt")
df3 <- merge(df2, tcga, by = "samples") 
df3 <- merge(df, df3, by = c("samples", "clusterByFibroCount"))
ggplot(df3, aes(x = status, y = YAP1_5SA.hetero)) + 
  geom_boxplot() + 
  stat_compare_means(method = 't.test') + 
  labs(title = "Cancer Heterogeneity by YAP1 CN", 
       subtitle = "YAP1 Copy Number 
  inferred from 11q22") + 
  theme_bw()
ggsave("plot/Heterogenetity by YAP1 Copy Number.png", 
       width = 3, height = 3.5, units = "in", dpi = 600) 


temp <- meta %>% dplyr::select(samples, clusterByFibroCount) %>% distinct()
df2 <- merge(temp, df, by.x = "samples", by.y = "row.names")
ggplot(df2, aes(x = clusterByFibroCount, y = YAP1_5SA)) + geom_boxplot() + 
  stat_compare_means(method = 't.test') + labs(title = "heterogeneity by YAP1 5SA")


tcga <- read.delim("~/Desktop/table.tsv")
tcga$chromosome <- unname(unlist(lapply(strsplit(tcga$Cytoband, split = "q"), function(x){x[[1]]})))
tcga <- tcga %>% dplyr::filter(chromosome == "11") %>% dplyr::filter(Spearman.s.Correlation > 0.5)

test <- AddModuleScore(test, features = list(tcga$Correlated.Gene), name = "11q22")

case <- test@meta.data %>% 
  dplyr::group_by(samples) %>% 
  summarise(copyNumber = mean(`11q221`))

case$status <- ifelse(case$copyNumber > 0, "amp", "wt")
df2 <- merge(df2, case, by = "samples") 
ggplot(df2, aes(x = status, y = YAP1_5SA)) + geom_boxplot() + 
  stat_compare_means(method = 't.test') + labs(title = "heterogeneity by YAP1 5SA")


cancer2 <- cancer %>% 
  dplyr::group_by(samples) %>% 
  summarise(YAP_5SA1_avg = mean(YAP_5SA1))
cancer3 <- merge(cancer2, case, by = "samples")
ggplot(cancer3, aes(x = status, y = YAP_5SA1_avg)) + geom_boxplot() + 
  stat_compare_means(method = 't.test') + labs(title = "YAP1 5SA score")

ggplot(cancer, aes(x = status, y = YAP_5SA1_avg)) + geom_boxplot() + 
  stat_compare_means(method = 't.test') + labs(title = "YAP1 5SA score")



## YAP1 5SA average per patient between high stiffness and low stiffness groups 
## YAP1 5SA heterogeneity score between high stiffness and low stiffness groups 
## YAP1 5SA heterogeneity between amp YAP1 and wt YAP1




