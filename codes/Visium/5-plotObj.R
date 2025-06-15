library(Seurat)
library(dplyr)
library(magrittr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(stringr)
library(scales)
library(viridis)

setwd("~/Desktop/YAP1/")


## Merged Sample Level ---------------------------------------------------------
meta <- read.delim("data/merged_visium_xenium2.meta", row.names = 1) %>% 
  dplyr::arrange(YAP1_5SA)
meta$YAP1_5SA.colour <- viridis::inferno(n = nrow(meta))

plot_df <- meta %>% dplyr::select(sample, Stiffness, YAP1_5SA.hetero) %>% distinct()
ggplot(plot_df, aes(x = Stiffness, y = YAP1_5SA.hetero)) + 
  geom_boxplot() + stat_compare_means(vjust = 3, method = "t.test") + theme_bw() + 
  labs(title = "Cancer Cell Heterogeneity", subtitle = "Heterogeneity by YAP1-5SA") + xlab("Fibroblast Stiffness") + ylab("YAP1-5SA") 

plot_df <- meta %>% dplyr::select(sample, Stiffness, YAP1_5SA) %>% 
  dplyr::group_by(sample, Stiffness) %>% 
  summarise(mean_YAP1_5SA = mean(YAP1_5SA, na.rm = TRUE), .groups = 'drop')
ggplot(plot_df, aes(x = Stiffness, y = mean_YAP1_5SA)) + 
  geom_boxplot() + stat_compare_means(vjust = 3, method = "t.test") + theme_bw() + 
  labs(title = "YAP1-5SA Activation Score") + xlab("Fibroblast Stiffness") + ylab("YAP1-5SA") 




## Each Sample Level -----------------------------------------------------------
GSM_cohort <- plot_df[grepl("GSM", plot_df$sample), ]

high_stiffness <- GSM_cohort$sample[GSM_cohort$Stiffness == "high_stiff"]
low_stiffness <- GSM_cohort$sample[GSM_cohort$Stiffness == "low_stiff"]

i <- low_stiffness[[6]]

for(i in c(high_stiffness, low_stiffness)){
  
  setwd(file.path("/Users/polly_hung/Desktop/YAP1/data/spatial/GSE211956", i))
  spatial <- readRDS("output/labelled2.rds")
  
  ## Use Viridis from all the samples, the normalized version 
  submeta <- spatial@meta.data 
  viridis <- meta %>% dplyr::filter(sample == i) %>% 
    dplyr::select(sample, YAP1_5SA.colour) 
  rownames(viridis) <- gsub(i, "", rownames(viridis))
  rownames(viridis) <- gsub("_", "", rownames(viridis))
  viridis$sample <- NULL
  viridis <- viridis[rownames(submeta), ]
  submeta <- cbind(submeta, viridis)
  
  ## colour non-cancer cells
  submeta <- submeta %>% dplyr::mutate(viridis = case_when(cell.types == "Fibroblast" ~ "#A7CF5D", 
                                                           cell.types == "Ovarian.cancer.cell" ~ viridis,  
                                                           TRUE ~ "darkgrey"))
  
  ## Add Spatial Coordinates 
  submeta <- submeta[rownames(spatial@images$slice1@coordinates), ]
  submeta$imagerow <- spatial@images$slice1@coordinates$imagerow
  submeta$imagecol <- spatial@images$slice1@coordinates$imagecol
  
  ## Plot 
  n = table(submeta$cell.types)[["Ovarian.cancer.cell"]]
  ggplot(submeta, aes(x = imagerow, y = imagecol, color = viridis)) + 
    geom_point(size =1.5) + 
    scale_color_identity() + theme_void() + coord_flip()+ scale_x_reverse() 
  ggsave("plots/overlayed-YAP1-5SA.png", width = 5, height = 5, units = "in", dpi = 600)
  
  ## save the plot
  write.table(submeta, "output/plot_setting.meta", sep = "\t", quote = F, row.names = T, col.names = T)
}


## Create a Colour Bar? 
merged <- readRDS("~/Desktop/YAP1/data/merged_visium_xenium2.rds")
meta <- meta[rownames(merged@meta.data), ]
merged@meta.data <- meta

FeaturePlot(merged, feature = "YAP1_5SA", ) + 
  scale_color_viridis(option = "inferno", direction = 1)



## YAP1-5SA score and fibroblast percentage 
# GSM <- meta %>% dplyr::filter(grepl("GSM", sample))
GSM <- meta %>% 
  dplyr::filter(cell.types == "Ovarian.cancer.cell") %>% 
  dplyr::group_by(sample) %>% summarise(mean(YAP1_5SA))

fibro <- as.data.frame.matrix(table(merged$sample, merged$cell.types))
fibro$All.Cells <- rowSums(fibro)
fibro$Fibroblast.Pct <- 100*(fibro$Fibroblast / fibro$All.Cells)
fibro <- fibro[GSM$sample, ]

plot_df <- data.frame(sample = GSM$sample, 
                      YAP1_5SA_mean = GSM$`mean(YAP1_5SA)`, 
                      fibroblast_pct = fibro$Fibroblast.Pct, 
                      group = ifelse(grepl("GSM", GSM$sample), "visium", "xenium"))


ggplot(plot_df, aes(x = fibroblast_pct, y = YAP1_5SA_mean, colour = group)) + 
  geom_point() + 
  geom_text(aes(label = plot_df$sample), vjust = -1, size = 3) + 
  geom_smooth(method = "lm", se = FALSE, colour = "#5F8B4C") + 
  xlim(-10, 120) + ylim(-0.14, 0.12) + 
  theme_bw() + 
  labs(title = "YAP1-5SA score in Cancer Cells Positively Correlates 
       with Increased Fibroblast Percentage in Tissue Section") + 
  xlab("Percent of Fibroblast in Tissue Section") + 
  ylab("Mean YAP1-5SA score in Cancer Cells")



## Each Sample Level -----------------------------------------------------------
GSM_cohort <- list.dirs("data/spatial/GSE211956", recursive = F)

for(i in GSM_cohort){
  
  setwd(file.path("/Users/polly_hung/Desktop/YAP1/", i))
  submeta <- read.delim("output/plot_setting.meta", row.names = 1)
  
  ## Plot 
  n = table(submeta$cell.types)[["Ovarian.cancer.cell"]]
  ggplot(submeta, aes(x = imagerow, y = imagecol, color = viridis)) + 
    geom_point(size =1.5) + 
    scale_color_identity() + theme_void() + coord_flip()+ scale_x_reverse() 
  ggsave("plots/overlayed-YAP1-5SA.png", width = 5, height = 5, units = "in", dpi = 600)
  
}












