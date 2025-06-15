library(Seurat)
library(dplyr)
library(magrittr)
library(tidyr)
library(ggplot2)
library(stringr)
library(CARD)
library(MuSiC)
library(RColorBrewer)

set.seed(20020208)

directories <- list.dirs("~/Desktop/YAP1/data/spatial/GSE211956", recursive = F)

get_cell_type <- function(x) {
  max_value <- max(x)
  max_indices <- which(x == max_value)
  max_names <- names(card_columns)[max_indices]
  
  # Check for ties and apply rules
  if (length(max_names) > 1) {
    if (any(grepl("Fibroblast", max_names))) {
      return("CARD.Fibroblast")
    } else if (any(grepl("cancer.cell", max_names))) {
      return(max_names[!grepl("cancer.cell", max_names)][1])  # Pick the first non-cancer cell type
    } else {
      return(sample(max_names, 1) )  # Randomly pick one of the tied cell types
    }
  }
  
  return(max_names[1])  # No ties, return the single max name
}


## Work on Color ----------------------------------------------------------
for(i in directories){
  
  setwd(i)
  spatial <- readRDS("output/labelled.rds")
  
  ## View the Distribution of CARD proportion 
  plot_df <- spatial@meta.data %>% select(starts_with("CARD.")) %>% 
    pivot_longer(cols = everything(), names_to = "Cell_Type", values_to = "Likelihood")
  ggplot(plot_df, aes(x = Cell_Type, y = Likelihood)) + geom_boxplot() + 
    labs(title = "Boxplot of Cell Type Likelihoods", x = "Cell Type", y = "Likelihood") +
    theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1))  
  ggsave("plots/dist-CARD-proportion.png", width = 4, height = 4, units = "in", dpi = 600)
  
  ## add label 
  card_cols <- grepl("CARD", colnames(spatial@meta.data))
  # spatial@meta.data[, card_cols] <- round(spatial@meta.data[, card_cols], digits = 1)
  card_columns <- spatial@meta.data %>% select(starts_with("CARD"))
  spatial@meta.data <- spatial@meta.data %>% mutate(cell_type = apply(card_columns, 1, get_cell_type))
  
  ## Plot 
  colour <- brewer.pal(9, name = "Pastel1")
  names(colour) <- c("CARD.Fibroblast", "CARD.Ovarian.cancer.cell", "CARD.Endothelial.cell",   
                     "CARD.T.cell", "CARD.B.cell", "CARD.Myeloid.cell", "CARD.Plasma.cell")
  SpatialDimPlot(spatial, group.by = "cell_type", pt.size.factor = 2.6, image.alpha = 0, cols = colour)
  ggsave("plots/card-deconvolute-mapped.png", width = 6, height = 4, units = "in", dpi = 600)

  ## Save RDS
  saveRDS(spatial, "output/labelled.rds")
}


## Extract metadata 
meta <- list()

meta <- lapply(directories, function(x){
  spatial <- readRDS(file.path(x, "output/labelled.rds"))
  spatial <- spatial@meta.data
  spatial$cells <- paste0(spatial$sample, "_", spatial$cells)
  rownames(spatial) <- spatial$cells
  return(spatial)
})

meta_df <- do.call(rbind, meta)
write.csv(meta_df, "data/spatial/collapsed/labelled_metadata.csv")
















