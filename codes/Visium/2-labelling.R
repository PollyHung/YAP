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

# ## Read in Reference 
# if(!file.exists("data/reference/CARD_REF.RData")){
#   sc_count <- readRDS("data/reference/GSE180661.rds")
#   sc_meta <- read.delim("data/reference/refDB/9606_map.tsv") %>% 
#     dplyr::select(cell_id, cell_type, patient_id) %>% 
#     group_by(cell_type) %>% slice_sample(prop = 0.4) %>% ungroup %>% 
#     dplyr::filter(cell_type != "Other") %>% as.data.frame() ## Slice 40%
#   rownames(sc_meta) <- sc_meta$cell_id
#   sc_count <- sc_count[, sc_meta$cell_id]
# } else {
#   load("data/reference/CARD_REF.RData")
# }

## Create CARD object for each sample -------------------------------------------
for(i in directories){
  setwd(i)
  
  ## Extract Spatial Objects 
  spatial <- readRDS("output/seurat.rds")
  spatial_count <- spatial@assays$Spatial$counts
  spatial_location <- spatial@images$slice1@coordinates %>% dplyr::select(row, col)
  colnames(spatial_location) <- c("x", "y")
  
  ## Build CARD
  CARD_obj = createCARDObject(sc_count = sc_count,
                              sc_meta = sc_meta,
                              spatial_count = spatial_count,
                              spatial_location = spatial_location,
                              ct.varname = "cell_type",
                              ct.select = unique(sc_meta$cell_type),
                              sample.varname = "patient_id",
                              minCountGene = 100,
                              minCountSpot = 5) 
  CARD_obj = CARD_deconvolution(CARD_object = CARD_obj)
  
  ## Plot 
  colors = c("grey","grey","grey","#F0027F", "grey", "#7FC97F", "#386CB0", "grey", "#FFD92F")
  p1 <- CARD.visualize.pie(proportion = CARD_obj@Proportion_CARD,
                           spatial_location = CARD_obj@spatial_location, 
                           colors = colors, radius = 0.7) 
  png("plots/card-deconvolute.png", width = 10, height = 10, units = "in", res = 300)
  print(p1)
  dev.off()
  
  p2 <- CARD.visualize.prop(proportion = CARD_obj@Proportion_CARD,        
                            spatial_location = CARD_obj@spatial_location, 
                            ct.visualize = c("Ovarian.cancer.cell", "Fibroblast"),                 
                            colors = c("lightblue","lightyellow","red"), 
                            NumCols = 4,                   
                            pointSize = 0.5)                            
  png("plots/card-deconvolute-cancer&fibro.png", width = 10, height = 5, units = "in", res = 300)
  print(p2)
  dev.off()
  
  ## save CARD 
  saveRDS(CARD_obj, "output/CARD.rds")
  
  ## Clean up 
  rm(list = setdiff(ls(), c("sc_count", "sc_meta", "directories")))
  gc()
}




## Work on Color ----------------------------------------------------------
i <- directories[[8]]
# for(i in directories){
  
  setwd(i)
  spatial <- readRDS("output/seurat.rds")
  CARD <- readRDS("output/CARD.rds")

  CARD <- CARD@Proportion_CARD %>% as.data.frame()
  CARD.round <- round(CARD, digits = 3)
  CARD.round <- CARD.round[rownames(spatial@meta.data), ]
  colnames(CARD.round) <- paste0("CARD.", colnames(CARD.round))

  colour <- brewer.pal(n = ncol(CARD.round), name = "Set3")
  names(colour) <- colnames(CARD.round)

  spatial@meta.data <- merge(spatial@meta.data, CARD.round, by = "row.names") %>%
    tibble::column_to_rownames(var = "Row.names")
  SpatialFeaturePlot(spatial, features = colnames(CARD.round), image.alpha = 0, 
                     alpha = 1, pt.size.factor = 2.4, keep.scale = "all")
  ggsave("plots/card-deconvolute-mapped-selected.png", width = 12, height = 12, units = "in", dpi = 600)
  
  saveRDS(spatial, "output/labelled.rds")
# }
















