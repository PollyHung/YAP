library(Seurat)
library(dplyr)
library(magrittr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(Matrix)
library(stringr)

setwd("~/Desktop/YAP1/")


## Organize the Seurat Objects 
files <- list.dirs("/Users/polly_hung/Desktop/YAP1/data/spatial/GSE211956", recursive = F)
files <- files[2:8]

for(i in files){
  
  setwd(file.path(i, "output"))
  
  # ## Create Directories 
  # if(!dir.exists("metadata")){(dir.create("metadata"))}
  # if(!dir.exists("expression")){dir.create("expression")}
  # if(!dir.exists("reduction")){dir.create("reduction")}
  # if(!dir.exists("rds")){dir.create("rds")}
  
  ## Read in RDS
  labelled2 <- readRDS("rds/labelled2.rds")
  labelled2_meta <- labelled2@meta.data

  ## Save Metadata
  writexl::write_xlsx(labelled2_meta, "metadata/metadata.xlsx")

  ## Save Counts
  counts <- labelled2@assays$SCT$counts
  log2_norm <- labelled2@assays$SCT$data
  scaled <- labelled2@assays$SCT$scale.data
  writeMM(counts, "expression/counts.mtx")
  writeMM(log2_norm, "expression/data.mtx")
  writeMM(log2_norm, "expression/scale.data.mtx")

  ## Save Reductions
  pca <- labelled2@reductions$pca
  saveRDS(pca, file = "reduction/pca.rds")
  umap <- labelled2@reductions$umap
  saveRDS(pca, file = "reduction/umap.rds")
}






























