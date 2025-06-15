library(dplyr)
library(magrittr)
library(Seurat)
library(harmony)

setwd("~/Desktop/Collection/spatial-seq/10xXenium/SeuratObj")

spatial <- readRDS("ST_Test2_HGSC2_so.rds")

assay5 <- as(spatial[["RNA"]], Class = "Assay5")
spatial <- CreateSeuratObject(assay5, meta.data = spatial@meta.data)
spatial <- JoinLayers(spatial)

spatial[["percent.mt"]] <- PercentageFeatureSet(spatial, pattern = "^MT-")
spatial[["mitoRatio"]] <- spatial@meta.data$percent.mt/100
spatial[["log10GenesPerUMI"]] <- log10(spatial$nFeature_RNA)/log10(spatial$nCount_RNA)

spatial <- FindVariableFeatures(spatial, selection.method = "vst", nfeatures = nrow(spatial))

regev_lab_cell_cycle_genes <- read.delim("~/Desktop/regev_lab_cell_cycle_genes.txt", header=F) %>% unlist
s.genes <- intersect(rownames(spatial), regev_lab_cell_cycle_genes[1:43]) 
g2m.genes <- intersect(rownames(spatial), regev_lab_cell_cycle_genes[44:97])  
spatial <- CellCycleScoring(object = spatial,
                            s.features = s.genes,
                            g2m.features = g2m.genes,
                            set.ident = TRUE)

spatial <- RunPCA(spatial, features = VariableFeatures(object = spatial))

spatial <- harmony::RunHarmony(spatial,
                               group.by.vars = "samples",
                               reduction = "pca",
                               assay.use = "RNA",
                               reduction.save = "harmony",
                               lambda = 1)
spatial <- FindNeighbors(spatial, dims = 1:20, reduction = "harmony")
spatial <- FindClusters(object = spatial, resolution = c(0.1, 0.3, 0.5, 0.7), verbose = TRUE)
spatial <- RunUMAP(spatial, dims = 1:20, reduction="harmony")

saveRDS(object = spatial, file = "ST_Test2_HGSC2_so.rds")



