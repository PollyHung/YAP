library(Seurat)
library(dplyr)
library(magrittr)
library(tidyr)
library(ggplot2)
library(stringr)

## Load and Build Seurat Object
id = list.dirs("data/spatial/GSE211956/", full.names = F, recursive = F)
lst = list()
for (i in 1:length(id)){
  temp = Seurat::Read10X(data.dir = paste0("data/spatial/GSE211956/", id[i], "/count_matrix"))
  lst[[i]] = CreateSeuratObject(counts = temp)
}
spatial = merge(lst[[1]], y = c(lst[[2]],lst[[3]],lst[[4]],lst[[5]],lst[[6]],lst[[7]],lst[[8]]),
                 add.cell.ids = c(id[1],id[2],id[3],id[4],id[5],id[6],id[7],id[8]),
                 project = "HGSOC_spatial")
rm(lst)
rm(temp)

## Some small modifications
spatial = JoinLayers(spatial)
spatial$orig.ident = unname(unlist(lapply(strsplit(rownames(spatial@meta.data), split = "_"), function(x){x[[1]]})))
Idents(spatial) = spatial$orig.ident

## Add QC Metrics
spatial[["percent.mito"]] = PercentageFeatureSet(spatial, pattern = "^MT-")
spatial[["percent.ribo"]] = PercentageFeatureSet(spatial, pattern = "^RPL|^RPS|^MRP-")
spatial[["log10GenesPerUMI"]] = log10(spatial$nFeature_RNA)/log10(spatial$nCount_RNA)
VlnPlot(spatial, features = c("nCount_RNA", "nFeature_RNA", "percent.ribo", "percent.mito"), ncol = 1, pt.size = 0)
ggsave("plot/merged_spatial/preprocess/QC.png", width = 4, height = 12, units = "in", dpi = 600)
pie(table(spatial$orig.ident))

## Update Metadata, add Sample Names
spatial$cells = rownames(spatial)
spatial$samples = spatial$orig.ident

## Normalize, no subset because we're using spatial data
spatial = NormalizeData(spatial, normalization.method = "LogNormalize", scale.factor = 10000)
spatial = FindVariableFeatures(spatial, selection.method = "vst", nfeatures = 2000)
p1 = VariableFeaturePlot(spatial)
LabelPoints(plot = p1, points = head(VariableFeatures(spatial), 10), repel = TRUE)

## Add Cell Cycle
regev_lab_cell_cycle_genes = read.delim("data/regev_lab_cell_cycle_genes.txt", header=F) %>% unlist
s.genes = regev_lab_cell_cycle_genes[1:43]
g2m.genes = regev_lab_cell_cycle_genes[44:97]
spatial = CellCycleScoring(object = spatial,
                            s.features = s.genes,
                            g2m.features = g2m.genes,
                            set.ident = TRUE)
cell_cycles = table(spatial$orig.ident, spatial$Phase) %>% as.data.frame.matrix()
cell_cycles = as.data.frame(t(apply(cell_cycles, 1, function(x) x / sum(x) * 100)))
pheatmap::pheatmap(cell_cycles, cluster_cols = F, cluster_rows = F)

## Scale, run PCA, run UMAP, and find Cluster 
spatial = ScaleData(spatial, features = rownames(spatial), 
                    vars.to.regress = c("nCount_RNA", "percent.ribo", "percent.mito",
                                        "S.Score", "G2M.Score"))
spatial = RunPCA(spatial, features = VariableFeatures(object = spatial))

## Run Harmony 
spatial = harmony::RunHarmony(spatial, 
                               group.by.vars = "samples", 
                               reduction = "pca", 
                               assay.use = "RNA", 
                               reduction.save = "harmony", 
                               lambda = 1)
spatial = FindNeighbors(spatial, dims = 1:15, reduction = "harmony")
spatial = FindClusters(object = spatial, resolution = c(0.1, 0.3, 0.5, 0.7), verbose = TRUE)
spatial = RunUMAP(spatial, dims = 1:15, reduction="harmony")

## Save Object 
saveRDS(spatial, "spatial.rds")



