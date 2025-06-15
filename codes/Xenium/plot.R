library(dplyr)
library(magrittr)
library(Seurat)
library(Matrix)
library(ggplot2)
library(EBImage)
library(scales)
library(future.apply)
library(progressr)

set.seed(123)
options(future.globals.maxSize = 1e+09) # Set to 1 GB

# setwd("/home/polly_hung/spatial-seq")

## Processing ------------------------------------------------------------------
# seurat <- readRDS("ST_Discovery_so.rds")
meta <- read.csv("~/Desktop/test.csv", row.names = 1)
sample <- unique(meta$samples)[[1]]

plotEBI <- function(sample){
  
  ## Read in Files 
  meta_sub <- meta %>% dplyr::filter(samples == sample)
  seg_path <- paste0("~/Desktop/Collection/spatial-seq/10xXenium/Test/segments/_whole-cell_03.csv")
  
  ## Build empty EBI
  cellseg = read.csv(seg_path)
  colnames(cellseg) <- unlist(lapply(colnames(cellseg), function(x) {
    strsplit(x, split = "X0.")[[1]][2]
  }))
  colnames(cellseg)[1] <- "0"
  cellmask <- cellseg
  cellmask[cellmask != 0] <- 0
  if("black" == "white") {
    cellmask[cellmask == 0] <- 1
  } else {
    cellmask[cellmask == 0] <- 0
  }
  cellmask <- EBImage::Image(as.matrix(cellmask))
  cellmask <- EBImage::channel(cellmask, "rgb")
  
  ## Add Cell Colours 
  plan(multisession) 
  cell_index_map <- setNames(as.integer(gsub(".*_c", "", rownames(meta_sub))), rownames(meta_sub))
  rgb_matrix <- t(col2rgb(meta_sub$YAP1_colour)) / 255
  rownames(rgb_matrix) <- rownames(meta_sub)
  
  ## use subsets 
  subset1 <- rownames(meta_sub)[1:100]
  
  with_progress({
    p <- progressor(along = subset1)
    future_lapply(subset1, function(cell_id) {
      p() 
      cell_num <- cell_index_map[[cell_id]]
      rgb_vals <- rgb_matrix[cell_id, ]
      cell_mask <- cellseg == cell_num
      cellmask@.Data[,,1][cell_mask] <- rgb_vals[1]
      cellmask@.Data[,,2][cell_mask] <- rgb_vals[2]
      cellmask@.Data[,,3][cell_mask] <- rgb_vals[3]
      NULL
    }, future.seed = TRUE)
  })
  
  EBImage::writeImage(cellmask, paste0("~/Desktop/", sample, ".png"))
}


cell_2_rgb[["cancer1"]] = c(unname(col2rgb("#344CB7")[,1]))
cell_2_rgb[["cancer2"]] = c(unname(col2rgb("#6173C7")[,1]))
cell_2_rgb[["cancer3"]] = c(unname(col2rgb("#8E9BD7")[,1]))
cell_2_rgb[["cancer4"]] = c(unname(col2rgb("#BBC3E7")[,1]))
cell_2_rgb[["cancer5"]] = c(unname(col2rgb("#E8EBF7")[,1]))
cell_2_rgb[["cancer6"]] = c(unname(col2rgb("#FCE7E7")[,1]))
cell_2_rgb[["cancer7"]] = c(unname(col2rgb("#F7B7B7")[,1]))
cell_2_rgb[["cancer8"]] = c(unname(col2rgb("#F18888")[,1]))
cell_2_rgb[["cancer9"]] = c(unname(col2rgb("#EC5858")[,1]))
cell_2_rgb[["cancer10"]] = c(unname(col2rgb("#E72929")[,1]))