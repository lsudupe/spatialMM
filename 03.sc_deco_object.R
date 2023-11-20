## SCRIPT: Singel cell reference data preparing BM project

## 23.11.22 Laura Sudupe , git @lsudupe

#Libraries---------------------------------
library(Seurat)
library(scCustomize)
library(dplyr)
library(tidyverse)
library(base)
library(harmony)

## Data -- Load data from: https://apps.embl.de/nicheview/
load(file='./data/single-cell/RNAMagnetDataBundle/NicheData10x.rda')
single_cell_bonemarrow = UpdateSeuratObject(object = NicheData10x)
single_cell_bonemarrow@meta.data[["cell"]] <- single_cell_bonemarrow@active.ident

## Group data
a <- as.factor(single_cell_bonemarrow@meta.data[["cell"]])
levels(a) <- list(PSC  = "Ery/Mk prog.", 
                  PSC = "Neutro prog.",
                  PSC = "Mono prog.",
                  PSC = "Gran/Mono prog.",
                  PSC = "LMPPs",
                  Bcell = "large pre-B.",
                  PSC = "Mk prog.",
                  Erythroblasts = "Erythroblasts",
                  PSC = "Eo/Baso prog.",
                  Monocytes = "Monocytes",
                  PSC = "Ery prog.",
                  Bcell = "pro-B",
                  Tcell = "T cells",
                  Neutrophils = "Neutrophils",
                  MSC = "Adipo-CAR",
                  MSC = "Ng2+ MSCs",
                  NC = "Schwann cells",
                  MSC = "Osteoblasts",
                  Fibro = "Arteriolar fibro.",
                  EC = "Sinusoidal ECs",
                  MSC = "Osteo-CAR",
                  Bcell = "small pre-B.",
                  MSC = "Chondrocytes",
                  Fibro = "Endosteal fibro.",
                  Fibro = "Fibro/Chondro p.",
                  Fibro = "Stromal fibro.",
                  EC = "Arteriolar ECs",
                  Fibro = "Myofibroblasts",
                  MSC = "Smooth muscle",
                  DC = "Dendritic cells",
                  NK = "NK cells",
                  Bcell = "B cell"
)

#SMC, NC, Fibro_Chondro_p take out, 
single_cell_bonemarrow@meta.data[["ident"]] <- a
Seurat::Idents(object = single_cell_bonemarrow) <- single_cell_bonemarrow@meta.data[["ident"]]

single_cell_bonemarrow <- subset(x = single_cell_bonemarrow, idents = c("NC", "PSC", "Fibro", "MSC", "EC", "NK"), invert = TRUE)
single_cell_bonemarrow@meta.data[["ident"]] <- single_cell_bonemarrow@active.ident

## Save healthy bone marrow object
#saveRDS(single_cell_bonemarrow, "./objects/heterogeneity/single_cell_bonemarrow.rds")

## Add in HOUSE plasma cell data
PC_MM <- readRDS("./data/single-cell/PC/scRNA_MM_PC.rds")
PC_MM <- subset(x = PC_MM, idents = c("MM_MIC"))
PC_MM@meta.data[["ident"]] <- PC_MM@active.ident

## Merge data
single_cell_bonemarrow@meta.data[["orig.ident"]] <- "ref"
combined_sc <- merge(single_cell_bonemarrow, y = c( PC_MM), 
                     add.cell.ids = c("single_cell_bonemarrow", "PC_MM"), project = "BM")
combined_sc@meta.data[["split"]] <- combined_sc@active.ident
Seurat::Idents(object = combined_sc) <- combined_sc@meta.data[["orig.ident"]]
x <- combined_sc

## Integrate using harmony
harmony_i <- harmony %>%
  NormalizeData() %>%
  ScaleData() %>%
  FindVariableFeatures() %>%
  RunPCA(npcs = 20) %>%
  RunHarmony(assay.use="RNA",reduction = "pca", dims = 1:20, group.by.vars = "orig.ident") %>%
  FindNeighbors(reduction = "harmony", dims = 1:20) %>%
  FindClusters(resolution=0.7) %>%
  RunUMAP(reduction = "harmony", dims = 1:20, n.epochs = 1e3) 

## Save multiple myeloma plasma cells
saveRDS(harmony_i, "./objects/sc/integrated/integrated_sc_harmony.rds")




