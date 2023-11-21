## SCRIPT: Deconvolution using CARD areas subgroups BM project ST object ONLY HEALTHY

## 28.03.23 Laura Sudupe , git @lsudupe
#https://github.com/YingMa0107/CARD/

#Libraries---------------------------------
library(Seurat)
library(ggplot2)
library(dplyr)
library(CARD)
library(gtools)
library(scatterpie)
library(ggcorrplot)
library(STutility)


## Data
se <- readRDS("./objects/sp/integrated/se.rds")
single_cell_bonemarrow <- readRDS("./objects/heterogeneity/single_cell_bonemarrow.rds")
single_cell_bonemarrow@meta.data[["orig.ident"]] <- "ref"

## Prepare data
#single cell
sub_list <- levels(single_cell_bonemarrow@meta.data[["ident"]])
single_cell_bonemarrow_counts <- single_cell_bonemarrow@assays[["RNA"]]@counts
single_cell_bonemarrow_meta <- single_cell_bonemarrow@meta.data

#spatial
Idents(object = se) <- "name"
name <- unique(se@meta.data[["name"]])
objects <- c()

for (i in name){
  a <- SubsetSTData(se, idents = i)
  ## add object to list
  objects[[length(objects) + 1]] <- a
  
}
names(objects) <- name

## select healthy
objects <- objects[3:4]

## Card
for (i in 1:length(objects)){
  # spatial
  a <- objects[[i]]
  a_count <- a@assays[["SCT"]]@data
  image.info <- a@tools[["Staffli"]]@meta.data
  a_location <- image.info[,1:2]
  
  ###################CARD object creation
  CARD_obj_a = createCARDObject(
    sc_count = single_cell_bonemarrow_counts,
    sc_meta = single_cell_bonemarrow_meta,
    spatial_count = a_count,
    spatial_location = a_location,
    ct.varname = "ident",
    ct.select = unique(single_cell_bonemarrow_meta$ident),
    sample.varname = "metadata....experiment..",
    minCountGene = 100,
    minCountSpot = 5)
  
  ###################Deco
  CARD_a = CARD_deconvolution(CARD_object = CARD_obj_a)
  CARD_a_proportions <- CARD_a@Proportion_CARD
  
  ###################Add proportions to object
  CARD_a
  pro <- as.data.frame(CARD_a@Proportion_CARD)
  cell.types <- colnames(pro)
  for (i in 1:length(cell.types)){
    d <- cell.types[i]
    e <- pro[[d]]
    a@meta.data[[d]] <- e
  }

  ## save object
  saveRDS(a,file = paste0("./objects/card/last/healthy/all/",names(objects[i]),"_subgroup_ST.rds"))
  ## save card results
  saveRDS(CARD_a,file = paste0("./objects/card/last/healthy/all/",names(objects[i]),"_CARD_obj_subgroup_ST.rds"))
}




