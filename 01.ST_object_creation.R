## SCRIPT: Spatial object creation STutility visium data BM project

## 19.02.23 Laura Sudupe , git @lsudupe

## Libraries
library(Seurat)
library(tidyverse)
library(STutility)


## Define proyect info
ANALYSIS_ID <- "visium_bmn"

## Define project paths
DIR_WD <- getwd()
DIR_ROOT <- file.path(getwd())  
DIR_DATA <- file.path(DIR_ROOT, "data/data/")


## Read data
## samples
samples <- list.files(path = DIR_DATA, recursive = TRUE, pattern = "filtered",full.names = TRUE)
imgs <- list.files(path = DIR_DATA, pattern = "tissue_hires_image.png", recursive = TRUE, full.names = TRUE)
spotfiles <- list.files(path = DIR_DATA, pattern = "tissue_positions_list.csv", recursive = TRUE, full.names = TRUE)
json <- list.files(path = DIR_DATA, pattern = "scalefactors_json.json", recursive = TRUE, full.names = TRUE)
name <- dir(path = DIR_DATA)
#condition <- c("MM", "MM", "MM", "MM", "control", "control", "control", "MM", "MM")


## Create the infotable
infoTable <- data.frame(samples = samples, imgs = imgs, spotfiles = spotfiles, json = json, 
                        name = name)

infoTable <- subset(infoTable, name != "M1_tib_1A" & name != "M3_tib_2A" & grepl("^M", name))

## Create the object
se <- InputFromTable(infoTable,
                     platform =  "Visium")

## Add image info
st.object <- GetStaffli(se)
st.object
se <- LoadImages(se, time.resolve = FALSE)

#se <- ManualAnnotation(se)
#saveRDS(se, "./objects/sp/integrated/se.rds")

## Read the data
se <- readRDS("./objects/sp/integrated/se.rds")

## Image orientation modification
transforms <- list("3" = list("angle" = 90))
se.rotate90 <- WarpImages(se, transforms)
transforms <- list("2" = list("angle" = 90))
se.rotate <- WarpImages(se.rotate90, transforms)
#saveRDS(se.rotate, "./objects/sp/integrated/se.rds")

## Clean bone spots
#se <- readRDS("./objects/sp/integrated/se.rds")
Seurat::Idents(object = se) <- se@meta.data[["labels"]]
se <- SubsetSTData(se, idents = "bone_marrow")

## Normalization step
se <- SCTransform(se, vars.to.regress = "name")
#saveRDS(se_h, "./objects/sp/integrated/se.rds")
#se <- readRDS("./objects/sp/integrated/se.rds")


## QC plots
pdf(file.path("./results/ST/nFeatures_RNA.pdf"))
ST.FeaturePlot(se, features = "nFeature_RNA", cols = c("lightgray", "mistyrose", "red", "dark red", "black"), ncol = 2, pt.size = 1.1)
dev.off()

pdf(file.path("./results/ST/area.pdf"))
FeatureOverlay(se, features = "labels", sampleids = 1:6, ncols = 2,pt.size = 0.7)
dev.off()






