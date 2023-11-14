## SCRIPT: QC of the combinedurat spatial objects BMN project

## 24.06.22 Laura Sudupe , git @lsudupe

#Libraries---------------------------------
library(Seurat)
library(dplyr)
library(ggplot2)
library(BuenColors)

## Read data
combined <- readRDS("./objects/sp/integrated/combined.rds")

##Visualization
# Visualize the number of spots counts per sample
pdf(file.path("./results/QC",filename = "Number of spot per sample.pdf"))
combined@meta.data%>% 
  ggplot(aes(x=name, fill=name)) + 
  geom_bar() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("NSpots")
dev.off()

# Visualize the distribution of genes detected per spot via boxplot
pdf(file.path("./results/QC",filename = "genes detected per spot boxplot.pdf"))
combined@meta.data %>% 
  ggplot(aes(x=name, y=log10(nFeature_RNA), fill=name)) + 
  geom_boxplot() + 
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("Ngenes vs Npots")
dev.off()

# Visualize the distribution of genes detected per spot via histogram
pdf(file.path("./results/QC",filename = "genes detected per spot histogram.pdf"))
combined@meta.data %>% 
  ggplot(aes(color=name, x=nFeature_RNA, fill=name)) + 
  geom_density(alpha = 0.2) + 
  theme_classic() +
  scale_x_log10() + 
  geom_vline(xintercept = 300)
dev.off()

# Visualize the number UMIs/transcripts per cell
pdf(file.path("./results/QC",filename = "number UMIs sati transcripts per spot.pdf"))
combined@meta.data %>% 
  ggplot(aes(color=name, x=nCount_RNA, fill=name)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  ylab("Cell density") +
  geom_vline(xintercept = 500)
dev.off()


## Spatial plots
color <- jdb_palette("brewer_spectra")

feature.list <- c("nCount_RNA", "nFeature_RNA")

## Divide by sample
Idents(object = combined) <- "name"
name <- unique(combined@meta.data[["name"]])
objects <- c()

for (i in name){
  a <- SubsetSTData(combined, idents = i)
  ## add object to list
  objects[[length(objects) + 1]] <- a
  
}
names(objects) <- name
## combined separate list
list2env(objects,envir=.GlobalEnv)

for (i in 1:length(objects)){
  a <- objects[[i]]
  b <- c(min(a@meta.data[["nFeature_RNA"]]), max(a@meta.data[["nFeature_RNA"]]))
  p1 <- FeatureOverlay(a, features = c("nFeature_RNA"), ncols = 1, pt.size = 1.5, value.scale = "all") +
        scale_fill_gradientn(colours = color,
                        breaks = b,
                        limits = b)
  pdf(paste("./results/QC/",names(objects[i]),"_feature_spatial.pdf",sep=""))
  print(p1)
  dev.off()
}

for (i in 1:length(objects)){
  a <- objects[[i]]
  b <- c(min(a@meta.data[["nCount_RNA"]]), max(a@meta.data[["nCount_RNA"]]))
  p1 <- FeatureOverlay(a, features = c("nCount_RNA"), ncols = 1, pt.size = 1.5, value.scale = "all") +
    scale_fill_gradientn(colours = color,
                         breaks = b,
                         limits = b)
  pdf(paste("./results/QC/",names(objects[i]),"_count_spatial.pdf",sep=""))
  print(p1)
  dev.off()
}

combined.meta <- combined@meta.data


pdf(paste("./results/QC/features_violinplot.pdf",sep=""))
VlnPlot(object = combined, features = 'nFeature_RNA', split.by = 'name')
dev.off()

pdf(paste("./results/QC/count_violinplot.pdf",sep=""))
VlnPlot(object = combined, features = 'nCount_RNA', split.by = 'name')
dev.off()



## Filter 
filtered_combined <- SubsetSTData(combined, expression =
                            (nCount_RNA >= 200 & nCount_RNA <= 100000) & 
                            (nFeature_RNA >= 150 & nFeature_RNA <= 10000))

## Save combined
saveRDS(filtered_combined,"./objects/sp/combined_filtered.rds")
