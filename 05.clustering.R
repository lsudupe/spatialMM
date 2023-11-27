## SCRIPT: Clustering deconvolution results BM project

## 11.12.22 Laura Sudupe , git @lsudupe

#Libraries---------------------------------
library(Seurat)
library(ggplot2)
library(dplyr)
library(clustertend)
library(factoextra)
library(STutility)
library(corrplot)
library(tidyverse)
library(RColorBrewer)


## Data
se <- readRDS("./objects/sc/integrated/se_deco.rds")

meta <- se@meta.data
types <- meta[,12:18]


## Clustering
hc3 <- eclust(types, k=7, FUNcluster="hclust", hc_metric="euclidean", hc_method = "ward.D2")

## Add clusters to spatial data
a <- as.data.frame(as.factor(hc3[["cluster"]]))
se <- AddMetaData(se, metadata=a, col.name = "clustering")

## Save object
#saveRDS(se, "./objects/heterogeneity/se_hierarchical.rds")

## Plot spatial
b <- SetIdent(se, value = se@meta.data[["clustering"]])
pdf(file.path("./results/endogram/st",filename = "both_spatial_hierarchical.pdf"))
print(FeatureOverlay(b, features = "clustering", sampleids = 1:6, ncols = 2,pt.size = 0.7))
dev.off()

## Plot correlation
matrix <- data.matrix(types, rownames.force = NA)
M <- cor(matrix)

pdf(file.path("./results/endogram/st",filename = "both_cor.pdf"))
print(corrplot(M, method = 'square', title="MM samples cell type correlation",
               col=colorRampPalette(c("blue","white","red"))(100),
               order = 'original', type = 'lower', diag = FALSE,
               mar=c(0,0,1,0)))
dev.off()

## Plot mayor cell type
#coordinates
coor <- se@tools[["Staffli"]]@meta.data
types["clustering"] <- as.vector(se@meta.data[["clustering"]])

df <- types
df["x_coord"] <- as.vector(coor$x) 
df["y_coord"] <- as.vector(coor$y) 
df["sample"] <- as.vector(meta$name) 

#convert factor columns to numeric
df <- df %>%
  mutate(across(where(is.factor), as.numeric))

#filter values below the 20% threshold for each row
df_filtered <- df %>%
  mutate(across(-c(clustering, x_coord, y_coord, sample), ~ifelse(. < 0.20, 0, .)))

#rescale non-zero values in each row to sum up to 1 and overwrite the original values
df_rescaled <- df_filtered %>%
  rowwise() %>%
  mutate(across(-c(clustering, x_coord, y_coord, sample), 
                ~ifelse(. != 0, . / sum(c_across(-c(clustering, x_coord, y_coord, sample))[c_across(-c(clustering, x_coord, y_coord, sample)) != 0]), 0))) %>%
  ungroup()

#reshape the data frame to a long format for plotting
df_new <- df_rescaled
df_new$x_coord <- NULL
df_new$y_coord <- NULL
df_new$sample <- NULL

df_long <- df_new %>%
  pivot_longer(cols = -clustering,
               names_to = "cell_type",
               values_to = "percentage") %>%
  filter(percentage > 0)

#sum percentage by cell type and cluster
df_aggregated <- df_long %>%
  group_by(clustering, cell_type) %>%
  summarise(total_percentage = sum(percentage)) %>%
  ungroup()

#normalize the percentaje for the sum to be 1
df_normalized <- df_aggregated %>%
  group_by(clustering) %>%
  mutate(normalized_percentage = total_percentage / sum(total_percentage)) %>%
  ungroup()

#color
cells_order <- c("Bcell", "DC", "Erythroblasts","MM_MIC", "Monocytes","Neutrophils","Tcell")
cell_type_colors <- brewer.pal(length(cells_order), "RdBu")
cell_type_colors <-  c("#ef8a62" ,"#ffd1b6" ,"#faeae0" ,"#b2182b" ,"#d1e5f0" ,"#67a9cf" ,"#2166ac")
cell_type_color_map <- setNames(cell_type_colors, cells_order)

#plot
ggplot(df_normalized, aes(x = clustering, y = cell_type, size = normalized_percentage, color = cell_type)) +
  geom_point(alpha = 1) +
  # El resto de tu código de ggplot
  labs(title = "Major Cell Types Driving Each Cluster",
       x = "Cluster",
       y = "Cell Type",
       size = "Percentage") +
  theme_minimal() +
  scale_y_discrete(limits = rev(cells_order)) +
  scale_color_manual(name = "Cell Type", values = cell_type_color_map, limits = cells_order) +
  scale_size_continuous(range = c(1, 5)) +  # Ajusta esto según sea necesario
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major = element_line(color = "grey80", size = 0.05),
    panel.grid.minor = element_blank()
  ) 

## Plot piechart
# Separate the dataframe by sample
df_list <- split(df_rescaled, df_rescaled$sample)

for (i in 1:length(df_list)){
  data <- df_list[[i]]
  #reshape the data frame to a long format for plotting
  data$sample <- NULL
  data$clustering <- NULL
  colnames(data)[which(names(data) == "x_coord")] <- "x"
  colnames(data)[which(names(data) == "y_coord")] <- "y"
  
  #take out coordinates
  num_cell_types <- ncol(data) - 2
  
  data_long <- data %>% 
    tidyr::pivot_longer(cols = 1:num_cell_types, names_to = "cell_type", values_to = "percentage")
  
  data_polar <- data_long %>%
    group_by(x, y) %>%
    mutate(end = cumsum(2 * pi * percentage),
           start = lag(end, default = 0)) %>%
    ungroup()
  
  #filter out values equal to or below 0
  data_polar <- data_polar %>%
    filter(percentage > 0)
  #spot size
  num_points <- 100
  radius <- 0.73 
  
  #create pie
  data_pie <- tidyr::expand_grid(x = unique(data_polar$x),
                                 y = unique(data_polar$y),
                                 point = seq_len(num_points + 1)) %>%
    left_join(data_polar, by = c("x", "y")) %>%
    mutate(angle = ifelse(point == 1, 0, start + (end - start) * (point - 2) / (num_points - 1)),
           x_plot = x + radius * cos(angle) * (point != 1),
           y_plot = y + radius * sin(angle) * (point != 1))
  
  #cell types and colors
  cells_order <- c("Bcell", "DC", "Erythroblasts","MM_MIC", "Monocytes","Neutrophils","Tcell")
  cell_type_colors <-  c("#ef8a62" ,"#ffd1b6" ,"#faeae0" ,"#b2182b" ,"#d1e5f0" ,"#67a9cf" ,"#2166ac")
  cell_type_color_map <- setNames(cell_type_colors, cells_order)
  
  pie_chart_plot <- ggplot(data_pie, aes(x = x_plot, y = y_plot, group = interaction(x, y, cell_type))) +
    geom_polygon(aes(fill = cell_type), color = "white") +
    coord_fixed() +
    scale_y_reverse() +
    theme_void() +
    theme(legend.position = "right") +
    labs(fill = "Cell Type") +
    scale_fill_manual(values = cell_type_color_map)
  
  pdf(paste("./results/piechart/", names(df_list[i]),"_piechart_20percent.pdf",sep=""))
  print(pie_chart_plot)
  dev.off()
  
}






