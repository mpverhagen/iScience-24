library(Seurat)

# plotting and data science packages
library(tidyverse)
library(cowplot)
library(patchwork)

# install this package, which allows us to compute distance between the spots
library(proxy)

# set random seed for reproducibility
set.seed(12345)

path <- "/mnt/data/Projects/SW480_spheres/Out/Seurat/"

#### Neighbourhood Analysis ####
library(ggplot2)
library(dplyr)
library(sf)


Spatial_CRC_integrated_metadata <- read.csv("/mnt/data/Projects/SW480_spheres/Out/Spatial_CRC_integrated_metadata.csv", row.names=1)
Palette <- c("#FFD4D4", "#CFB997", "#E3F6FF", "#8DCBE6", "#EF9A53", "#E3ACF9", "#EAE0DA", "lightgrey", "#ABC270", "#658864", "#CD0404")
Palette_2 <- c("#FFD4D4", "#CFB997",  "#8DCBE6", "lightgrey", "#ABC270",  "#CD0404")


# Spatial plot example
Spatial_CRC_integrated_metadata %>% filter(Sample == "CRC1") %>%
  ggplot(aes(x = col, y = row, color = Annotation)) + geom_point(size = 1) + facet_wrap(~ Sample, ncol = 5, scales = "free") +
  scale_color_manual(values = Palette_2) + theme_classic()


Sample.take <- names(table(Spatial_CRC_integrated_metadata$Sample))
Spatial_out <- data.frame()

for (i in 1:length(Sample.take)){
  print(Sample.take[i])
  
  SF_dat <- st_as_sf(Spatial_CRC_integrated_metadata %>% dplyr::filter(Sample == Sample.take[i]), coords = c("row","col"))
  
  # Make a honey comb grids  
  area_honeycomb_grid = st_make_grid(SF_dat, c(4, 4), what = "polygons", square = FALSE)
  
  # To sf and add grid ID
  honeycomb_grid_sf = st_sf(area_honeycomb_grid) %>%
    mutate(grid_id = paste(Sample.take[i], 1:length(lengths(area_honeycomb_grid)), sep = "-"))   # add grid ID
  
  Spatial_temp <- SF_dat %>% 
    st_join(honeycomb_grid_sf) %>%
    filter(!duplicated(Row.names))
  
  Spatial_out <- rbind(Spatial_out, Spatial_temp)
  rm(Spatial_temp)
}

# Plot Neighboorhoods
Spatial_out %>% filter(Sample == "CRC1") %>%
  ggplot(aes(x = imagecol, y = imagerow, color = Annotation)) + geom_point(size = 1) + facet_wrap(~ Sample, ncol = 5, scales = "free") +
  scale_color_manual(values = Palette_2) + theme_classic() + 
  Spatial_out %>% filter(Sample == "CRC1") %>%
  ggplot(aes(x = imagecol, y = imagerow, color = grid_id)) + geom_point(size = 1) + facet_wrap(~ Sample, ncol = 5, scales = "free") + theme_classic() + 
  Seurat::NoLegend() + scale_color_manual(values = rep(c("orange", "blue", "red", "brown","grey", "black", "magenta", "pink", "darkgreen","white", "yellow", "green",  "darkblue"), times = 1000))


### Neighbourhood Seurat ####
library(Seurat)
library(ggplot2)
library(dplyr)
library(reshape2)
setwd("/mnt/data/Projects/SW480_spheres/")
Palette <- c("#FFD4D4", "#CFB997", "#E3F6FF", "#8DCBE6", "#EF9A53", "#E3ACF9", "#EAE0DA", "lightgrey", "#ABC270", "#658864", "#CD0404")
Palette_2 <- c("#FFD4D4", "#CFB997",  "#8DCBE6", "lightgrey", "#ABC270",  "#CD0404")

Spatial_CRC_integrated <- readRDS("/mnt/data/Projects/SW480_spheres/Out/Seurat/Spatial_merge_integrated.rds")

Spatial_CRC_integrated <- AddMetaData(Spatial_CRC_integrated, metadata = data.frame(Spatial_out[,21:22]))


# Merge  data by Neighbourhood
#NN_dat <- AverageExpression(Spatial_CRC_integrated, group.by = "grid_id", assays = "integrated", slot = "data")

NN_dat <- AggregateExpression(Spatial_CRC_integrated, group.by = "grid_id", assays = "Spatial", slot = "counts")
NN_dat <- NN_dat$Spatial

# Only keep genes with at least 10 counts in total data sets
NN_dat <- NN_dat[rowSums(NN_dat) > 10,]

NN_metadat <- data.frame(Spatial_CRC_integrated@meta.data %>% group_by(grid_id, Sample) %>% summarize(Nspots = n()))
NN_metadat <- left_join(NN_metadat,Spatial_CRC_integrated@meta.data %>% group_by(grid_id, Sample, Celltype) %>% dcast(grid_id ~ Celltype ))
rownames(NN_metadat) <- NN_metadat$grid_id

# Create Seurat Object
NN.Seur <- CreateSeuratObject(counts = NN_dat,
                              meta.data = NN_metadat)


NN.Seur <- NN.Seur %>% NormalizeData() %>%  FindVariableFeatures() %>% ScaleData() %>% RunPCA() %>% RunUMAP(1:30)


DimPlot(NN.Seur, group.by = "Sample")
DimPlot(NN.Seur, group.by = "Nspots")

saveRDS(NN.Seur, "/mnt/data/Projects/SW480_spheres/Out/Seurat/Spatial_Neighbourhood.rds")

### Neighbourhood Rpca ####
library(Seurat)
library(ggplot2)
library(dplyr)
setwd("/mnt/data/Projects/SW480_spheres/")
Palette <- c("#FFD4D4", "#CFB997", "#E3F6FF", "#8DCBE6", "#EF9A53", "#E3ACF9", "#EAE0DA", "lightgrey", "#ABC270", "#658864", "#CD0404")
Palette_2 <- c("#FFD4D4", "#CFB997",  "#8DCBE6", "lightgrey", "#ABC270",  "#CD0404")


NN.Seur <- readRDS("/mnt/data/Projects/SW480_spheres/Out/Seurat/Spatial_Neighbourhood.rds")

DimPlot(NN.Seur, reduction = "umap", group.by = "Sample", label = TRUE, repel = TRUE) + NoLegend()


# split the dataset into a list of two seurat objects (stim and CTRL)
CRC.list <- SplitObject(NN.Seur, split.by = "Sample")

# normalize and identify variable features for each dataset independently
CRC.list <- lapply(X = CRC.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

# select features that are repeatedly variable across datasets for integration run PCA on each
# dataset using these features
features <- SelectIntegrationFeatures(object.list = CRC.list)
CRC.list <- lapply(X = CRC.list, FUN = function(x) {
  x <- ScaleData(x, features = features, verbose = TRUE)
  x <- RunPCA(x, features = features, verbose = FALSE, approx=FALSE)
})

CRC.anchors <- FindIntegrationAnchors(object.list = CRC.list, anchor.features = features, reduction = "rpca")

# this command creates an 'integrated' data assay
NN.Seur_integrated <- IntegrateData(anchorset = CRC.anchors, k.weight = 25, verbose = TRUE)

# preprocess RNA assay
DefaultAssay(NN.Seur_integrated) <- "RNA"
NN.Seur_integrated <- NN.Seur_integrated %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData()

# specify that we will perform downstream analysis on the corrected data note that the
# original unmodified data still resides in the 'RNA' assay
DefaultAssay(NN.Seur_integrated) <- "integrated"

# Run the standard workflow for visualization and clustering
NN.Seur_integrated <- ScaleData(NN.Seur_integrated, verbose = FALSE)
NN.Seur_integrated <- RunPCA(NN.Seur_integrated, npcs = 50, verbose = FALSE)
NN.Seur_integrated <- FindNeighbors(NN.Seur_integrated, reduction = "pca", dims = 1:30)
NN.Seur_integrated <- RunUMAP(NN.Seur_integrated, reduction = "pca", dims = 1:50, min.dist = 0.2, n.neighbors = 100, spread = 2)
NN.Seur_integrated <- FindClusters(NN.Seur_integrated, resolution = 0.3)

DimPlot(NN.Seur_integrated, reduction = "umap", group.by = "seurat_clusters", label = TRUE, repel = TRUE)
DimPlot(NN.Seur_integrated, reduction = "umap", group.by = "Sample")

DefaultAssay(NN.Seur_integrated) <- "RNA"

saveRDS(NN.Seur_integrated, "/mnt/data/Projects/SW480_spheres/Out/Seurat/Spatial_Neighbourhood_integrated.rds")

### Neighbourhood evaluation ####
library(Seurat)
library(ggplot2)
library(dplyr)
library(reshape2)
setwd("/mnt/data/Projects/SW480_spheres/")
Palette <- c("#FFD4D4", "#CFB997", "#E3F6FF", "#8DCBE6", "#EF9A53", "#E3ACF9", "#EAE0DA", "lightgrey", "#ABC270", "#658864", "#CD0404")
Palette_2 <- c("#FFD4D4", "#CFB997",  "#8DCBE6", "lightgrey", "#ABC270",  "#CD0404")

NN.Seur_integrated <- readRDS("/mnt/data/Projects/SW480_spheres/Out/Seurat/Spatial_Neighbourhood_integrated.rds")

DimPlot(NN.Seur_integrated, reduction = "umap", group.by = "seurat_clusters", label = TRUE, repel = TRUE)

FeaturePlot(NN.Seur_integrated, features = c("EPCAM", "VIM"))

Pct.dat <- NN.Seur_integrated@meta.data %>% 
  dplyr::select(-nFeature_RNA, -nCount_RNA, -grid_id, -Sample, -Nspots, -orig.ident, -integrated_snn_res.0.3, -integrated_snn_res.0.2, -Niche) %>%
  reshape2::melt(id.vars = "seurat_clusters", value.name = "Nspots") %>%
  group_by(variable, seurat_clusters) %>% summarise(N = sum(Nspots)) %>%
  mutate(pct = round(N/sum(N)*100, 2))

Pct.dat.heat <- Pct.dat %>% dplyr::select(-N) %>% dcast(variable ~ seurat_clusters)
rownames(Pct.dat.heat) <- Pct.dat.heat$variable
Pct.dat.heat <- Pct.dat.heat[,-1]

pheatmap::pheatmap(Pct.dat.heat, angle_col = 0, display_numbers = T, treeheight_row = 15, treeheight_col = 15,
                   border_color = "black", number_color = "white",
                   color  = viridis::inferno(10)[2:8],
                   cluster_cols = F)


# Annotation
NN.Seur_integrated$Niche <- "Other"
NN.Seur_integrated$Niche[NN.Seur_integrated$seurat_clusters == "0"] <- "Liver niche"
NN.Seur_integrated$Niche[NN.Seur_integrated$seurat_clusters == "1"] <- "Tumor core"
NN.Seur_integrated$Niche[NN.Seur_integrated$seurat_clusters == "2"] <- "Tumor front"
NN.Seur_integrated$Niche[NN.Seur_integrated$seurat_clusters == "3"] <- "Stroma"
NN.Seur_integrated$Niche[NN.Seur_integrated$seurat_clusters == "4"] <- "Immune niche"
NN.Seur_integrated$Niche[NN.Seur_integrated$seurat_clusters == "5"] <- "Colon epithelium"
NN.Seur_integrated$Niche[NN.Seur_integrated$seurat_clusters == "6"] <- "Endothelial niche"
NN.Seur_integrated$Niche[NN.Seur_integrated$seurat_clusters == "7"] <- "Colon epithelium"
NN.Seur_integrated$Niche[NN.Seur_integrated$seurat_clusters == "8"] <- "Liver epithelium"
NN.Seur_integrated$Niche[NN.Seur_integrated$seurat_clusters == "9"] <- "Smooth muscle"
NN.Seur_integrated$Niche[NN.Seur_integrated$seurat_clusters == "10"] <- "ENS"


Pct.dat <- NN.Seur_integrated@meta.data %>% 
  dplyr::select(-nFeature_RNA, -nCount_RNA, -grid_id, -Sample, -Nspots, -orig.ident, -integrated_snn_res.0.3, -integrated_snn_res.0.2, -seurat_clusters) %>%
  reshape2::melt(id.vars = "Niche", value.name = "Nspots") %>%
  group_by(variable, Niche) %>% summarise(N = sum(Nspots)) %>%
  mutate(pct = round(N/sum(N)*100, 2))

Pct.dat.heat <- Pct.dat %>% dplyr::select(-N) %>% dcast(variable ~ Niche)
rownames(Pct.dat.heat) <- Pct.dat.heat$variable
Pct.dat.heat <- Pct.dat.heat[,-1]

pheatmap::pheatmap(Pct.dat.heat, angle_col = 90, display_numbers = T, treeheight_row = 15, treeheight_col = 15,
                   border_color = "black", number_color = "white",
                   color  = viridis::inferno(10)[2:8],
                   cluster_cols = F)

Palette_niche <- c("lightgrey", "#00B2CA", "#F865B0", "#1D4E89", "#C5AFA0", "#E9BCB7", "#285943", "#94A187","#984447", "#DF2935")
DimPlot(NN.Seur_integrated, group.by = "Niche", label = F, cols = Palette_niche)

saveRDS(NN.Seur_integrated, "/mnt/data/Projects/SW480_spheres/Out/Seurat/Spatial_Neighbourhood_integrated.rds")
write.csv(NN.Seur_integrated@meta.data, "/mnt/data/Projects/SW480_spheres/Out/Spatial_Neighbourhood_Annotated_meta.data.csv")

### Back to Spatial Spot Seurat ####
library(Seurat)
library(ggplot2)
library(dplyr)
library(reshape2)
library(sf)


Spatial_CRC_integrated_metadata <- read.csv("/mnt/data/Projects/SW480_spheres/Out/Spatial_CRC_integrated_metadata.csv", row.names=1)
Palette <- c("#FFD4D4", "#CFB997", "#E3F6FF", "#8DCBE6", "#EF9A53", "#E3ACF9", "#EAE0DA", "lightgrey", "#ABC270", "#658864", "#CD0404")
Palette_2 <- c("#FFD4D4", "#CFB997",  "#8DCBE6", "lightgrey", "#ABC270",  "#CD0404")
Palette_niche <- c("lightgrey", "#00B2CA", "#F865B0", "#1D4E89", "#C5AFA0", "#E9BCB7", "#285943", "#94A187","#984447", "#DF2935")


Sample.take <- names(table(Spatial_CRC_integrated_metadata$Sample))
Spatial_out <- data.frame()

for (i in 1:length(Sample.take)){
  print(Sample.take[i])
  
  SF_dat <- st_as_sf(Spatial_CRC_integrated_metadata %>% dplyr::filter(Sample == Sample.take[i]), coords = c("row","col"))
  
  # Make a honey comb grids  
  area_honeycomb_grid = st_make_grid(SF_dat, c(4, 4), what = "polygons", square = FALSE)
  
  # To sf and add grid ID
  honeycomb_grid_sf = st_sf(area_honeycomb_grid) %>%
    mutate(grid_id = paste(Sample.take[i], 1:length(lengths(area_honeycomb_grid)), sep = "-"))   # add grid ID
  
  Spatial_temp <- SF_dat %>% 
    st_join(honeycomb_grid_sf) %>%
    filter(!duplicated(Row.names))
  
  Spatial_out <- rbind(Spatial_out, Spatial_temp)
  rm(Spatial_temp)
}

# Plot Neighboorhoods
Spatial_out %>% filter(Sample == "CRC1") %>%
  ggplot(aes(x = imagecol, y = imagerow, color = Annotation)) + geom_point(size = 1) + facet_wrap(~ Sample, ncol = 5, scales = "free") +
  scale_color_manual(values = Palette_2) + theme_classic() + 
  Spatial_out %>% filter(Sample == "CRC1") %>%
  ggplot(aes(x = imagecol, y = imagerow, color = grid_id)) + geom_point(size = 1) + facet_wrap(~ Sample, ncol = 5, scales = "free") + theme_classic() + 
  Seurat::NoLegend() + scale_color_manual(values = rep(c("orange", "blue", "red", "brown","grey", "black", "magenta", "pink", "darkgreen","white", "yellow", "green",  "darkblue"), times = 1000))


Spatial_CRC_integrated <- readRDS("/mnt/data/Projects/SW480_spheres/Out/Seurat/Spatial_merge_integrated.rds")
# Add to Seurat object
Spatial_CRC_integrated <- AddMetaData(Spatial_CRC_integrated, metadata = data.frame(Spatial_out[,21:22]))

Spatial_Neighbourhood_Annotated <- read.csv( "/mnt/data/Projects/SW480_spheres/Out/Spatial_Neighbourhood_Annotated_meta.data.csv", row.names = 1)
colnames(Spatial_Neighbourhood_Annotated)[21] <- "Niche_integrated_snn_res.0.3"
Spatial_Neighbourhood_Annotated <- Spatial_Neighbourhood_Annotated[, c("grid_id","Niche_integrated_snn_res.0.3", "Niche")]

new_meta.data <- left_join(Spatial_CRC_integrated@meta.data, Spatial_Neighbourhood_Annotated, by = "grid_id")
rownames(new_meta.data) <- new_meta.data$Row.names
new_meta.data <- new_meta.data[,24:25]

# Add to Seurat object
Spatial_CRC_integrated <- AddMetaData(Spatial_CRC_integrated, metadata = data.frame(new_meta.data))




Spatial_CRC_integrated@meta.data %>% filter(Sample == "CRC1") %>%
  ggplot(aes(x = imagecol, y = imagerow, color = Annotation)) + geom_point(size = 1) + facet_wrap(~ Sample, ncol = 5, scales = "free") +
  scale_color_manual(values = Palette_2) + theme_classic() + 
  Spatial_CRC_integrated@meta.data %>% filter(Sample == "CRC1") %>%
  ggplot(aes(x = imagecol, y = imagerow, color = grid_id)) + geom_point(size = 1) + facet_wrap(~ Sample, ncol = 5, scales = "free") + theme_classic() + 
  Seurat::NoLegend() + scale_color_manual(values = rep(c("orange", "blue", "red", "brown","grey", "black", "magenta", "pink", "darkgreen","white", "yellow", "green",  "darkblue"), times = 1000)) +
  Spatial_CRC_integrated@meta.data %>% filter(Sample == "CRC1") %>%
  ggplot(aes(x = imagecol, y = imagerow, color = Niche)) + geom_point(size = 1) + facet_wrap(~ Sample, ncol = 5, scales = "free") + theme_classic() + 
  scale_color_manual(values = Palette_niche)

Spatial_CRC_integrated@meta.data %>%
  ggplot(aes(x = imagecol, y = imagerow, color = Niche)) + geom_point(size = 1) + facet_wrap(~ Sample, ncol = 5, scales = "free") + theme_classic() + 
  scale_color_manual(values = Palette_niche) + theme_void()

saveRDS(Spatial_CRC_integrated, "/mnt/data/Projects/SW480_spheres/Out/Seurat/Spatial_merge_integrated_with_niches.rds")

### START: Load spatial niche: subcluster tumor####
library(Seurat)
library(ggplot2)
library(dplyr)
library(reshape2)
#library(sf)
Palette <- c("#FFD4D4", "#CFB997", "#E3F6FF", "#8DCBE6", "#EF9A53", "#E3ACF9", "#EAE0DA", "lightgrey", "#ABC270", "#658864", "#CD0404")
Palette_sub <- c("#FFD4D4", "#CFB997", "#E3F6FF", "#8DCBE6", "#EF9A53", "#E3ACF9", "#EAE0DA", "lightgrey", "#ABC270", "#658864", "#1B998B", "firebrick", "orange", "#93827F","#BCD8B7",   "blue", "magenta", "pink")



Palette_2 <- c("#FFD4D4", "#CFB997",  "#8DCBE6", "lightgrey", "#ABC270",  "#CD0404")
Palette_niche <- c("lightgrey", "#00B2CA", "#F865B0", "#1D4E89", "#C5AFA0", "#E9BCB7", "#285943", "#94A187","#984447", "#DF2935")

Spatial_CRC_integrated <- readRDS("/mnt/data/Projects/SW480_spheres/Out/Seurat/Spatial_merge_integrated_with_niches.rds")

DimPlot(Spatial_CRC_integrated, group.by = "Celltype", cols = Palette)

DefaultAssay(Spatial_CRC_integrated) <- "integrated"
Spatial_CRC_integrated <- SetIdent(Spatial_CRC_integrated, value = "Celltype")

Spatial_CRC_integrated <- FindSubCluster(
  Spatial_CRC_integrated,
  "Tumor",
  graph.name = "integrated_nn",
  subcluster.name = "sub.cluster",
  resolution = 0.25,
  algorithm = 1
)
DefaultAssay(Spatial_CRC_integrated) <- "Spatial"

table(Spatial_CRC_integrated$sub.cluster, Spatial_CRC_integrated$Celltype)

DimPlot(Spatial_CRC_integrated, group.by = "sub.cluster", cols = Palette_sub)

Cluster_cor <- AverageExpression(Spatial_CRC_integrated, assays = "integrated", group.by = "sub.cluster")
#pheatmap::pheatmap(cor(Cluster_cor$integrated), number_color = "black", display_numbers = T)

# Cluster 5 is similar to Cluster 0 and 2, Cluster 6 is similar to 1 and 3
# Group small clusters with the largest similar cluster
Spatial_CRC_integrated$sub.cluster[Spatial_CRC_integrated$sub.cluster == "Tumor_5"] <- "Tumor_0"
Spatial_CRC_integrated$sub.cluster[Spatial_CRC_integrated$sub.cluster == "Tumor_6"] <- "Tumor_1"

DimPlot(Spatial_CRC_integrated, group.by = "sub.cluster", cols = Palette_sub)

cells.show <- WhichCells(Spatial_CRC_integrated, idents = "Tumor")
cells.show <- intersect(cells.show,colnames(Spatial_CRC_integrated)[Spatial_CRC_integrated@reductions$umap@cell.embeddings[,2] > -2])

DimPlot(Spatial_CRC_integrated, group.by = "sub.cluster", cols = Palette_sub[11:15], cells = cells.show)

Spatial_CRC_integrated$Tumor <- "Stroma"
Spatial_CRC_integrated$Tumor[Spatial_CRC_integrated$Celltype == "Tumor"] <- "Tumor"

# Overview per sample
Spatial_CRC_integrated@meta.data %>% dplyr::select(sub.cluster, Sample, Study, Tumor) %>% group_by(Study, Sample, sub.cluster, Tumor) %>% 
  filter(sub.cluster != "Poor quality") %>%
  summarize(N = n()) %>%
  ggplot(aes(x = Sample, y = N, fill = sub.cluster)) + geom_bar(stat = "identity", position = "fill", width = 0.75) + theme_classic() + RotatedAxis() +
  xlab("") + facet_grid(Tumor~ Study, scales = "free_x", space = "free") + theme(strip.background = element_blank(), panel.border = element_rect(colour = "black", fill=NA, size=1)) + 
  ylab("Fraction per sample") + scale_fill_manual("Celltype",values = Palette_sub[-8]) +
  scale_y_continuous(labels = function(x) paste0(x*100, "%"),expand = c(0,0))

# Compare between liver met
Spatial_CRC_integrated$Site <- "Primary"
Spatial_CRC_integrated$Site[grep("LM", Spatial_CRC_integrated$Sample)] <- "LiverMet"
Spatial_CRC_integrated$Site <- factor(Spatial_CRC_integrated$Site, levels = c("Primary", "LiverMet"))

Spatial_CRC_integrated@meta.data %>% dplyr::select(sub.cluster, Site, Study, Tumor) %>% group_by(sub.cluster, Site, Tumor) %>% 
  filter(sub.cluster != "Poor quality") %>%
  summarize(N = n()) %>%
  ggplot(aes(x = Site, y = N, fill = sub.cluster)) + geom_bar(stat = "identity", position = "fill", width = 0.75) + theme_classic() + RotatedAxis() +
  xlab("") + facet_grid(~Tumor, scales = "free_x", space = "free") + theme(strip.background = element_blank(), panel.border = element_rect(colour = "black", fill=NA, size=1)) + 
  ylab("Fraction per sample") + scale_fill_manual("Celltype",values = Palette_sub[-8]) 



Spatial_CRC_integrated@meta.data %>% dplyr::select(sub.cluster, Sample, Study, Tumor) %>%
  filter(sub.cluster != "Poor quality", Tumor == "Tumor") %>%
  mutate(Sample = factor(Sample), sub.cluster = factor(sub.cluster)) %>%
  count(Sample, sub.cluster, .drop = FALSE) %>%
  group_by(Sample )  %>% mutate(freq = n/sum(n)) %>%
  group_by(sub.cluster) %>% summarize(Mean = mean(freq), SD = sd(freq), Total_Samples = sum(n>0))

### Evaluate Tumor sub in niches ####
library(Seurat)
library(ggplot2)
library(dplyr)
library(reshape2)
#library(sf)
#Markers.all <- FindAllMarkers(Spatial_CRC_integrated, logfc.threshold = 0.5, only.pos = T)
#write.csv(Markers.all, "/mnt/data/Projects/SW480_spheres/Out/Spatial_Markers_Celltype.csv")
#Spatial_CRC_integrated <- SetIdent(Spatial_CRC_integrated, value = "sub.cluster")
#Markers.all.sub <- FindAllMarkers(Spatial_CRC_integrated, logfc.threshold = 0.5, only.pos = T)
#write.csv(Markers.all.sub, "/mnt/data/Projects/SW480_spheres/Out/Spatial_Markers_Tumor_subcluster.csv")
Markers.all.sub <- read.csv("/mnt/data/Projects/SW480_spheres/Out/Spatial_Markers_Tumor_subcluster.csv", row.names=1)

Markers.tumor.sub <- Markers.all.sub %>% filter(p_val_adj < 0.05, cluster %in% c("Tumor_0", "Tumor_1","Tumor_2", "Tumor_3","Tumor_4","Tumor_5"))


GOI <- c("ASCL2","HOXA2", "HOXA3", "LGR5", "MYC", "VIM", "SPARC","RUNX2", "ZEB1" ,"CD44","ESRP1", "EMP1","MAL2","TACSTD2", "EPCAM","MUC2", "SPINK4","TFF2","AGR2" , "NOTUM","CCL20","LCN2", "CD24","REG3A", "RPS24")


Cluster_exp <- AverageExpression(Spatial_CRC_integrated, assays = "Spatial", features = GOI, group.by = "sub.cluster", slot = "data")
cal_z_score <- function(x){(x - mean(x)) / sd(x)}

Cluster_exp_z <- data.frame(t(apply(Cluster_exp$Spatial[,11:15], 1, cal_z_score)))

# Tumor marker plot
pheatmap::pheatmap(Cluster_exp_z, border_color = "black", cluster_cols = F, cluster_rows = F,
                   color  = viridis::inferno(10), angle_col = 0)

# new colors
pheatmap::pheatmap(Cluster_exp_z, border_color = "white", cluster_cols = F, cluster_rows = F, angle_col = 0)

# Make niche overview

# Filter out grid_id that include poor quality spots
Pct.dat <- Spatial_CRC_integrated@meta.data %>% 
  dplyr::select(sub.cluster, Niche, grid_id) %>% 
  filter(! grid_id %in% unique(Spatial_CRC_integrated$grid_id[Spatial_CRC_integrated$sub.cluster == "Poor quality"])) %>%
  group_by(sub.cluster, Niche) %>% summarise(N = n()) %>%
  mutate(pct = round(N/sum(N)*100, 2))


Pct.dat.heat <- Pct.dat %>% dplyr::select(-N) %>% dcast(sub.cluster ~ Niche) %>% replace(is.na(.), 0)
rownames(Pct.dat.heat) <- Pct.dat.heat$sub.cluster
Pct.dat.heat <- Pct.dat.heat[,-1]

pheatmap::pheatmap(Pct.dat.heat, angle_col = 90, display_numbers = T, treeheight_row = 15, treeheight_col = 15,
                   border_color = "black", number_color = "white",
                   color  = viridis::inferno(10)[2:8],
                   cluster_cols = T,
                   cluster_rows = F, gaps_row = 9)

# new colors
rownames(Pct.dat.heat)[10:14] <- c("CSC", "EMT", "HRC", "Secretoty", "Inflammatory")

pheatmap::pheatmap(Pct.dat.heat, angle_col = 90, display_numbers = T, treeheight_row = 15, treeheight_col = 0,
                   border_color = "black", number_color = "black",
                   cluster_cols = T,
                   cluster_rows = F, gaps_row = 9)

Spatial_CRC_integrated$tissue <- "Primary"
Spatial_CRC_integrated$tissue[grep("LM", Spatial_CRC_integrated$Sample)] <- "Liver met"
Spatial_CRC_integrated$tissue <- factor(Spatial_CRC_integrated$tissue, levels = c("Primary", "Liver met"))

DimPlot(Spatial_CRC_integrated, group.by = "sub.cluster", cols = Palette_sub, split.by = "tissue")

# Compare tumor sub clusters between primary and liver met
Spatial_CRC_integrated@meta.data %>% group_by(tissue, sub.cluster) %>% summarize(N = n()) %>%
  filter(sub.cluster %in% c("Tumor_0","Tumor_1","Tumor_2","Tumor_3","Tumor_4" )) %>%
  mutate(freq = N/sum(N)*100) %>%
  ggplot(aes(x = tissue, y = freq, fill = sub.cluster)) + geom_bar(stat = "identity", width = 0.7) + 
  theme_classic() + xlab("") + ylab("% of Tumor") + scale_fill_manual(values = Palette_sub[11:15])


# Compare neighbourhoods
Nhood_dat <- Spatial_CRC_integrated@meta.data %>% 
  dplyr::select(sub.cluster, Niche, grid_id) %>% 
  filter(! grid_id %in% unique(Spatial_CRC_integrated$grid_id[Spatial_CRC_integrated$sub.cluster == "Poor quality"])) %>%
  filter(grid_id %in% unique(Spatial_CRC_integrated$grid_id[Spatial_CRC_integrated$sub.cluster %in% c("Tumor_0", "Tumor_1", "Tumor_2", "Tumor_3","Tumor_4")]))

Nhood_dat %>% group_by(Niche, sub.cluster) %>% summarize(N = n()) %>%
  filter(sub.cluster %in% c("Tumor_0", "Tumor_1", "Tumor_2", "Tumor_3","Tumor_4")) %>%
  ggplot(aes(x = Niche, y = N, fill = sub.cluster)) + geom_bar(stat = "identity", position = "fill") +
  theme_classic() +  RotatedAxis() + scale_fill_manual(values = Palette_sub[11:15]) + ylab("Fraction of Tumor")


#### Single Examples: Spatial Plots ####

#Spatial_metadata <- readr::read_csv("/mnt/data/Projects/SW480_spheres/Out/Spatial_metadata_26oct.csv")
Spatial_metadata <- readr::read_csv("/mnt/data/Projects/SW480_spheres/Out/Spatial_metadata_09jan24.csv")

Palette_sub <- c("#FFD4D4", "#CFB997", "#E3F6FF", "#8DCBE6", "#EF9A53", "#E3ACF9", "#EAE0DA", "lightgrey", "#ABC270", "#658864", "#1B998B", "firebrick", "orange", "#93827F","#BCD8B7",   "blue", "magenta", "pink")

Palette_col <- c("firebrick", "red", "#7C96AB","orange", "#1B998B","#FFD4D4", "#CFB997", "#E3F6FF","#8DCBE6", "#EF9A53","#E3ACF9","#EAE0DA","lightgrey","#ABC270","#658864" ,"#CD0404")

Palette_col <- c("firebrick", "red", "grey","orange", "#1B998B","#FFD4D4", "#CFB997", "#E3F6FF","#8DCBE6", "#EF9A53","#E3ACF9","#EAE0DA","lightgrey","#ABC270","#658864" ,"#CD0404")


# Load Sample: CRC3
CRC3 <- readRDS("/mnt/data/Projects/SW480_spheres/Data/Qi_Spatial/CRC3_ST_SeuratObject.rds.gz")

# Add Meta data
meta_3 <- data.frame(Spatial_metadata[Spatial_metadata$Sample == "CRC3",])
rownames(meta_3) <- gsub("_1_1_1", "", meta_3$Row.names)
CRC3 <- AddMetaData(CRC3, metadata = meta_3)
SpatialDimPlot(CRC3, group.by = "sub.cluster") + scale_fill_manual(values = Palette_col[c(1,2,4:6, 8:9,11:15)])
SpatialDimPlot(CRC3, group.by = "DefineTypes") 
SpatialDimPlot(CRC3, group.by = "DefineTypes", alpha = 0) 
CRC3_Tumor <- subset(CRC3, subset = Celltype == "Tumor")

SpatialDimPlot(CRC3_Tumor, group.by = "sub.cluster") + scale_fill_manual(values = Palette_sub[11:15])
SpatialDimPlot(CRC3_Tumor, group.by = "sub.cluster", alpha = 0) + scale_fill_manual(values = Palette_col[c(1,2,4:6, 8:9,11:15)])
SpatialFeaturePlot(CRC3_Tumor, features = c("LGR5", "ZEB1", "EMP1"))

# Load Sample: CRC5
CRC5 <- Load10X_Spatial(data.dir = "/mnt/data/Projects/SW480_spheres/Data/Wu_Spatial/ST/ST-colon1")

# Add Meta data
meta_5 <- data.frame(Spatial_metadata[Spatial_metadata$Sample == "CRC5",])
rownames(meta_5) <- gsub("_1_1_1_1_2_1", "", meta_5$Row.names)
CRC5 <- AddMetaData(CRC5, metadata = meta_5)

# Remove hepatocytes
CRC5 <- subset(CRC5, subset = Celltype != "Hepatocytes")
SpatialDimPlot(CRC5, group.by = "Celltype") + scale_fill_manual(values = Palette_col[c(1,2,4:6, 8:9,11:15)])

# Subset tumor
CRC5_Tumor <- subset(CRC5, subset = Celltype == "Tumor")

SpatialDimPlot(CRC5_Tumor, group.by = "sub.cluster") + scale_fill_manual(values = Palette_sub[11:15])
SpatialDimPlot(CRC5_Tumor, group.by = "sub.cluster", alpha = 0) + scale_fill_manual(values = Palette_col[c(1,2,4:6, 8:9,11:15)])
SpatialFeaturePlot(CRC5_Tumor, features = c("LGR5", "ZEB1", "EMP1"))


# Load Sample: CRC8
CRC8 <- Load10X_Spatial(data.dir = "/mnt/data/Projects/SW480_spheres/Data/Wu_Spatial/ST/ST-colon4")

# Add Meta data
meta_8 <- data.frame(Spatial_metadata[Spatial_metadata$Sample == "CRC8",])
rownames(meta_8) <- gsub("_2_1_1_2_1", "", meta_8$Row.names)
CRC8 <- AddMetaData(CRC8, metadata = meta_8)

# Remove hepatocytes
CRC8 <- subset(CRC8, subset = Celltype != "Hepatocytes")
SpatialDimPlot(CRC8, group.by = "Celltype") + scale_fill_manual(values = Palette_col[c(1:6, 8:9,11:15)])

# Subset tumor
CRC8_Tumor <- subset(CRC8, subset = Celltype == "Tumor")
SpatialDimPlot(CRC8_Tumor, group.by = "sub.cluster") + scale_fill_manual(values = Palette_sub[11:15])
SpatialDimPlot(CRC8_Tumor, group.by = "sub.cluster", alpha = 0) + scale_fill_manual(values = Palette_col[c(1,2,4:6, 8:9,11:15)])

# Load Sample: CRC9_R1
CRC9_R1 <- Load10X_Spatial(data.dir = "/mnt/data/Projects/SW480_spheres/Data/Valdeolivas_Spatial/SN048_A121573_Rep1/SN048_A121573_Rep1")

# Add Meta data
meta_9R1 <- data.frame(Spatial_metadata[Spatial_metadata$Sample == "CRC9_R1",])
rownames(meta_9R1) <- gsub("_1_1_1_1_1_1_1_2", "", meta_9R1$Row.names)
CRC9_R1 <- AddMetaData(CRC9_R1, metadata = meta_9R1)

# Remove hepatocytes
CRC9_R1 <- subset(CRC9_R1, subset = Celltype != "Hepatocytes")
SpatialDimPlot(CRC9_R1, group.by = "Celltype") + scale_fill_manual(values = Palette_col[c(1:6, 8:9,11:15)])

# Subset tumor
CRC9_R1_Tumor <- subset(CRC9_R1, subset = Celltype == "Tumor")
SpatialDimPlot(CRC9_R1_Tumor, group.by = "sub.cluster") + scale_fill_manual(values = Palette_sub[11:15])
SpatialDimPlot(CRC9_R1_Tumor, group.by = "sub.cluster", alpha = 0) + scale_fill_manual(values = Palette_col[c(1,2,4:6, 8:9,11:15)])
