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

#### Merge Spatial Transcriptomics Sets ####
Spatial_Qi <- readRDS(paste(path,"/Spatial_Seur_Qi.rds", sep =""))
Spatial_Wu <- readRDS(paste(path, "/Spatial_Seur_Wu.rds", sep = ""))
Spatial_Valdeolivas <- readRDS(paste(path, "Spatial_Seur_Valdeolivas.rds", sep = ""))

# Merge Objects
Spatial_CRC <- merge(Spatial_Qi, Spatial_Wu)
Spatial_CRC <- merge(Spatial_CRC, Spatial_Valdeolivas)
rm(list = c("Spatial_Qi", "Spatial_Wu", "Spatial_Valdeolivas"))
gc()

Spatial_CRC <- Spatial_CRC %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA() %>%
  FindNeighbors(dims = 1:30) %>%
  FindClusters(verbose = TRUE) %>%
  RunUMAP(dims = 1:30)

# plot UMAP embedding
DimPlot(Spatial_CRC, label=TRUE, repel = TRUE, reduction = "umap", group.by = "Sample") + NoLegend()


# make a dataframe containing the image coordinates for each sample
image_df <- do.call(rbind, lapply(names(Spatial_CRC@images), function(x){
  Spatial_CRC@images[[x]]@coordinates
}))

# merge the image_df with the Seurat metadata
new_meta <- merge(Spatial_CRC@meta.data[,-c(6:10)], image_df, by='row.names')

# fix the row ordering to match the original seurat object
rownames(new_meta) <- new_meta$Row.names
ix <- match(as.character(colnames(Spatial_CRC)), as.character(rownames(new_meta)))
new_meta <- new_meta[ix,]

# add the new metadata to the seurat object
Spatial_CRC@meta.data <- new_meta


Spatial_CRC@meta.data %>% group_by(Sample) %>% summarize(N = n())
Spatial_CRC$Sample <- factor(Spatial_CRC$Sample, levels = c("CRC1", "CRC2", "CRC3", "CRC4",
                                                            "CRC5", "CRC5_LM", "CRC6", "CRC6_LM",
                                                            "CRC7", "CRC7_LM", "CRC8", "CRC8_LM",
                                                            "CRC9_R1", "CRC9_R2", "CRC10_R1", "CRC10_R2",
                                                            "CRC11_R1", "CRC11_R2", "CRC12_R1", "CRC12_R2",
                                                            "CRC13_R1", "CRC13_R2", "CRC14_R1", "CRC14_R2",
                                                            "CRC15_R1", "CRC15_R2"
))
Spatial_CRC$Study <- factor(Spatial_CRC$Study, levels = c("Qi", "Wu", "Valdeolivas"))

# Spatial ggplot2
Spatial_CRC@meta.data %>% 
  ggplot(aes(x = imagecol, y = imagerow, color = Annotation)) + geom_point(size = 0.1) + facet_wrap(~ Sample, ncol = 5, scales = "free") +
  theme_void()

# QC plot I
Spatial_CRC@meta.data %>% dplyr::select(nCount_Spatial, nFeature_Spatial, Sample, Study) %>% group_by(Study, Sample) %>% 
  summarize(nCount_Spatial = mean(nCount_Spatial),nFeature_Spatial =median(nFeature_Spatial)) %>%
  ggplot(aes(x = Sample, y = nCount_Spatial)) + geom_bar(stat = "identity") + theme_classic() + RotatedAxis() +
  xlab("") + facet_grid(~ Study, scales = "free_x", space = "free") + theme(strip.background = element_blank()) + 
  ylab("Mean Reads per Spot")

# QC plot II
Spatial_CRC@meta.data %>% dplyr::select(nCount_Spatial, nFeature_Spatial, Sample, Study) %>% group_by(Study, Sample) %>% 
  summarize(nCount_Spatial = mean(nCount_Spatial),nFeature_Spatial =median(nFeature_Spatial)) %>%
  ggplot(aes(x = Sample, y = nFeature_Spatial)) + geom_bar(stat = "identity") + theme_classic() + RotatedAxis() +
  xlab("") + facet_grid(~ Study, scales = "free_x", space = "free") + theme(strip.background = element_blank()) + 
  ylab("Median Features per Spot")

#### Rpca integration of Spatial Data ####
# split the dataset into a list of two seurat objects (stim and CTRL)
CRC.list <- SplitObject(Spatial_CRC, split.by = "Sample")

# normalize and identify variable features for each dataset independently
CRC.list <- lapply(X = CRC.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

# select features that are repeatedly variable across datasets for integration run PCA on each dataset using these features
features <- SelectIntegrationFeatures(object.list = CRC.list)
CRC.list <- lapply(X = CRC.list, FUN = function(x) {
  x <- ScaleData(x, features = features, verbose = TRUE)
  x <- RunPCA(x, features = features, verbose = FALSE)
})

CRC.anchors <- FindIntegrationAnchors(object.list = CRC.list, anchor.features = features, reduction = "rpca")

# Create an 'integrated' data assay
Spatial_CRC_integrated <- IntegrateData(anchorset = CRC.anchors)

# preprocess RNA assay
DefaultAssay(Spatial_CRC_integrated) <- "Spatial"
Spatial_CRC_integrated <- Spatial_CRC_integrated %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData()

DefaultAssay(Spatial_CRC_integrated) <- "integrated"

# Run the standard workflow for visualization and clustering
Spatial_CRC_integrated <- ScaleData(Spatial_CRC_integrated, verbose = FALSE)
Spatial_CRC_integrated <- RunPCA(Spatial_CRC_integrated, npcs = 50, verbose = FALSE)
Spatial_CRC_integrated <- FindNeighbors(Spatial_CRC_integrated, reduction = "pca", dims = 1:30)
Spatial_CRC_integrated <- RunUMAP(Spatial_CRC_integrated, reduction = "pca", dims = 1:50, min.dist = 0.2, n.neighbors = 100, spread = 2)
Spatial_CRC_integrated <- FindClusters(Spatial_CRC_integrated, resolution = 0.2)

DimPlot(Spatial_CRC_integrated, reduction = "umap", group.by = "seurat_clusters", label = TRUE, repel = TRUE)

# Marker based annotation
DefaultAssay(Spatial_CRC_integrated) <- "Spatial"
GOI <- c("ASCL2", "EPCAM", "VIL1",
         "MMP11","POSTN","THBS2","SPINK4", "OLFM4", "MUC2", "VIM","VWF","SPARCL1","PTGDS",
         "IGLC2","IGHG3","IGHG4",  "ACTG2", "DES", "MYH11",
         "TF","ALDOB", "TTR","CXCL8", "CXCL5", "CXCL1",  "COL1A1", "COL3A1", "COL1A2",
         "VIP", "GAL", "UCHL1")

DotPlot(Spatial_CRC_integrated, group.by = "seurat_clusters", features = GOI) + coord_flip() + xlab("") + ylab("") +
  scale_color_viridis_c()


# Annotation
Spatial_CRC_integrated$Celltype <- "Other"
Spatial_CRC_integrated$Celltype[Spatial_CRC_integrated$seurat_clusters == "0"] <- "Poor quality"
Spatial_CRC_integrated$Celltype[Spatial_CRC_integrated$seurat_clusters == "1"] <- "Tumor"
Spatial_CRC_integrated$Celltype[Spatial_CRC_integrated$seurat_clusters == "2"] <- "Stromal"
Spatial_CRC_integrated$Celltype[Spatial_CRC_integrated$seurat_clusters == "3"] <- "Normal epithelium"
Spatial_CRC_integrated$Celltype[Spatial_CRC_integrated$seurat_clusters == "4"] <- "Endothelial"
Spatial_CRC_integrated$Celltype[Spatial_CRC_integrated$seurat_clusters == "5"] <- "Immune"
Spatial_CRC_integrated$Celltype[Spatial_CRC_integrated$seurat_clusters == "6"] <- "Smooth muscle"
Spatial_CRC_integrated$Celltype[Spatial_CRC_integrated$seurat_clusters == "7"] <- "Hepatocytes"
Spatial_CRC_integrated$Celltype[Spatial_CRC_integrated$seurat_clusters == "8"] <- "Macrophages"
Spatial_CRC_integrated$Celltype[Spatial_CRC_integrated$seurat_clusters == "9"] <- "Myofibroblast"
Spatial_CRC_integrated$Celltype[Spatial_CRC_integrated$seurat_clusters == "10"] <- "Neuron"

# Annotation
Spatial_CRC_integrated$Annotation <- "Other"
Spatial_CRC_integrated$Annotation[Spatial_CRC_integrated$Celltype %in% c("Endothelial","Stromal", "Smooth muscle", "Myofibroblast")] <- "Stroma"
Spatial_CRC_integrated$Annotation[Spatial_CRC_integrated$Celltype %in% c("Hepatocytes","Normal epithelium")] <- "Epithelium"
Spatial_CRC_integrated$Annotation[Spatial_CRC_integrated$Celltype %in% c("Poor quality")] <- "Poor quality"
Spatial_CRC_integrated$Annotation[Spatial_CRC_integrated$Celltype %in% c("Immune", "Macrophages")] <- "Immune"
Spatial_CRC_integrated$Annotation[Spatial_CRC_integrated$Celltype %in% c("Neuron")] <- "ENS"
Spatial_CRC_integrated$Annotation[Spatial_CRC_integrated$Celltype %in% c("Tumor")] <- "Tumor"

DotPlot(Spatial_CRC_integrated, group.by = "Celltype", features = GOI) + coord_flip() + xlab("") + ylab("") +
  scale_color_viridis_c() + RotatedAxis()

# Palette
Palette <- c("#FFD4D4", "#CFB997", "#E3F6FF", "#8DCBE6", "#EF9A53", "#E3ACF9", "#EAE0DA", "lightgrey", "#ABC270", "#658864", "#CD0404")
Palette_2 <- c("#FFD4D4", "#CFB997",  "#8DCBE6", "lightgrey", "#ABC270",  "#CD0404")

DimPlot(Spatial_CRC_integrated, group.by = "Celltype", cols = Palette)

# Spatial ggplot2
Spatial_CRC_integrated@meta.data %>% 
  ggplot(aes(x = imagecol, y = imagerow, color = Celltype)) + geom_point(size = 0.1) + facet_wrap(~ Sample, ncol = 5, scales = "free") +
  theme_void() + scale_color_manual(values = Palette)

# Overview of annotation per celltype
Spatial_CRC_integrated@meta.data %>% dplyr::select(Celltype, Sample, Study) %>% group_by(Study, Sample, Celltype) %>% 
  summarize(N = n()) %>%
  ggplot(aes(x = Sample, y = N, fill = Celltype)) + geom_bar(stat = "identity", position = "fill") + theme_classic() + RotatedAxis() +
  xlab("") + facet_grid(~ Study, scales = "free_x", space = "free") + theme(strip.background = element_blank()) + 
  ylab("Distribution of Cell types") + scale_fill_manual(values = Palette)

# Overview of annotation per annotation
Spatial_CRC_integrated@meta.data %>% dplyr::select(Annotation, Sample, Study) %>% group_by(Study, Sample, Annotation) %>% 
  summarize(N = n()) %>%
  ggplot(aes(x = Sample, y = N, fill = Annotation)) + geom_bar(stat = "identity", position = "fill") + theme_classic() + RotatedAxis() +
  xlab("") + facet_grid(~ Study, scales = "free_x", space = "free") + theme(strip.background = element_blank()) + 
  ylab("Distribution of Cell types") + scale_fill_manual(values = Palette_2)

# Overview of annotation per celltype: filter poor quality
Spatial_CRC_integrated@meta.data %>% dplyr::select(Celltype, Sample, Study) %>% group_by(Study, Sample, Celltype) %>% filter(Celltype != "Poor quality") %>%
  summarize(N = n()) %>%
  ggplot(aes(x = Sample, y = N, fill = Celltype)) + geom_bar(stat = "identity", position = "fill") + theme_classic() + RotatedAxis() +
  xlab("") + facet_grid(~ Study, scales = "free_x", space = "free") + theme(strip.background = element_blank()) + 
  ylab("Distribution of Cell types") + scale_fill_manual(values = Palette[-8])

# Overview of annotation per celltype: filter poor quality
Spatial_CRC_integrated@meta.data %>% dplyr::select(Annotation, Sample, Study) %>% group_by(Study, Sample, Annotation) %>% filter(Annotation != "Poor quality") %>%
  summarize(N = n()) %>%
  ggplot(aes(x = Sample, y = N, fill = Annotation)) + geom_bar(stat = "identity", position = "fill") + theme_classic() + RotatedAxis() +
  xlab("") + facet_grid(~ Study, scales = "free_x", space = "free") + theme(strip.background = element_blank()) + 
  ylab("Distribution of Cell types") + scale_fill_manual(values = Palette_2[-4])



Spatial_CRC_integrated$Sample <- factor(Spatial_CRC_integrated$Sample, levels = c("CRC1", "CRC2", "CRC3", "CRC4",
                                                                                  "CRC5", "CRC5_LM", "CRC6", "CRC6_LM",
                                                                                  "CRC7", "CRC7_LM", "CRC8", "CRC8_LM",
                                                                                  "CRC9_R1", "CRC9_R2", "CRC10_R1", "CRC10_R2",
                                                                                  "CRC11_R1", "CRC11_R2", "CRC12_R1", "CRC12_R2",
                                                                                  "CRC13_R1", "CRC13_R2", "CRC14_R1", "CRC14_R2",
                                                                                  "CRC15_R1", "CRC15_R2"
))
Spatial_CRC_integrated$Study <- factor(Spatial_CRC_integrated$Study, levels = c("Qi", "Wu", "Valdeolivas"))

FeaturePlot(Spatial_CRC_integrated, features = c("nCount_Spatial", "nFeature_Spatial"))
# Visualization
p1 <- DimPlot(Spatial_CRC_integrated, reduction = "umap", group.by = "Sample")
p2 <- DimPlot(Spatial_CRC_integrated, reduction = "umap", group.by = "Annotation", label = TRUE,
              repel = TRUE)
p1 + p2


#saveRDS(Spatial_CRC_integrated, "/mnt/data/Projects/SW480_spheres/Out/Seurat/Spatial_merge_integrated.rds")

#### Spatial plotting ####
library(Seurat)
library(ggplot2)
library(dplyr)

setwd("/mnt/data/Projects/SW480_spheres/")
Palette <- c("#FFD4D4", "#CFB997", "#E3F6FF", "#8DCBE6", "#EF9A53", "#E3ACF9", "#EAE0DA", "lightgrey", "#ABC270", "#658864", "#CD0404")
Palette_2 <- c("#FFD4D4", "#CFB997",  "#8DCBE6", "lightgrey", "#ABC270",  "#CD0404")
Palette_3 <- c("lightgrey", "#CD0404")

Spatial_CRC_integrated <- readRDS("/mnt/data/Projects/SW480_spheres/Out/Seurat/Spatial_merge_integrated.rds")

Spatial_CRC_integrated$Group <- "Other"
Spatial_CRC_integrated$Group[Spatial_CRC_integrated$Annotation == "Tumor" & Spatial_CRC_integrated@reductions$umap@cell.embeddings[,2] > -2] <- "Tumor"

DimPlot(Spatial_CRC_integrated, group.by = "Celltype", cols = Palette)
DimPlot(Spatial_CRC_integrated, group.by = "Annotation", cols = Palette_2)
DimPlot(Spatial_CRC_integrated, group.by = "Group", cols = Palette_3)
DimPlot(Spatial_CRC_integrated, group.by = "Sample")

FeatureScatter(Spatial_CRC_integrated, feature1 = "nCount_Spatial", feature2 = "nFeature_Spatial", group.by = "Celltype", cols = Palette)

# Make Cell_type2
Spatial_CRC_integrated$CelltypeII <- Spatial_CRC_integrated$Celltype
DimPlot(Spatial_CRC_integrated, group.by = "CelltypeII", cols = Palette)
DimPlot(Spatial_CRC_integrated, group.by = "Celltype", cols = Palette)

# Annotation
Spatial_CRC_integrated$Annotation <- "Other"
Spatial_CRC_integrated$Annotation[Spatial_CRC_integrated$CelltypeII %in% c("Endothelial","Stromal", "Smooth muscle", "Myofibroblast")] <- "Stroma"
Spatial_CRC_integrated$Annotation[Spatial_CRC_integrated$CelltypeII %in% c("Hepatocytes","Normal epithelium")] <- "Epithelium"
Spatial_CRC_integrated$Annotation[Spatial_CRC_integrated$CelltypeII %in% c("Poor quality")] <- "Poor quality"
Spatial_CRC_integrated$Annotation[Spatial_CRC_integrated$CelltypeII %in% c("Immune", "Macrophages")] <- "Immune"
Spatial_CRC_integrated$Annotation[Spatial_CRC_integrated$CelltypeII %in% c("Neuron")] <- "ENS"
Spatial_CRC_integrated$Annotation[Spatial_CRC_integrated$CelltypeII %in% c("Tumor")] <- "Tumor"



# Spatial ggplot2
Spatial_CRC_integrated@meta.data %>% 
  ggplot(aes(x = imagecol, y = imagerow, color = Annotation)) + geom_point(size = 0.1) + facet_wrap(~ Sample, ncol = 5, scales = "free") +
  theme_void() + scale_color_manual(values = Palette_2)

# Overview of annotation per celltype
Spatial_CRC_integrated@meta.data %>% dplyr::select(CelltypeII, Sample, Study) %>% group_by(Study, Sample, CelltypeII) %>% 
  summarize(N = n()) %>%
  ggplot(aes(x = Sample, y = N, fill = CelltypeII)) + geom_bar(stat = "identity", position = "fill") + theme_classic() + RotatedAxis() +
  xlab("") + facet_grid(~ Study, scales = "free_x", space = "free") + theme(strip.background = element_blank()) + 
  ylab("Distribution of Cell types") + scale_fill_manual(values = Palette)

# Overview of annotation per annotation
Spatial_CRC_integrated@meta.data %>% dplyr::select(Annotation, Sample, Study) %>% group_by(Study, Sample, Annotation) %>% 
  summarize(N = n()) %>%
  ggplot(aes(x = Sample, y = N, fill = Annotation)) + geom_bar(stat = "identity", position = "fill") + theme_classic() + RotatedAxis() +
  xlab("") + facet_grid(~ Study, scales = "free_x", space = "free") + theme(strip.background = element_blank()) + 
  ylab("Distribution of Cell types") + scale_fill_manual(values = Palette_2) + geom_hline(yintercept = 0.2)

# Overview of annotation per celltype: filter poor quality
Spatial_CRC_integrated@meta.data %>% dplyr::select(Annotation, Sample, Study) %>% group_by(Study, Sample, Annotation) %>% filter(Annotation != "Poor quality") %>%
  summarize(N = n()) %>%
  ggplot(aes(x = Sample, y = N, fill = Annotation)) + geom_bar(stat = "identity", position = "fill") + theme_classic() + RotatedAxis() +
  xlab("") + facet_grid(~ Study, scales = "free_x", space = "free") + theme(strip.background = element_blank()) + 
  ylab("Distribution of Cell types") + scale_fill_manual(values = Palette_2[-4]) + geom_hline(yintercept = 0.2)


write.csv(Spatial_CRC_integrated@meta.data, "/mnt/data/Projects/SW480_spheres/Out/Spatial_CRC_integrated_metadata.csv")

