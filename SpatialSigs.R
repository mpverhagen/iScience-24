#### MAGIC Imputation #####
setwd("/mnt/data/")
library(Seurat)
library(ggplot2)
library(dplyr)
reticulate::py_discover_config(required_module = "magic")
reticulate::import("magic")
library(Rmagic)
# Load data
Palette_sub <- c("#FFD4D4", "#CFB997", "#E3F6FF", "#8DCBE6", "#EF9A53", "#E3ACF9", "#EAE0DA", "lightgrey", "#ABC270", "#658864", "#1B998B", "firebrick", "orange", "#93827F","#BCD8B7",   "blue", "magenta", "pink")

Spatial_CRC_tumor <- readRDS("/mnt/data/Projects/SW480_spheres/Out/Seurat/Spatial_CRC_tumor_jan24.rds")
Spatial_CRC_tumor <- Spatial_CRC_tumor %>% NormalizeData() %>% FindVariableFeatures() %>% ScaleData()
DimPlot(Spatial_CRC_tumor, group.by = "sub.cluster", cols = Palette_sub[11:15])

# Take counts and normalize
data <- library.size.normalize(t(Spatial_CRC_tumor@assays$integrated@data))
if (TRUE) {data <- sqrt(data)}

# Run magic
Magic_out <- Rmagic::magic(data, n.jobs = 10, knn = 15)
Magic_out <- Magic_out$result

# Magic-imputed, z-scores, linear scaled [0,1] signatures
cal_z_score <- function(x){(x - mean(x)) / sd(x)}
range01 <- function(x){(x-min(x))/(max(x)-min(x))}

# HRC
Sig_dat <- Magic_out[,colnames(Magic_out) %in% HRC_signature]
Sig_dat <- data.frame(apply(t(Sig_dat), 1, cal_z_score))
Spatial_CRC_tumor$HRC_Magic <- range01(rowMeans(Sig_dat, na.rm = TRUE))

# EMT
Sig_dat <- Magic_out[,colnames(Magic_out) %in% unique(Epcam_paper_signatures$EMT_Nanostring_Mes[!is.na(Epcam_paper_signatures$EMT_Nanostring_Mes)])]
Sig_dat <- data.frame(apply(t(Sig_dat), 1, cal_z_score))
Spatial_CRC_tumor$EMT_Magic <- range01(rowMeans(Sig_dat, na.rm = TRUE))

# LGR5
LGR5_sig <- c("LGR5", "SMOC2", "ASCL2", "NOTCH1","OLFM4", "PROM1", "ALCAM", "POLR1A", "STMN1")
Sig_dat <- Magic_out[,colnames(Magic_out) %in% LGR5_sig]
Sig_dat <- data.frame(apply(t(Sig_dat), 1, cal_z_score))
Spatial_CRC_tumor$CSC_Magic <- range01(rowMeans(Sig_dat, na.rm = TRUE))

FeaturePlot(Spatial_CRC_tumor, features = c("EMT_Magic", "HRC_Magic", "CSC_Magic"), ncol = 1) &
  scale_colour_gradientn(colours = rev(RColorBrewer::brewer.pal(n = 11, name = "RdBu"))) & theme_void()


