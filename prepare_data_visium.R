library(NMF)
library(Seurat)
library(Seurat)
library(SeuratData)
library(ggplot2)
library(patchwork)
library(dplyr)
rm(list = ls())

# 
# load(paste0("G:/My Drive/Dissertation/Spatial/data/LuLabSpatialData/LuLabSpatial_data0_",sample,".RData"))
# load(paste0("G:/My Drive/Dissertation/Spatial/data/LuLabSpatialData/LuLabSpatial_seuObj_",sample,".RData"))
sample <- "11"
visiumDataPath <- "C:/Users/whitealj/Downloads/More10Xdata/More10Xdata/Visium_Mouse_Olfactory_Bulb"
visiumDataFileName <- "Visium_Mouse_Olfactory_Bulb_filtered_feature_bc_matrix.h5"

#no update
visiumData <- Load10X_Spatial(
  data.dir = visiumDataPath,
  filename = visiumDataFileName,
  assay = 'Spatial',
  slice = 'A1',
  filter.matrix = TRUE,
  to.upper = FALSE
)
visiumData <- SCTransform(visiumData, assay = "Spatial", verbose = FALSE)

x <- visiumData@assays$Spatial@counts
coords <- visiumData@images$A1@coordinates[, 2:3]
dim(x)
dim(coords)


saveRDS(x,paste0("G:/My Drive/Dissertation/Spatial/spatial_analysis/data_08_16_2022/data_",sample,"_before.rds"))
saveRDS(coords,paste0("G:/My Drive/Dissertation/Spatial/spatial_analysis/data_08_16_2022/data_",sample,"_before_coords.rds"))

x <- visiumData@assays$SCT@counts
coords <- visiumData@images$A1@coordinates[, 2:3]
dim(x)
dim(coords)

saveRDS(x,paste0("G:/My Drive/Dissertation/Spatial/spatial_analysis/data_08_16_2022/data_",sample,"_after.rds"))
saveRDS(coords,paste0("G:/My Drive/Dissertation/Spatial/spatial_analysis/data_08_16_2022/data_",sample,"_after_coords.rds"))
