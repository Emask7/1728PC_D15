library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(cowplot)

# cohort 1 (anti-CD8): 36181 36201 35310 35289
# cohort 2 (control):  36051 36178 35309 35965

# load the datasets and create integrated dataset
  bn36181.data <- Read10X(data.dir = "./data/Bn_36181/")
  bn36181 <- CreateSeuratObject(counts = bn36181.data, project = "1728PC",
                                min.cells = 3, min.features = 200)
  bn36181[["bn"]] <- "bn36181"
  bn36181[["cohort"]] <- "anti-CD8"

  bn36201.data <- Read10X(data.dir = "./data/Bn_36201/")
  bn36201 <- CreateSeuratObject(counts = bn36201.data, project = "1728PC",
                                min.cells = 3, min.features = 200)
  bn36201[["bn"]] <- "bn36201"
  bn36201[["cohort"]] <- "anti-CD8"

  bn36310.data <- Read10X(data.dir = "./data/Bn_36310/")
  bn36310 <- CreateSeuratObject(counts = bn36310.data, project = "1728PC",
                                min.cells = 3, min.features = 200)
  bn36310[["bn"]] <- "bn36310"
  bn36310[["cohort"]] <- "anti-CD8"

  bn35289.data <- Read10X(data.dir = "./data/Bn_35289/")
  bn35289 <- CreateSeuratObject(counts = bn35289.data, project = "1728PC",
                                min.cells = 3, min.features = 200)
  bn35289[["bn"]] <- "bn35289"
  bn35289[["cohort"]] <- "anti-CD8"

  bn36051.data <- Read10X(data.dir = "./data/Bn_36051/")
  bn36051 <- CreateSeuratObject(counts = bn36051.data, project = "1728PC",
                                min.cells = 3, min.features = 200)
  bn36051[["bn"]] <- "bn36051"
  bn36051[["cohort"]] <- "control"

  bn36178.data <- Read10X(data.dir = "./data/Bn_36178/")
  bn36178 <- CreateSeuratObject(counts = bn36178.data, project = "1728PC",
                                min.cells = 3, min.features = 200)
  bn36178[["bn"]] <- "bn36178"
  bn36178[["cohort"]] <- "control"

  bn35309.data <- Read10X(data.dir = "./data/Bn_35309/")
  bn35309 <- CreateSeuratObject(counts = bn35309.data, project = "1728PC",
                                min.cells = 3, min.features = 200)
  bn35309[["bn"]] <- "bn35309"
  bn35309[["cohort"]] <- "control"

  bn35965.data <- Read10X(data.dir = "./data/Bn_35965/")
  bn35965 <- CreateSeuratObject(counts = bn35965.data, project = "1728PC",
                                min.cells = 3, min.features = 200)
  bn35965[["bn"]] <- "bn35965"
  bn35965[["cohort"]] <- "control"


# not sure how to organize bn.list yet
  bn.list <- list(bn36181 = bn36181, bn36201 = bn36201,
                  bn36310 = bn36310, bn35289 = bn35289,
                  bn36051 = bn36051, bn36178 = bn36178,
                  bn35309 = bn35309, bn35965 = bn35965)

  bn.list <- list(antiCD8 = c(bn36181, bn36201, bn36310, bn35289),
                  control = c(bn36051, bn36178, bn35309, bn35965))


# Normalize datasets individually by SCTransform(), instead of NormalizeData()
bn.list <- lapply(X = bn.list, FUN = SCTransform)
features <- SelectIntegrationFeatures(object.list = bn.list, nfeatures = 3000)

# Run the PrepSCTIntegration() function prior to identifying anchors
bn.list <- PrepSCTIntegration(object.list = bn.list, anchor.features = features)

# Find integration anchors and integrate the data
immune.anchors <- FindIntegrationAnchors(object.list = bn.list, normalization.method = "SCT",
                                         anchor.features = features)
immune.combined.sct <- IntegrateData(anchorset = immune.anchors, normalization.method = "SCT")

# Perform an integrated analysis
# specify that we will perform downstream analysis on the corrected data note that the
# original unmodified data still resides in the 'RNA' assay
DefaultAssay(immune.combined.sct) <- "integrated"

# Run the standard workflow for visualization and clustering
immune.combined.sct <- ScaleData(immune.combined.sct, verbose = FALSE)
immune.combined.sct <- RunPCA(immune.combined.sct, verbose = FALSE)

immune.combined.sct.UMAP <- RunUMAP(immune.combined.sct, reduction = "pca", dims = 1:30)
immune.combined.sct.UMAP <- FindNeighbors(immune.combined.sct.UMAP, reduction = "pca", dims = 1:30)
immune.combined.sct.UMAP <- FindClusters(immune.combined.sct.UMAP, resolution = 0.5)

immune.combined.sct.tSNE <- RunTSNE(immune.combined.sct)
immune.combined.sct.tSNE <- FindNeighbors(immune.combined.sct.tSNE, reduction = "pca", dims = 1:30)
immune.combined.sct.tSNE <- FindClusters(immune.combined.sct.tSNE, resolution = 0.5)

# Visualization
p1 <- DimPlot(immune.combined.sct.UMAP, reduction = "umap", group.by = "bn")
p2 <- DimPlot(immune.combined.sct.UMAP, reduction = "umap", label = TRUE, repel = TRUE)
p1 + p2

p1 <- DimPlot(immune.combined.sct.tSNE, reduction = "tsne", group.by = "bn")
p2 <- DimPlot(immune.combined.sct.tSNE, reduction = "tsne", label = TRUE, repel = TRUE)
p1 + p2

DimPlot(immune.combined.sct.UMAP, reduction = "umap", split.by = "bn")

# Identify conserved cell type markers
# For performing differential expression after integration, we use the original data
DefaultAssay(immune.combined.sct.UMAP) <- "RNA"
ident.0.markers <- FindConservedMarkers(immune.combined.sct.UMAP, ident.1 = 0, grouping.var = "bn", verbose = FALSE)
head(ident.0.markers)


fcm <- function(x) {
  FindConservedMarkers(immune.combined.sct.UMAP, ident.1 = x, grouping.var = "bn", verbose = TRUE)
}

conserved.markers <- list(
  ident.0 = fcm(0), ident.1 = fcm(1), ident.2 = fcm(2), ident.3 = fcm(3),
  ident.4 = fcm(4), ident.5 = fcm(5), ident.6 = fcm(6), ident.7 = fcm(7),
  ident.8 = fcm(8), ident.9 = fcm(9), ident.10 = fcm(10), ident.11 = fcm(11),
  ident.12 = fcm(12), ident.13 = fcm(13), ident.14 = fcm(14),
  ident.15 = fcm(15), ident.16 = fcm(16)
)


rm(x)
rm(y)
rm(z)
rm(conserved.markers.arr)
rm(m)


library(openxlsx)

wb <- createWorkbook("Conserved Markers.xlsx")

for (n in 0:16) {
  wsName <- stringi::stri_join("cluster ", n)
  addWorksheet(wb, wsName)
}

writeData(wb, "cluster 0", conserved.markers$ident.0, rowNames = TRUE)
writeData(wb, "cluster 1", conserved.markers$ident.1, rowNames = TRUE)
writeData(wb, "cluster 2", conserved.markers$ident.2, rowNames = TRUE)
writeData(wb, "cluster 3", conserved.markers$ident.3, rowNames = TRUE)
writeData(wb, "cluster 4", conserved.markers$ident.4, rowNames = TRUE)
writeData(wb, "cluster 5", conserved.markers$ident.5, rowNames = TRUE)
writeData(wb, "cluster 6", conserved.markers$ident.6, rowNames = TRUE)
writeData(wb, "cluster 7", conserved.markers$ident.7, rowNames = TRUE)
writeData(wb, "cluster 8", conserved.markers$ident.8, rowNames = TRUE)
writeData(wb, "cluster 9", conserved.markers$ident.9, rowNames = TRUE)
writeData(wb, "cluster 10", conserved.markers$ident.10, rowNames = TRUE)
writeData(wb, "cluster 11", conserved.markers$ident.11, rowNames = TRUE)
writeData(wb, "cluster 12", conserved.markers$ident.12, rowNames = TRUE)
writeData(wb, "cluster 13", conserved.markers$ident.13, rowNames = TRUE)
writeData(wb, "cluster 14", conserved.markers$ident.14, rowNames = TRUE)
writeData(wb, "cluster 15", conserved.markers$ident.15, rowNames = TRUE)
writeData(wb, "cluster 16", conserved.markers$ident.16, rowNames = TRUE)

saveWorkbook(wb, "Conserved Markers.xlsx", overwrite = TRUE)










# need to identify clusters before we can do this part
# immune.combined.sct.UMAP <- RenameIdents(immune.combined.sct.UMAP, `0` = "CD14 Mono", `1` = "CD4 Naive T", `2` = "CD4 Memory T",
#                                 `3` = "CD16 Mono", `4` = "B", `5` = "CD8 T", `6` = "NK", `7` = "T activated", `8` = "DC", `9` = "B Activated",
#                                 `10` = "Mk", `11` = "pDC", `12` = "Eryth", `13` = "Mono/Mk Doublets", `14` = "HSPC")
# DimPlot(immune.combined.sct.UMAP, label = TRUE)

