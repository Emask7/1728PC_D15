library(dplyr)
library(foreach)
library(Seurat)
library(patchwork)
library(ggplot2)
library(cowplot)

# cohort 1 (anti-CD8): 36181 36201 35310 35289
# cohort 2 (control):  36051 36178 35309 35965

# load the datasets and create integrated dataset
  bn36181 <- Read10X(data.dir = "./data/Bn_36181/")
  rownames(bn36181) <- gsub("_", "-", rownames(bn36181))
  bn36181 <- CreateSeuratObject(counts = bn36181, project = "1728PC",
                                min.cells = 3, min.features = 200)
  bn36181[["bn"]] <- "bn36181"
  bn36181[["cohort"]] <- "anti-CD8"

  bn36201 <- Read10X(data.dir = "./data/Bn_36201/")
  rownames(bn36201) <- gsub("_", "-", rownames(bn36201))
  bn36201 <- CreateSeuratObject(counts = bn36201, project = "1728PC",
                                min.cells = 3, min.features = 200)
  bn36201[["bn"]] <- "bn36201"
  bn36201[["cohort"]] <- "anti-CD8"

  bn35310 <- Read10X(data.dir = "./data/Bn_35310/")
  rownames(bn35310) <- gsub("_", "-", rownames(bn35310))
  bn35310 <- CreateSeuratObject(counts = bn35310, project = "1728PC",
                                min.cells = 3, min.features = 200)
  bn35310[["bn"]] <- "bn35310"
  bn35310[["cohort"]] <- "anti-CD8"

  bn35289 <- Read10X(data.dir = "./data/Bn_35289/")
  rownames(bn35289) <- gsub("_", "-", rownames(bn35289))
  bn35289 <- CreateSeuratObject(counts = bn35289, project = "1728PC",
                                min.cells = 3, min.features = 200)
  bn35289[["bn"]] <- "bn35289"
  bn35289[["cohort"]] <- "anti-CD8"

  bn36051 <- Read10X(data.dir = "./data/Bn_36051/")
  rownames(bn36051) <- gsub("_", "-", rownames(bn36051))
  bn36051 <- CreateSeuratObject(counts = bn36051, project = "1728PC",
                                min.cells = 3, min.features = 200)
  bn36051[["bn"]] <- "bn36051"
  bn36051[["cohort"]] <- "control"

  bn36178 <- Read10X(data.dir = "./data/Bn_36178/")
  rownames(bn36178) <- gsub("_", "-", rownames(bn36178))
  bn36178 <- CreateSeuratObject(counts = bn36178, project = "1728PC",
                                min.cells = 3, min.features = 200)
  bn36178[["bn"]] <- "bn36178"
  bn36178[["cohort"]] <- "control"

  bn35309 <- Read10X(data.dir = "./data/Bn_35309/")
  rownames(bn35309) <- gsub("_", "-", rownames(bn35309))
  bn35309 <- CreateSeuratObject(counts = bn35309, project = "1728PC",
                                min.cells = 3, min.features = 200)
  bn35309[["bn"]] <- "bn35309"
  bn35309[["cohort"]] <- "control"

  bn35965 <- Read10X(data.dir = "./data/Bn_35965/")
  rownames(bn35965) <- gsub("_", "-", rownames(bn35965))
  bn35965 <- CreateSeuratObject(counts = bn35965, project = "1728PC",
                                min.cells = 3, min.features = 200)
  bn35965[["bn"]] <- "bn35965"
  bn35965[["cohort"]] <- "control"

# not sure how to organize data.list yet
  cohort1 <- list(bn36181 = bn36181, bn36201 = bn36201,
                  bn35310 = bn35310, bn35289 = bn35289)

  cohort2 <- list(bn36051 = bn36051, bn36178 = bn36178,
                  bn35309 = bn35309, bn35965 = bn35965)

  # test <- list(bn36181 = bn36181, bn36201 = bn36201)

# Function to run Seurat integration analysis
  run_seurat <- function(data.list) {
    # Normalize datasets individually by SCTransform(), instead of NormalizeData()
      data.list <- lapply(X = data.list, FUN = SCTransform, method = "glmGamPoi")

    # Select integration features
      features <- SelectIntegrationFeatures(object.list = data.list, nfeatures = 3000)

    # Reciprocal PCA
      data.list <- lapply(X = data.list, FUN = function(x) {
        x <- ScaleData(x, features = features, verbose = FALSE)
        x <- RunPCA(x, features = features, verbose = FALSE)
      })

    # Run the PrepSCTIntegration() function prior to identifying anchors
      data.list <- PrepSCTIntegration(object.list = data.list, anchor.features = features)
      data.list <- lapply(X = data.list, FUN = RunPCA, features = features)

    # Find integration anchors and integrate the data
      anchors <- FindIntegrationAnchors(object.list = data.list, normalization.method = "SCT",
                                        anchor.features = features, dims = 1:30, reduction = "rpca", k.anchor = 20)
      data.combined <- IntegrateData(anchorset = anchors, normalization.method = "SCT", dims = 1:30)

    # Return integrated data
      data.combined
  }

# test.combined <- run_seurat(test)

  cohort1.combined <- run_seurat(cohort1)
  cohort1.combined <- RunPCA(cohort1.combined, verbose = FALSE)
  cohort1.combined <- RunUMAP(cohort1.combined, reduction = "pca", dims = 1:30)

  cohort2.combined <- run_seurat(cohort2)
  cohort2.combined <- RunPCA(cohort2.combined, verbose = FALSE)
  cohort2.combined <- RunUMAP(cohort2.combined, reduction = "pca", dims = 1:30)

  p1 <- DimPlot(cohort1.combined, reduction = "umap", group.by = "bn")
  p1

  p2 <- DimPlot(cohort2.combined, reduction = "umap", group.by = "bn")
  p2

# test.combined <- RunPCA(test.combined, verbose = FALSE)
# test.combined <- RunUMAP(test.combined, reduction = "pca", dims = 1:30)
#
# p1 <- DimPlot(test.combined, reduction = "umap", group.by = "bn")
# p2 <- DimPlot(test.combined, reduction = "umap", group.by = "seurat_annotations", label = TRUE, repel = TRUE)
# p1
# p2
