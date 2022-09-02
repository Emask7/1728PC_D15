bn.list2 <- list(bn36181 = bn36181, bn36051 = bn36051)

find.mt <- function(x) {
  x <- PercentageFeatureSet(x, pattern = "^MT[0-9]", col.name = "percent.mt")
  x <- SCTransform(x, vars.to.regress = "percent.mt", verbose = FALSE)
  x
}

bn.list2 <- foreach(a = bn.list2) %do% (find.mt(a))
rm(a)

# Normalize datasets individually by SCTransform(), instead of NormalizeData()
features2 <- SelectIntegrationFeatures(object.list = bn.list2, nfeatures = 3000)

# Run the PrepSCTIntegration() function prior to identifying anchors
bn.list2 <- PrepSCTIntegration(object.list = bn.list2, anchor.features = features2)

# Find integration anchors and integrate the data
immune.anchors2 <- FindIntegrationAnchors(object.list = bn.list2, normalization.method = "SCT",
                                         anchor.features = features2)
save.image(file='myEnvironment.RData')

##### here #####
immune.combined.sct2 <- IntegrateData(anchorset = immune.anchors2, normalization.method = "SCT")
save.image(file='myEnvironment.RData')
