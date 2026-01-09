# InverseSC integration with Seurat
#
# This script shows how to use InverseSC from R via reticulate

library(Seurat)
library(reticulate)

# Import InverseSC
isc <- import("inverse_sc")

# Load your Seurat object
seu <- readRDS("your_seurat_object.rds")

# Convert to AnnData
adata <- isc$bridge$from_seurat(seu)

# Fit inverse model
isc$pp$fit_inverse_model(
  adata,
  n_epochs = 100L,
  n_programs = 20L
)

# Convert back to Seurat
seu <- isc$bridge$to_seurat(adata)

# Now seu has:
# - Z_true assay with inferred expression
# - Program_0, Program_1, ... in metadata
# - Uncertainty estimates

# Continue with standard Seurat workflow on inferred expression
DefaultAssay(seu) <- "Z_true"
seu <- FindNeighbors(seu)
seu <- FindClusters(seu)
seu <- RunUMAP(seu, dims = 1:30)

# Visualize
DimPlot(seu, reduction = "umap", group.by = "seurat_clusters")

# Check program activity
FeaturePlot(seu, features = c("Program_0", "Program_1", "Program_2"))

print("Done! Inferred expression in seu@assays$Z_true")
