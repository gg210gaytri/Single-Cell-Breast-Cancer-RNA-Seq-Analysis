# script to perform standard workflow steps to analyze single cell RNA-Seq data  
# Of Mammalian Breast_Cancer_3p_LT_filtered_feature_bc_matrix.h5 from 10xgenomics
# load libraries
library(Seurat)
library(tidyverse)
library(hdf5r)
library(ggrepel)
library(Rtsne)


# Load the dataset and Read the sparse matrix from an HDF5 file
nsclc.sparse.m <- Read10X_h5(filename = '/home/gaytri/Downloads/Breast_Cancer_3p_LT_filtered_feature_bc_matrix.h5')
str(nsclc.sparse.m)
cts <-  nsclc.sparse.m

# Initialize the Seurat object with the raw (non-normalized data).
nsclc.seurat.obj <- CreateSeuratObject(counts = cts, project = "BreastC", min.cells = 3, min.features = 200)
str(nsclc.seurat.obj)
nsclc.seurat.obj

# View metadata
View(nsclc.seurat.obj@meta.data)

# Calculate percentage of mitochondrial reads
nsclc.seurat.obj[["percent.mt"]] <- PercentageFeatureSet(nsclc.seurat.obj, pattern = "^MT-")
View(nsclc.seurat.obj@meta.data)

# Visualize QC metrics
VlnPlot(nsclc.seurat.obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
FeatureScatter(nsclc.seurat.obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
  geom_smooth(method = 'lm')

# Filtering: Subset based on specified criteria
nsclc.seurat.obj <- subset(nsclc.seurat.obj, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & 
                             percent.mt < 5)

# Perform normalization
nsclc.seurat.obj <- NormalizeData(nsclc.seurat.obj)
str(nsclc.seurat.obj)


# Find variable features using variance stabilization transformation (vst)
nsclc.seurat.obj <- FindVariableFeatures(nsclc.seurat.obj, assay = "RNA", selection.method = "vst", nfeatures = 200, span = 0.9)
nsclc.seurat.obj[["percent.mt"]] <- PercentageFeatureSet(nsclc.seurat.obj, pattern = "^MT-")
View(nsclc.seurat.obj@meta.data)


# Identify the 20 most highly variable genes
top20 <- head(VariableFeatures(nsclc.seurat.obj), 20)
plot1 <- VariableFeaturePlot(nsclc.seurat.obj)
LabelPoints(plot = plot1, points = top20, repel = TRUE, xnudge = 0, ynudge = 0, max.overlaps = 100)

``
unlabeled_points <- top20[!top20 %in% LabelPoints(plot = plot1, points = top20, repel = TRUE, xnudge = 0, ynudge = 0, max.overlaps = 100)$label]
print(unlabeled_points)

any(is.infinite(top20))
print(top20)

# Scaling
all.genes <- rownames(nsclc.seurat.obj)
nsclc.seurat.obj <- ScaleData(nsclc.seurat.obj, features = all.genes)
str(nsclc.seurat.obj)

# Run Principal Component Analysis (PCA)
nsclc.seurat.obj <- RunPCA(nsclc.seurat.obj, features = VariableFeatures(object = nsclc.seurat.obj), npcs = 10)

# Visualize PCA results
print(nsclc.seurat.obj[["pca"]], dims = 1:5, nfeatures = 5)
DimHeatmap(nsclc.seurat.obj, dims = 1, cells = 22, balanced = TRUE)

# Determine dimensionality of the data # Plot elbow plot
ElbowPlot(nsclc.seurat.obj)

# Clustering # Find neighbors for clustering
nsclc.seurat.obj <- FindNeighbors(nsclc.seurat.obj, dims = 1:10)

# # Find clusters with specified resolution
nsclc.seurat.obj <- FindClusters(nsclc.seurat.obj, resolution = 0.5)  # Choose an appropriate resolution
View(nsclc.seurat.obj@meta.data)

# Visualize clusters
DimPlot(nsclc.seurat.obj, group.by = "RNA_snn_res.0.5", label = TRUE)

# Setting identity of clusters
Idents(nsclc.seurat.obj)
Idents(nsclc.seurat.obj) <- "RNA_snn_res.0.1"  # Adjust to the chosen resolution
Idents(nsclc.seurat.obj)

# Display available assays or layers
print(names(nsclc.seurat.obj@assays))


# Non-linear dimensionality reduction with t-SNE
# Extract expression matrix
expression_matrix <- GetAssayData(object = nsclc.seurat.obj, assay = "RNA")
# Convert expression matrix to matrix
expression_matrix <- as.matrix(expression_matrix)
# Remove duplicate rows from the expression matrix
expression_matrix <- expression_matrix[!duplicated(expression_matrix), ]
# Non-linear dimensionality reduction with t-SNE
nsclc.seurat.obj <- RunTSNE(object = nsclc.seurat.obj, dims = 1:10, perplexity = 7)
# Plot the t-SNE results
DimPlot(nsclc.seurat.obj, reduction = "tsne")
