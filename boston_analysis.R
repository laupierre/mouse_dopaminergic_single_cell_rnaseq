## see https://dataverse.harvard.edu/dataset.xhtml?persistentId=doi:10.7910/DVN/QB5CC8

library (Seurat)

neurons <- readRDS ("neurons.rds")

# update the old seurat 2 object from the raw.data slot
neurons5 <- CreateSeuratObject(counts = neurons@raw.data, project = "dopaminergic", min.cells = 3, min.features = 200)
neurons5@meta.data <- neurons@meta.data
neurons <- neurons5

head(neurons@meta.data, n=5)
table (neurons@meta.data$Curated_cellTypeLabels)

Idents (neurons) <- "Curated_cellTypeLabels"
length (WhichCells(neurons, idents = "Dopaminergic Neurons"))
#444

neurons <-  subset(neurons, subset = Th > 1 & percent.mito < 10)
# 389


# Normalize the data
neurons <- NormalizeData(neurons, normalization.method = "LogNormalize", scale.factor = 10000)

# Identification of 2000 highly variable features 
neurons <- FindVariableFeatures(neurons, selection.method = "vst", nfeatures = 2000)

# Get the 10 most highly variable genes
top10 <- head(VariableFeatures(neurons), 10)

# plot the variable features with and without labels
plot1 <- VariableFeaturePlot(neurons)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot2

# Scaling the data on the 2000 higlhy variable features
neurons <- ScaleData (neurons)

# perform PCA on the scaled data
neurons <- RunPCA (neurons, features = VariableFeatures(object = neurons))

VizDimLoadings(neurons, dims = 1:2, reduction = "pca")

# See how many PCA dimensions are visible (7 are possibly visible)
DimHeatmap (neurons, dims = 1:10, cells = 350, balanced = TRUE)
ElbowPlot(neurons)
# 10 dimensions is flat already

# Cluster the cells
neurons <- FindNeighbors(neurons, dims = 1:10)
neurons <- FindClusters(neurons, resolution = 0.5)







