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

neurons <-  subset(neurons, cells = WhichCells(neurons, idents = "Dopaminergic Neurons"))

table (neurons@meta.data$Curated_cellTypeLabels)
# Dopaminergic Neurons 
#                  360 

# Normalize the data
neurons <- NormalizeData(neurons, normalization.method = "LogNormalize", scale.factor = 10000)

# Identification of 5000 highly variable features 
neurons <- FindVariableFeatures(neurons, selection.method = "vst", nfeatures = 5000)

# Get the 10 most highly variable genes
top10 <- head(VariableFeatures(neurons), 10)

# plot the variable features with and without labels
plot1 <- VariableFeaturePlot(neurons)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot2

# Scaling the data of the highly variable features
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
neurons <- FindClusters(neurons)

# See the clusters (four clusters are made)
head(Idents(neurons), 5)

# Visualize in UMAP plot
neurons <- RunUMAP(neurons, dims = 1:10)
DimPlot(neurons, reduction = "umap")

# See expression of genes
# FeaturePlot will display the normalized data (from the @data slot)
FeaturePlot(neurons, features = c("Th", "Slc6a3", "Vip", "Snca"))


# Based on Th+, Slc6a3- and Vip-, to perform differential expression on cells grouped by the expression of these genes
# first create a new set of cell identities based on the expression of the gene/s

Idents(neurons, WhichCells(object = neurons, expression = Slc6a3 < 0.5 & Vip < 1, slot = 'data')) <- 'DN'
Idents(neurons, WhichCells(object = neurons, expression = Slc6a3 > 0.5, slot = 'data')) <- 'SP1'
Idents(neurons, WhichCells(object = neurons, expression = Vip > 1, slot = 'data')) <- 'SP2'
Idents(neurons, WhichCells(object = neurons, expression = Slc6a3 > 0.5 & Vip > 1, slot = 'data')) <- 'DP'

# verify the expression level
FeatureScatter(object = neurons, feature1 = "Slc6a3", feature2 = "Vip")

table (Idents (neurons))
# DP SP2 SP1  DN 
#144 108  58  50 

neurons@meta.data$mygroup <- Idents(neurons) 

# See the new identity
DimPlot(neurons, reduction = "umap")



## Finding differentially expressed features
# genes <- FindMarkers(neurons, ident.1 = 'DN', ident.2 = 'DP')

# Find markers for the DN group
# the fold change column will be named according to the logarithm base (eg, "avg_log2FC"), or if using the scale.data slot "avg_diff"
all.markers <- FindAllMarkers(object = neurons, only.pos =TRUE, return.thresh = 0.01, logfc.threshold = 0.25, min.pct = 0.25)
all.markers1 <- all.markers
head (all.markers1)


# Markers for DN population
all.markers <- all.markers[all.markers$cluster == "DN", ]
head (all.markers)
table (all.markers$p_val_adj < 0.05)
#FALSE  TRUE 
#   89     6 

res <- all.markers

# Annotation
library (org.Mm.eg.db)

#columns(org.Mm.eg.db)

symbols <- res$gene
res2 <- mapIds(org.Mm.eg.db, symbols, 'GENENAME', 'SYMBOL')

idx <- match (res$gene, names (res2))
res$Description <- as.vector (res2) [idx]
res <- res[order (res$p_val_adj), ]
res <- na.omit (res)
head (res)

library (openxlsx)

write.xlsx (res, "dopaminergic marker DAT and Vip double negative cells.xlsx", rowNames=F)


# See expression of canonical genes
FeaturePlot(neurons, features = c("Slc6a3", "Vip"), blend=TRUE, pt.size = 0.6, blend.threshold = 0.5)


# See expression of top 50 genes
library (ggpubr)

for (i in (1:50)) {
print (i)
p1 <- FeaturePlot(neurons, features = c("Slc6a3", all.markers$gene[i]), blend=TRUE, pt.size = 0.6)
p2 <- FeaturePlot(neurons, features = c("Vip", all.markers$gene[i]), blend=TRUE, pt.size = 0.6)
p3 <- ggarrange (p1, p2, nrow=2)
ggplot2::ggsave (paste (all.markers$gene[i], "dopaminergic_screen.pdf", sep="_"), p3, height=8, width=10)
}


# Heatmap of top 30 genes in each group
library (dplyr)

all.markers1 %>%
    group_by(cluster) %>%
    top_n(n = 20, wt = avg_log2FC) -> top20

DoHeatmap(neurons, features = top20$gene, slot= "scale.data")

## Warning: following features were omitted as they were not found in the scale.data slot for the RNA assay





