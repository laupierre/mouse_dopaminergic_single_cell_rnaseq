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

neurons <-  subset(neurons, subset = Th > 1)
# 389

# Normalize the data
neurons <- NormalizeData(neurons, normalization.method = "LogNormalize", scale.factor = 10000)
