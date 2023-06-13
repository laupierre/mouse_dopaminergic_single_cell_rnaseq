## see https://dataverse.harvard.edu/dataset.xhtml?persistentId=doi:10.7910/DVN/QB5CC8

library (Seurat)

neurons <- readRDS ("neurons.rds")

# update the old seurat object
neurons <- UpdateSeuratObject(object = neurons)

head(neurons@meta.data, n=5)
table (neurons@meta.data$Curated_cellTypeLabels)
