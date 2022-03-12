library(monocle3)
library(Seurat)

#Load dataset
load("C:/Users/Administrator/Desktop/Organoid-new6.RData")
DefaultAssay(Organoid.combined) <- "RNA"

#Create cds and preprocess
data <- GetAssayData(Organoid.combined, assay = 'RNA', slot = 'counts')
cell_metadata <- Organoid.combined@meta.data
cell_metadata$celltype <- Organoid.combined@active.ident
gene_annotation <- data.frame(gene_short_name = rownames(data))
rownames(gene_annotation) <- rownames(data)
cds <- new_cell_data_set(data,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_annotation)
rm(cell_metadata)
rm(data)
rm(gene_annotation)
cds <- preprocess_cds(cds, num_dim = 18)
cds <- align_cds(cds = cds,alignment_group = "orig.ident")
cds <- reduce_dimension(cds, preprocess_method = "PCA")
cds.embed <- cds@int_colData$reducedDims$UMAP
int.embed <- Embeddings(Organoid.combined, reduction = "umap")
rm(Organoid.combined)
int.embed <- int.embed[rownames(cds.embed),]
cds@int_colData$reducedDims$UMAP <- int.embed

cds <- cluster_cells(cds)
cds <- learn_graph(cds)
cds <- order_cells(cds)

