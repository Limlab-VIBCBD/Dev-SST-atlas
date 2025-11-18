########################
# Load required packages
########################
library(Seurat)
library(SeuratDisk)
library(dittoSeq)
library(ggplot2)
library(gtools)
library(dplyr)
library(R.utils)
library(stringr)
library(scales)
source("Seurat_Utils.R") 

###########
# Load data
###########
data <- read.csv("cell_by_gene.csv",header = TRUE, row.names = 1, colClasses = c("character", rep("numeric", ncol(read.csv("cell_by_gene.csv", nrows = 1)) - 1)))
count_matrix <- as.matrix(data)
count_matrix <- t(count_matrix)
count_matrix <- as(count_matrix, "dgCMatrix")
seurat_object <- CreateSeuratObject(counts = count_matrix)
csv_data <- read.csv("cell_metadata.csv", row.names = 1, colClasses = c("character", rep("numeric", ncol(read.csv("cell_metadata.csv", nrows = 1)) - 1)))
seurat_cell_ids <- Cells(seurat_object)
if(!all(seurat_cell_ids %in% rownames(csv_data))) {
  warning("Some cell IDs in Seurat object are missing from CSV data")
} else {
  message("All cells are present in the CSV data")
}
# Add positional information
metadata_to_add <- csv_data[seurat_cell_ids, ]
seurat_object <- AddMetaData(object = seurat_object,metadata = metadata_to_add)

#########################
# Filter and process data
#########################
Seurat_object_subset <- subset(seurat_object, subset = nCount_RNA > 0 & nFeature_RNA > 0)
Seurat_object_subset$sample <- "P5"
rm(count_matrix, csv_data, data, metadata_to_add, seurat_cell_ids, seurat_object)
gc()
Seurat_object_subset <- NormalizeData(Seurat_object_subset)
Seurat_object_subset <- FindVariableFeatures(Seurat_object_subset, nfeatures = nrow(Seurat_object_subset)) #Use all 550 genes
Seurat_object_subset <- ScaleData(Seurat_object_subset,vars.to.regress = c("nFeature_RNA"), model.use = 'linear', block.size= 5000) #Change block.size because we have ~100k cells and 550 genes

######
# PCA
######
Seurat_object_subset <- RunPCA(Seurat_object_subset, npcs = 60)
ElbowPlot(Seurat_object_subset, ndims = 60)
data.use.Seurat_object_subset <- PrepDR(object = Seurat_object_subset, genes.use = VariableFeatures(object = Seurat_object_subset), use.imputed = F, assay.type = "RNA")
nPCs.Seurat_object_subset <- PCA_estimate_nPC(data.use.Seurat_object_subset, whereto = "npc.RDS", k = 10, by.nPC = 1, from.nPC = 19, to.nPC = 23)

###########################
# Compute clusters and UMAP
###########################
Seurat_object_subset <- RunUMAP(Seurat_object_subset, dims = 1:nPCs.Seurat_object_subset)
Seurat_object_subset <- FindNeighbors(Seurat_object_subset, dims = 1:nPCs.Seurat_object_subset)
Seurat_object_subset <- FindClusters(Seurat_object_subset, resolution = 0.8)
#Save final object
saveRDS(Seurat_object_subset, file = "Merfish.rds")

############
# MapMyCells
############
counts_with_genes <- as.matrix( Merfish@assays$RNA@layers$counts)
rownames(counts_with_genes) <- rownames( Merfish@assays$RNA)
colnames(counts_with_genes) <- colnames( Merfish)
### Create h5ad file
obj<-CreateSeuratObject(counts_with_genes)
SaveH5Seurat(obj, filename = "Merfish_counts.h5Seurat", overwrite = TRUE)
Convert("Merfish_counts.h5Seurat", dest = "h5ad")
# Run MapMyCells in ABI webpage
# Import output
mapmycell<-read.csv('mapmycell.csv',comment.char="#", header = TRUE,colClasses = c("character", rep(NA, 9)))
mapmycell$cell_id <- gsub("^X", "", mapmycell$cell_id)
mapmycell$cell_id <- gsub("\\.", "-", mapmycell$cell_id)
rownames(mapmycell) <- mapmycell$cell_id
# Add MapMyCell annotation
Merfish$class_name <- mapmycell[colnames(Merfish),"class_name"] #Repeat for all the other parameters of interest
saveRDS(Merfish, file = "Merfish.rds")

################
# Label transfer
################
## Label transfer from Atlas_V2 - Major
transfer.anchors <- FindTransferAnchors(reference = Atlas_v2, query =  Merfish, dims = 1:20,reference.reduction = "pca", features=intersect(rownames(Atlas_v2), rownames(Merfish)))
predictions <- TransferData(anchorset = transfer.anchors, refdata = Atlas_v2$Major, dims = 1:20)
 Merfish$Major <- predictions$predicted.id
 Merfish$Major_score <- predictions$prediction.score.max
## Label transfer from Atlas_V2 - Minor
predictions <- TransferData(anchorset = transfer.anchors, refdata = Atlas_v2$final_clusters_renamed, dims = 1:20)
 Merfish$Minor <- predictions$predicted.id
 Merfish$Minor_score <- predictions$prediction.score.max
saveRDS(Merfish, file = " Merfish.rds")
