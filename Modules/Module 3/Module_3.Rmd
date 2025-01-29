---
title: "Atlas_integration_and_filtering (Module 3)"
output: .txt, .rds and .pdf documents
---

## Required installation
```{r}
if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")

pacman::p_load(
  Seurat, dittoSeq, ggplot2, gtools, dplyr, R.utils, colorspace, tidyr,
  purrr, Matrix, patchwork, BiocManager
)

devtools::install_github("crazyhottommy/scclusteval")
BiocManager::install("scDblFinder")
BiocManager::install("MAST")
BiocManager::install("dittoSeq")
```

## Load libraries
```{r}
library(Seurat)
library(dittoSeq)
library(ggplot2)
library(gtools)
library(dplyr)
library(dittoSeq)
library(scDblFinder)
library(R.utils)
library(colorspace)
library(scclusteval)
library(tidyr)
library(purrr)
library(Matrix)
library(patchwork)
library(MAST)
source("Seurat_Utils.R")
colors_ditto<-dittoColors()
names(colors_ditto)<-as.character(c(0:(length(colors_ditto)-1)))
```

### LOAD FILTERED SST+ CELLS FROM ALL SAMPLES
```{r}
# Load previously saved Seurat objects
base_atlas <- readRDS("integrated_all_to_send.rds")
E16_Lim4 <- readRDS("lim_E16_may24_sst_seu_sub_clean_labeled.rds")
E16_Lim4$sample <- "E16_Lim4"
P1_Lim3 <- readRDS('lim_P1_old_seu_sub_clean_no8.rds')
P1_Lim3$sample <- "P1_Lim3"
#P1_Lim5 <- readRDS("lim_P1_march24_SST_seu_sub_clean2_labeled.rds")
#P1_Lim5$sample <- "P1_Lim5"
#test_dataset_filtered_sst<-readRDS("WT1_WT23_clean_labeled_minor_clusters.rds")
#P5_WT1_Lim1 <- subset(test_dataset_filtered_sst, subset= orig.ident %in% "lim_P5_sorted_sept23_WT1")
#P5_WT1_Lim1$sample <- "P5_WT1_Lim1"
#test_dataset_filtered_sst<-readRDS("WT1_WT23_clean_labeled_minor_clusters.rds")
#P5_WT23_Lim2 <- subset(test_dataset_filtered_sst, subset= orig.ident %in% "lim_P5_sorted_sept23_WT23")
#P5_WT23_Lim2$sample <- "P5_WT23_Lim2"
#E18_Lippi<-readRDS("Lippi_lab_e18.5_seu_sst_final.rds")
#E18_Lippi$sample <- "E18_Lippi"
#P1_EMI014_Lim<-readRDS("lim_P1_EMI014_finalFILTERED_withDoublets_upto_clustering.rds")
#P1_EMI014_Lim$sample <- "P1_EMI014_Lim"
base_atlas$sample <- paste0("BaseAtlas_",base_atlas$orig.ident)

# Set the default assay to RNA for analysis
DefaultAssay(base_atlas) <- "RNA"
DefaultAssay(E16_Lim4) <- "RNA"
DefaultAssay(P1_Lim3) <- "RNA"
#DefaultAssay(P1_Lim5) <- "RNA"
#DefaultAssay(P5_WT1_Lim1) <- "RNA"
#DefaultAssay(P5_WT23_Lim2) <- "RNA"
#DefaultAssay(E18_Lippi) <- "RNA"
#DefaultAssay(P1_EMI014_Lim) <- "RNA"

# Split the base atlas into separate objects based on the sample type
integrated_all_list <- SplitObject(base_atlas, split.by = "sample")

# Add the cleaned dataset to the list of samples
integrated_all_list$E16_Lim4 <- E16_Lim4
integrated_all_list$P1_Lim3 <- P1_Lim3
#integrated_all_list$P1_Lim5 <- P1_Lim5
#integrated_all_list$P5_WT1_Lim1 <- P5_WT1_Lim1
#integrated_all_list$P5_WT23_Lim2 <- P5_WT23_Lim2
#integrated_all_list$P1_EMI014_Lim <- P1_EMI014_Lim
#integrated_all_list$E18_Lippi <- E18_Lippi
```

### Minor Label Annotation with Label Transfer Using Base Atlas as Reference
```{r}
# Normalize, find variable features, and scale data for the base atlas
base_atlas <- NormalizeData(base_atlas)
base_atlas <- FindVariableFeatures(base_atlas)
base_atlas <- ScaleData(base_atlas, vars.to.regress = c("nFeature_RNA",'percent.mt','ccDiff'))
base_atlas <- RunPCA(base_atlas,npcs = 50)

# Initialize empty vector to store minor labels
minor_label <- c()

# Perform label transfer
for (sample in c("E16_Lim4","P1_Lim3")) {
  integrated_all_list[[sample]] <- NormalizeData(integrated_all_list[[sample]])
  integrated_all_list[[sample]] <- FindVariableFeatures(integrated_all_list[[sample]])
  integrated_all_list[[sample]] <- ScaleData(integrated_all_list[[sample]], vars.to.regress = c("nFeature_RNA",'percent.mt','ccDiff'))
  integrated_all_list[[sample]] <- RunPCA(integrated_all_list[[sample]],npcs = 50)
  # Find transfer anchors between the base atlas (reference) and the samples
  transfer.anchors <- FindTransferAnchors(reference = base_atlas, query = integrated_all_list[[sample]], dims = 1:40, reference.reduction = "pca", features=intersect(rownames(base_atlas), rownames(integrated_all_list[[sample]])))
  # Transfer labels based on the reference
  predictions <- TransferData(anchorset = transfer.anchors, refdata = base_atlas$cluster_label_trained_with_all, dims = 1:40)
  # Add the predicted minor labels to the sample object
  integrated_all_list[[sample]]$minor_label_transferAnchors_BaseAtlas <- predictions$predicted.id
  integrated_all_list[[sample]]@meta.data[rownames(predictions),"minor_label_transferAnchors_BaseAtlas"] <- predictions$predicted.id
  # Append the predicted labels to the overall list of minor labels
  minor_pred <- predictions$predicted.id
  names(minor_pred) <- rownames(predictions)
  minor_label <- c(minor_label,minor_pred)
}
#saveRDS(minor_label,"Minor_label_labelTranfer_beforeDoubletsRemoval.RDS")

# Assign the reference-trained labels to the other samples
for (sample in c("BaseAtlas_E16","BaseAtlas_lim_P5_fixed_sorted","BaseAtlas_P1","BaseAtlas_P5") ) {
  integrated_all_list[[sample]]$minor_label_transferAnchors_BaseAtlas <- integrated_all_list[[sample]]$cluster_label_trained_with_all
}

rm(list = c("predictions", "transfer.anchors", "minor_label", "minor_pred", "E16_Lim4", "P1_Lim3"))
gc()
```

### Doublets Prediction Using Minor Labels
```{r}
# Initialize an empty vector to store doublet annotations
doublets_sceDblF <- c()

# Loop through all datasets in the integrated list and predict doublets
for (i in 1:length(integrated_all_list)){
  sceDblF <- scDblFinder(integrated_all_list[[i]]@assays$RNA@counts,dbr =0.07, clusters=integrated_all_list[[i]]$minor_label_transferAnchors_BaseAtlas)
  # Extract doublet annotations from the results
  doublets_anno <- as.vector(sceDblF@colData$scDblFinder.class)
  names(doublets_anno) <- row.names(sceDblF@colData)
  # Append the doublet annotations to the list
  doublets_sceDblF <- c(doublets_sceDblF,doublets_anno)
  # Add the doublet classification to the Seurat object metadata
  integrated_all_list[[i]]$doublets <- doublets_anno
  integrated_all_list[[i]]@meta.data[names(doublets_anno),"doublets"] <- doublets_anno
  # Subset the Seurat object to retain only the singlets
  integrated_all_list[[i]] <- subset(integrated_all_list[[i]], subset=doublets=="singlet")
}
#saveRDS(doublets_sceDblF,"doublets_sceDblFinder_allSamples.RDS")


rm(list = c("sceDblF", "doublets_sceDblF", "doublets_anno"))
gc()
```


### Perform Integration of Datasets
```{r}
# Normalize data, find variable features, and process each dataset
for(sample in names(integrated_all_list)){
  obj <- integrated_all_list[[sample]]
  obj <- NormalizeData(obj)
  obj <- FindVariableFeatures(obj)
  integrated_all_list[[sample]] <- obj
}

# Find integration anchors between datasets
i.anchors <- FindIntegrationAnchors(object.list = integrated_all_list, dims = 1:30, reduction = 'cca', scale = T, k.anchor = 5, k.filter = 100, k.score = 15, anchor.features = 3000)

# Integrate the data using the anchors
integrated_v2 <- IntegrateData(anchorset = i.anchors, dims = 1:30, normalization.method ='LogNormalize', k.weight=100)

DefaultAssay(integrated_v2) <- 'integrated'
integrated_v2 <- ScaleData(integrated_v2, verbose = FALSE, vars.to.regress = c("nFeature_RNA",'percent.mt','ccDiff'))
integrated_v2 <- RunPCA(integrated_v2, npcs = 100)
ElbowPlot(integrated_v2, ndims = 60)
#saveRDS(integrated_v2, file="Integrated_atlas_V2.RDS")

data.use.integrated <- PrepDR(object = integrated_v2, genes.use = VariableFeatures(object = integrated_v2), use.imputed = F, assay.type = "integrated")
path_data <- getwd()

# Check if the optimal nPCs files exist
if (!file.exists("optimal_nPCs_5_integrated.RDS")) { 
  # If it doesn't exist, calculate and save the nPCs for the integrated dataset
  nPCs.data.use5 <- PCA_estimate_nPC(data.use.integrated, 
                                     whereto = paste0(path_data, "/optimal_nPCs_5_integrated.RDS"), 
                                     k = 5, by.nPC = 5, from.nPC = 30, to.nPC = 50) #check ElbowPlot and adjust range acoordingly

} else {
  # Load precomputed nPCs if they exist
  nPCs.data.use5 <- readRDS("optimal_nPCs_5_integrated.RDS")
}

if (!file.exists("optimal_nPCs_integrated.RDS")) { 
  # If it doesn't exist, calculate and save the nPCs for the full integrated dataset
  nPCs.data.use <- PCA_estimate_nPC(data.use.integrated, 
                                    whereto = paste0(path_data, "/optimal_nPCs_integrated.RDS"), 
                                    k = 5, by.nPC = 1, from.nPC = nPCs.data.use5 - 5, to.nPC = nPCs.data.use5 + 5)
} else {
  # Load precomputed nPCs if they exist
  nPCs.data.use <- readRDS("optimal_nPCs_integrated.RDS")
}

integrated_v2 <- RunUMAP(integrated_v2, dims=1:nPCs.data.use)
integrated_v2 <- FindNeighbors(integrated_v2, dims=1:nPCs.data.use)
integrated_v2 <- FindClusters(integrated_v2, resolution = c(0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.5))
#saveRDS(integrated_v2, file= "Integrated_atlas_V2_cluster.RDS")


rm(list = c("i.anchors", "data.use.integrated", "integrated_all_list", "obj", "PrepDR", "PCA_estimate_nPC", "nPCs.data.use5", "nPCs.data.use"))
gc()
```

### PERFORM QC AND REMOVE STRESSED CELLS #########
```{r}
DefaultAssay(integrated_v2) <- 'RNA'
integrated_v2 <- ScaleData(integrated_v2)

# Extract UMAP embeddings (2D coordinates for visualization)
umap <- Embeddings(integrated_v2,reduction = "umap")

# Loop through a list of stress-related genes to visualize their expression on UMAP and export in pdf
for (gene in c("Hsp90b1","Hspa5","Rpn1","Mapk8","Rpl27a","Rpl13","Rpl15")){
  # Get the gene expression data
  data <- integrated_v2@assays$RNA@data[gene,]
  # Combine the UMAP coordinates with gene expression values
  umap_gene <- cbind(umap[colnames((integrated_v2@assays$RNA@data)),1:2], data)
  # Order the data by expression values (ascending order)
  umap_gene <- umap_gene[order(umap_gene[,3], decreasing = FALSE),]
  colnames(umap_gene) = c("umap_1","umap_2","Expression")
  umap_gene <- as.data.frame(umap_gene)
  pdf(paste("GeneExpression_",gene,"_umap.pdf",sep=""))
  #print(ggplot(umap_gene, aes(umap_1, umap_2)) +  geom_point(aes(colour = Expression), size=1) + scale_color_continuous_sequential(palette='Purple_Yellow') +     theme(panel.background = element_rect(fill='white', colour='black')) +  theme(legend.position="none"))
  print(ggplot(umap_gene, aes(umap_1, umap_2)) +  geom_point(aes(colour = Expression), size=1) + scale_color_continuous_sequential(palette='Purple_Yellow') + theme(panel.background = element_rect(fill='white', colour='black'))+ ggtitle( paste0(gene)) ) 
  dev.off()
}

# Create UMAP plots based on feature count (nFeature) and mitochondrial percentage (percent.mt)
dittoDimPlot(integrated_v2, "nFeature_RNA", reduction.use = "umap", min.color = "lightgrey", max.color = "blue")
dittoDimPlot(integrated_v2, "percent.mt", reduction.use = "umap", min.color = "lightgrey", max.color = "blue")


# Create UMAP plots for different clustering resolutions and export in pdf
pdf(paste0('umap_plot_atlas_v2_clusters.pdf'))
for(j in c(0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.5)){
  Idents(integrated_v2) <- paste0("integrated_snn_res.",j)
  print(DimPlot(integrated_v2, reduction = "umap",raster=FALSE, label=TRUE)  + ggtitle(paste0("Resolution: ",j)))
}  
dev.off()

# Set custom color scheme for minor labels (clusters). Use the names of your label transfer
integrated_v2$minor_label_transferAnchors_BaseAtlas <- factor(integrated_v2$minor_label_transferAnchors_BaseAtlas, levels = c("LRP_1","LRP_2","LRP3","LRP4","Martinotti_1","Martinotti_2","Martinotti_3","Non-Martinotti_1","Non-Martinotti_2","Non-Martinotti_3","Non-Martinotti_4","Non-Martinotti_5"))
ClusterCol <- c('#1E2D59', '#0C85B9', '#7CD5E2', '#4BC3AB','#61257A', '#C69BDF', '#D43DA7', '#C73A3A', '#F78173', '#FFB870', '#B59734', '#ECDD77')
names(ClusterCol) <- c("LRP_1","LRP_2","LRP3","LRP4","Martinotti_1","Martinotti_2","Martinotti_3","Non-Martinotti_1","Non-Martinotti_2","Non-Martinotti_3","Non-Martinotti_4","Non-Martinotti_5")

# Plot the minor label clusters with custom color schemes
dittoDimPlot(integrated_v2, "minor_label_transferAnchors_BaseAtlas", split.by = "sample", color.panel=ClusterCol)

# Plot the minor label clusters (without splitting by sample)
dittoDimPlot(integrated_v2, "minor_label_transferAnchors_BaseAtlas", color.panel=ClusterCol)

# Set the clustering identity to the specific resolution (0.4) and generate bar plots of samples per cluster
Idents(integrated_v2)<-'integrated_snn_res.0.4'
dittoBarPlot(integrated_v2, var = "sample", group.by = "integrated_snn_res.0.4")  + labs(title = NULL)

# Generate violin plots for expression of various genes by cluster (based on the integrated SNN resolution 0.4)
VlnPlot(integrated_v2, "Rpn1", pt.size = 0, group.by = "integrated_snn_res.0.4")
VlnPlot(integrated_v2, "Rpl15", pt.size = 0, group.by = "integrated_snn_res.0.4")
VlnPlot(integrated_v2, "Rpl17", pt.size = 0, group.by = "integrated_snn_res.0.4")
VlnPlot(integrated_v2, "Hsp90b1", pt.size = 0, group.by = "integrated_snn_res.0.4")
VlnPlot(integrated_v2, "Hspa5", pt.size = 0, group.by = "integrated_snn_res.0.4")
VlnPlot(integrated_v2, "Mapk8", pt.size = 0, group.by = "integrated_snn_res.0.4", cols = dittoColors())


# Add a new feature for ribosomal RNA percentage (ribo genes)
integrated_v2[["percent.ribo"]] <- PercentageFeatureSet(integrated_v2, pattern = "^Rp[Sl]") 

# Violin plot for ribosomal RNA percentage by cluster
VlnPlot(integrated_v2, "percent.ribo", pt.size = 0, group.by = "integrated_snn_res.0.4")

# Violin plot for mitochondrial RNA percentage by cluster
VlnPlot(integrated_v2, "percent.mt", pt.size = 0, group.by = "integrated_snn_res.0.4")

# Violin plots for total RNA count and feature count by cluster
VlnPlot(integrated_v2, "nCount_RNA", pt.size = 0, group.by = "integrated_snn_res.0.4") + ylim(0,30000)
VlnPlot(integrated_v2, "nFeature_RNA", pt.size = 0, group.by = "integrated_snn_res.0.4")

# Find marker genes for all clusters using MAST (statistical test) and save results
Idents(integrated_v2)<-'integrated_snn_res.0.4'
MarkersAll <- FindAllMarkers(integrated_v2,only.pos = TRUE,logfc.threshold = 1,min.pct = 0.4, test.use = 'MAST', latent.vars = 'sample')
#saveRDS(MarkersAll,"MarkersAll_res0.4_integrated_MAST.RDS")

# Generate a dot plot for marker genes across clusters
Idents(integrated_v2) <- "integrated_snn_res.0.4"
pdf('DotPlot_MarkerGenes_clusters_res0.4.pdf', width=17)
par(las=2)
print(DotPlot(integrated_v2, assay="RNA", features=unique(MarkersAll$gene),scale=TRUE) + ylab("Clusters") + scale_colour_gradient2(low = "blue2", mid = "gray90", high = "red2", midpoint=0) + theme(axis.text.x = element_text(angle = 90, hjust = 1)))



rm(list = c("MarkersAll", "integrated_v2", "umap", "umap_gene"))
gc()
```

#* ATLAS V2 on filtered cells (remove all filtered cells from time-point + cells from clusters 8,14,15 (res 0.4) from atlas first integration)
### Perform Integration again of filtered cells
```{r}
# Specify clusters you want to remove
#integrated_v2_filtered<-subset(integrated_v2, idents = c('8'), invert = T)
# - or -  Specify clusters you want to keep
integrated_v2_filtered<-subset(integrated_v2, idents = c('5'), invert = F)

DefaultAssay(integrated_v2_filtered) <- 'RNA'
integrated_all_list <- SplitObject(integrated_v2_filtered, split.by = "sample")

# Loop over each sample and normalize the data, then find variable features
for(sample in names(integrated_all_list)){
  
  # Normalize the data for each sample
  obj <- integrated_all_list[[sample]]
  obj <- NormalizeData(obj)
  obj <- FindVariableFeatures(obj)
  # Save the updated sample object back to the list
  integrated_all_list[[sample]] <- obj
}

# Find integration anchors between the different samples (using CCA for dimensionality reduction)
i.anchors <- FindIntegrationAnchors(object.list = integrated_all_list, dims = 1:30, reduction='cca',scale=T, k.anchor=5, k.filter=100, k.score=15,anchor.features=3000)
integrated_v2_filtered <- IntegrateData(anchorset = i.anchors, dims = 1:30, normalization.method = 'LogNormalize', k.weight=80)
DefaultAssay(integrated_v2_filtered) <- 'integrated'

# Identify mitochondrial genes (genes with prefix 'mt-') and ribosomal genes (genes with prefix 'Rp[Sl]')
mt_genes <- grep(pattern = "^mt-", x = rownames(integrated_v2_filtered@assays$integrated@data), value = TRUE)
ribo_genes <- grep(pattern = "^Rp[Sl]", x = rownames(integrated_v2_filtered@assays$integrated@data), value = TRUE)

integrated_v2_filtered <- ScaleData(integrated_v2_filtered, verbose = FALSE, vars.to.regress = c("nFeature_RNA",'percent.mt','ccDiff'))
integrated_v2_filtered <-RunPCA(integrated_v2_filtered, npcs = 100)
ElbowPlot(integrated_v2, ndims = 60)
#saveRDS(integrated_v2_filtered, file = "Integrated_atlas_V2_filtered.RDS")
data.use.integrated<- PrepDR(object = integrated_v2_filtered, genes.use = VariableFeatures(object = integrated_v2_filtered), use.imputed = F, assay.type = "integrated")
path_data <- getwd()

if (!file.exists(paste0(path_data, "/optimal_nPCs_5_integrated_filtered.RDS"))) {
  # If the file doesn't exist, calculate and save the nPCs for the integrated dataset
  nPCs.data.use5 <- PCA_estimate_nPC(data.use.integrated, 
                                     whereto = paste0(path_data, "/optimal_nPCs_5_integrated_filtered.RDS"), 
                                     k = 5, by.nPC = 5, from.nPC = 40, to.nPC = 60) # Check Elbow plot and set the range for optimal  PC

} else {
  # If the file exists, load the precomputed nPCs
  nPCs.data.use5 <- readRDS(paste0(path_data, "/optimal_nPCs_5_integrated_filtered.RDS"))
}

# Estimate the optimal number of PCs again, focusing on a smaller range around the previously estimated value
if (!file.exists(paste0(path_data, "/optimal_nPCs_integrated_filtered.RDS"))) { # Check if the second nPCs file exists for nPCs.data.use

  # If the file doesn't exist, calculate and save the nPCs for the integrated dataset
  nPCs.data.use <- PCA_estimate_nPC(data.use.integrated, 
                                    whereto = paste0(path_data, "/optimal_nPCs_integrated_filtered.RDS"), 
                                    k = 2, by.nPC = 1, from.nPC = nPCs.data.use5 - 3, to.nPC = nPCs.data.use5 + 3)
} else {
  # If the file exists, load the precomputed nPCs
  nPCs.data.use <- readRDS(paste0(path_data, "/optimal_nPCs_integrated_filtered.RDS"))
}

integrated_v2_filtered <- RunUMAP(integrated_v2_filtered, dims=1:nPCs.data.use)
integrated_v2_filtered <- FindNeighbors(integrated_v2_filtered, dims=1:nPCs.data.use)
integrated_v2_filtered <- FindClusters(integrated_v2_filtered, resolution = c(0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.5))

#saveRDS(integrated_v2_filtered, file = "Integrated_atlas_V2_filtered_clusters.RDS")

rm(list = c("i.anchors", "data.use.integrated", "obj", "PrepDR", "PCA_estimate_nPC", "integrated_all_list", "integrated_v2"))
gc()
```

### Perform QC #########
```{r}
DefaultAssay(integrated_v2_filtered) <- 'RNA'
integrated_v2_filtered <- ScaleData(integrated_v2_filtered)

# Extract the UMAP embeddings (coordinates) for visualization
umap<-Embeddings(integrated_v2_filtered,reduction = "umap")

# Iterate over a predefined list of genes to plot their expression on UMAP
for (gene in c("Hsp90b1","Hspa5","Rpn1","Mapk8","Rpl27a","Rpl13","Rpl15")){
  # Extract the expression data for the current gene
  data <- integrated_v2_filtered@assays$RNA@data[gene,]
  
  # Combine the UMAP coordinates with the gene expression data
  umap_gene <- cbind(umap[colnames((integrated_v2_filtered@assays$RNA@data)),1:2], data)
  
  # Sort the data by gene expression values (ascending order)
  umap_gene <- umap_gene[order(umap_gene[,3], decreasing=FALSE),]
  colnames(umap_gene) = c("umap_1","umap_2","Expression")
  umap_gene <- as.data.frame(umap_gene)
  #pdf(paste("GeneExpression_",gene,"_umap_filtered.pdf",sep=""))
  
  # Create a UMAP plot with gene expression overlaid
  #print(ggplot(umap_gene, aes(umap_1, umap_2)) +  geom_point(aes(colour = Expression), size=1) + scale_color_continuous_sequential(palette='Purple_Yellow') + theme(panel.background = element_rect(fill='white', colour='black')) +  theme(legend.position="none"))
  print(ggplot(umap_gene, aes(umap_1, umap_2)) +  geom_point(aes(colour = Expression), size=1) + scale_color_continuous_sequential(palette='Purple_Yellow') + theme(panel.background = element_rect(fill='white', colour='black'))+ ggtitle( paste0(gene)) ) 
  dev.off()
}

# Create a UMAP plot showing the distribution of the 'nFeature_RNA' (number of detected features per cell)
dittoDimPlot(integrated_v2_filtered, "nFeature_RNA", reduction.use = "umap", min.color = "lightgrey", max.color = "blue")

# Create UMAP plots for different clustering resolutions
pdf(paste0('umap_plot_atlas_v2_clusters_filtered.pdf'))
for(j in c(0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.5)){
  Idents(integrated_v2_filtered) <- paste0("integrated_snn_res.",j)
  print(DimPlot(integrated_v2_filtered, reduction = "umap",raster=FALSE, label=TRUE)  + ggtitle(paste0("Resolution: ",j)))
}  
dev.off()

#Define new clutering resolution
Idents(integrated_v2_filtered) <- 'integrated_snn_res.0.7'

# Create a bar plot showing the number of cells per sample, grouped by clusters at resolution 0.7
dittoBarPlot(integrated_v2_filtered, var = "sample", group.by = "integrated_snn_res.0.7")  + labs(title = NULL)

# Generate violin plots for the expression of specific genes across clusters at resolution 0.7
VlnPlot(integrated_v2_filtered, "Rpn1", pt.size = 0, group.by = "integrated_snn_res.0.7", cols = dittoColors())
VlnPlot(integrated_v2_filtered, "Hsp90b1", pt.size = 0, group.by = "integrated_snn_res.0.7", cols = dittoColors())
VlnPlot(integrated_v2_filtered, "Hspa5", pt.size = 0, group.by = "integrated_snn_res.0.7", cols = dittoColors())
VlnPlot(integrated_v2_filtered, "Mapk8", pt.size = 0, group.by = "integrated_snn_res.0.7", cols = dittoColors())

# Calculate the percentage of ribosomal genes
integrated_v2_filtered[["percent.ribo"]] <- PercentageFeatureSet(integrated_v2_filtered, pattern = "^Rp[Sl]") ## NB not all samples have ribo genes: for these samples perc.ribo will be 0

# Generate violin plots for the percentage of ribosomal genes across samples, grouped by clusters at resolution 0.7
VlnPlot(integrated_v2_filtered, "percent.ribo", pt.size = 0, group.by = "integrated_snn_res.0.7")

# Generate violin plots for the percentage of mitochondrial genes (those starting with 'mt-') across samples
VlnPlot(integrated_v2_filtered, "percent.mt", pt.size = 0, group.by = "integrated_snn_res.0.7")

# Generate violin plots for the total RNA count across samples and clusters
VlnPlot(integrated_v2_filtered, "nCount_RNA", pt.size = 0, group.by = "integrated_snn_res.0.7") + ylim(0,30000)

# Generate violin plots for the number of features detected per cell across clusters
VlnPlot(integrated_v2_filtered, "nFeature_RNA", pt.size = 0, group.by = "integrated_snn_res.0.7")


rm(list = c("umap", "umap_gene", "ClusterCol", "colors_ditto", "data", "i", "j", "mt_genes", "gene", "path_data", "ribo_genes", "sample", "nPCs.data.use", "nPCs.data.use5"))
gc()
```

### Re-annotate filtered dataset from Base_Atlas######
```{r}
DefaultAssay(integrated_v2_filtered) <- 'RNA'

integrated_all_list <- SplitObject(integrated_v2_filtered, split.by = "sample")

base_atlas <- subset(integrated_v2_filtered, subset= sample %in% c("BaseAtlas_E16","BaseAtlas_lim_P5_fixed_sorted","BaseAtlas_P1","BaseAtlas_P5"))
DefaultAssay(base_atlas) <- 'RNA'
base_atlas <- NormalizeData(base_atlas)
base_atlas <- FindVariableFeatures(base_atlas)
base_atlas <- ScaleData(base_atlas, vars.to.regress = c("nFeature_RNA",'percent.mt','ccDiff'))
base_atlas<-RunPCA(base_atlas,npcs = 50)

# Initialize empty vectors to store predicted labels and scores
minor_label <- c()
major_label <- c()
pred_score  <- c()

# Loop through each sample, normalize, identify variable features, scale data, and run PCA
for (sample in c("E16_Lim4","P1_Lim3")) {
  
  # Normalize, find variable features, and scale the data for each sample
  integrated_all_list[[sample]] <- NormalizeData(integrated_all_list[[sample]])
  integrated_all_list[[sample]] <- FindVariableFeatures(integrated_all_list[[sample]])
  integrated_all_list[[sample]] <- ScaleData(integrated_all_list[[sample]], vars.to.regress = c("nFeature_RNA",'percent.mt','ccDiff'))
  integrated_all_list[[sample]] <- RunPCA(integrated_all_list[[sample]],npcs = 50)
  
  # Perform anchor finding for label transfer based on PCA reduction and intersected features
  transfer.anchors <- FindTransferAnchors(reference = base_atlas, query = integrated_all_list[[sample]], dims = 1:40,reference.reduction = "pca", features=intersect(rownames(base_atlas), rownames(integrated_all_list[[sample]])))
  
  # Transfer the predicted minor labels
  predictions <- TransferData(anchorset = transfer.anchors, refdata = base_atlas$cluster_label_trained_with_all, dims = 1:40)
  minor_pred <- predictions$predicted.id
  names(minor_pred) <- rownames(predictions)
  minor_label <- c(minor_label,minor_pred)
  
  # Capture the prediction score for each minor label prediction
  pred_score_tmp <- predictions$prediction.score.max
  names(pred_score_tmp) <- rownames(predictions)
  pred_score <- c(pred_score,pred_score_tmp)
  
  # Transfer the predicted major labels
  predictions <- TransferData(anchorset = transfer.anchors, refdata = base_atlas$major_cluster_label_trained_with_all, dims = 1:40)
  major_pred <- predictions$predicted.id
  names(major_pred) <- rownames(predictions)
  major_label <- c(major_label,major_pred)
}

#saveRDS(minor_label,"Minor_label_labelTranfer_afterFiltering.RDS")
#saveRDS(major_label,"Major_label_labelTranfer_afterFiltering.RDS")
#saveRDS(pred_score,"pred_score_minor_label.RDS")

# Assign minor and major labels to the meta-data of the integrated filtered dataset
integrated_v2_filtered$minor_label_transferAnchors_BaseAtlas <- integrated_v2_filtered$cluster_label_trained_with_all
integrated_v2_filtered@meta.data[names(minor_label),'minor_label_transferAnchors_BaseAtlas'] <- minor_label
integrated_v2_filtered$major_label_transferAnchors_BaseAtlas <- integrated_v2_filtered$major_cluster_label_trained_with_all
integrated_v2_filtered@meta.data[names(major_label),'major_label_transferAnchors_BaseAtlas'] <- major_label
saveRDS(integrated_v2_filtered, "Integrated_atlas_V2_filtered_annotated.RDS")

# Factorize the minor labels and create a custom color palette
integrated_v2_filtered$minor_label_transferAnchors_BaseAtlas <- factor(integrated_v2_filtered$minor_label_transferAnchors_BaseAtlas, levels = c("LRP_1","LRP_2","LRP3","LRP4","Martinotti_1","Martinotti_2","Martinotti_3","Non-Martinotti_1","Non-Martinotti_2","Non-Martinotti_3","Non-Martinotti_4","Non-Martinotti_5"))
ClusterCol <- c('#1E2D59', '#0C85B9', '#7CD5E2', '#4BC3AB','#61257A', '#C69BDF', '#D43DA7', '#C73A3A', '#F78173', '#FFB870', '#B59734', '#ECDD77')
names(ClusterCol) <- c("LRP_1","LRP_2","LRP3","LRP4","Martinotti_1","Martinotti_2","Martinotti_3","Non-Martinotti_1","Non-Martinotti_2","Non-Martinotti_3","Non-Martinotti_4","Non-Martinotti_5")

# Generate UMAP plots for minor labels
dittoDimPlot(integrated_v2_filtered,"minor_label_transferAnchors_BaseAtlas", split.by = "sample", color.panel=ClusterCol)
dittoDimPlot(integrated_v2_filtered,"minor_label_transferAnchors_BaseAtlas", color.panel=ClusterCol)

# Generate UMAP plots for major labels
ClusterCol <- c('#35C1D5', '#4A2884',  '#E05F36')
dittoDimPlot(integrated_v2_filtered,"major_label_transferAnchors_BaseAtlas", split.by = "sample", color.panel=ClusterCol)
DimPlot(integrated_v2_filtered, group.by="major_label_transferAnchors_BaseAtlas",  reduction = 'umap', cols=ClusterCol)
#saveRDS(integrated_v2_filtered,"Atlas_v2.RDS")

rm(list = c("base_atlas", "integrated_all_list", "integrated_v2_filtered"))
gc()
```

#* PERFORM INTEGRATION OF MAJOR SUBSETS (LRP, MC, nMC) #########
```{r}
setwd('/lustre/projects/lim_core/results/tgi/Atlas_v2/integration_filtered_atlas/')
integrated_v2_filtered <- readRDS("Integrated_atlas_V2_filtered_annotated.RDS")
DefaultAssay(integrated_v2_filtered) <- 'RNA'

# Subset the LRP group based on the major label and create a new directory for LRP subset
integrated_v2_LRP <- subset(integrated_v2_filtered, subset=major_label_transferAnchors_BaseAtlas %in% "LRP")
DefaultAssay(integrated_v2_LRP) <- "RNA"
dir.create('Atlas_v2_LRP_subset')
setwd('Atlas_v2_LRP_subset/')

# Remove mitochondrial and ribosomal genes from the LRP subset
mt_genes <- grep(pattern = "^mt-", x = rownames(integrated_v2_LRP@assays$RNA@counts), value = TRUE)
ribo_genes <- grep(pattern = "^Rp[sl]", x = rownames(integrated_v2_LRP@assays$RNA@counts), value = TRUE)
keep_features <- !grepl(paste0("^", mt_genes, collapse = "|"), rownames(integrated_v2_LRP@assays$RNA@counts))
integrated_v2_LRP <- subset(integrated_v2_LRP, features = rownames(integrated_v2_LRP@assays$RNA@counts)[keep_features])
keep_features <- !grepl(paste0("^", ribo_genes, collapse = "|"), rownames(integrated_v2_LRP@assays$RNA@counts))
integrated_v2_LRP <- subset(integrated_v2_LRP, features = rownames(integrated_v2_LRP@assays$RNA@counts)[keep_features])

# Split the LRP subset by sample
integrated_all_list <- SplitObject(integrated_v2_LRP, split.by = "sample")

# Normalize, find variable features, and scale the data for each sample in the LRP subset
for(sample in names(integrated_all_list)){
  obj <- integrated_all_list[[sample]]
  obj <- NormalizeData(obj)
  obj <- FindVariableFeatures(obj,nfeatures = 2000)
  integrated_all_list[[sample]]<-obj
}

i.anchors <- FindIntegrationAnchors(object.list = integrated_all_list, dims = 1:30,reduction='cca',scale=T,k.anchor=5,k.filter=100,k.score=15,anchor.features=3000)
integrated_v2_LRP <- IntegrateData(anchorset = i.anchors, dims = 1:30, normalization.method ='LogNormalize', k.weight=100)
DefaultAssay(integrated_v2_LRP) <- 'integrated'
integrated_v2_LRP <- ScaleData(integrated_v2_LRP,vars.to.regress = c("nFeature_RNA",'percent.mt','ccDiff'))
integrated_v2_LRP <- RunPCA(integrated_v2_LRP,npcs = 100)
data.use.integrated <- PrepDR(object = integrated_v2_LRP, genes.use = VariableFeatures(object = integrated_v2_LRP), use.imputed = F, assay.type = "integrated")
path_data <- getwd()

nPCs.data.use5 <- PCA_estimate_nPC(data.use.integrated, whereto=paste0(path_data,"/optimal_nPCs_5_integrated_filtered_LRP.RDS"), by.nPC=5,from.nPC = 10,to.nPC = 50) 
nPCs.data.use <- PCA_estimate_nPC(data.use.integrated, whereto=paste0(path_data,"/optimal_nPCs_integrated_filtered_LRP.RDS"), by.nPC=1,from.nPC = nPCs.data.use5-5,to.nPC = nPCs.data.use5+5) 

integrated_v2_LRP <- RunUMAP(integrated_v2_LRP, dims = 1:nPCs.data.use)
integrated_v2_LRP <- FindNeighbors(integrated_v2_LRP, dims=1:nPCs.data.use)
integrated_v2_LRP <- FindClusters(integrated_v2_LRP, resolution = c(0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.5))
#saveRDS(integrated_v2_LRP,file="Integrated_atlas_v2_FILTERED_LRP_cluster.RDS")

# Generate a UMAP plot of the data grouped by sample
DimPlot(integrated_v2_LRP,group.by = "sample",  reduction = 'umap')

# Generate UMAP plots for different cluster resolutions
pdf(paste0('umap_plot_atlas_v2_LRP_clusters.pdf'))
for(j in c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.5)){
  Idents(integrated_v2_LRP) <- paste0("integrated_snn_res.",j)
  print(DimPlot(integrated_v2_LRP, reduction = "umap",raster=FALSE, label=TRUE)  + ggtitle(paste0("Resolution: ",j)))
}  

# Generate a UMAP plot for nFeature_RNA expression
dittoDimPlot(integrated_v2_LRP, "nFeature_RNA", reduction.use = "umap", min.color = "lightgrey", max.color = "blue")

# Set the factor levels for minor labels
integrated_v2_LRP$minor_label_transferAnchors_BaseAtlas <- factor(integrated_v2_LRP$minor_label_transferAnchors_BaseAtlas, levels = c("LRP_1","LRP_2","LRP3","LRP4","Martinotti_1","Martinotti_2","Martinotti_3","Non-Martinotti_1","Non-Martinotti_2","Non-Martinotti_3","Non-Martinotti_4","Non-Martinotti_5"))

# Set cluster colors for visualization
ClusterCol <- c('#1E2D59', '#0C85B9', '#7CD5E2', '#4BC3AB','#61257A', '#C69BDF', '#D43DA7', '#C73A3A', '#F78173', '#FFB870', '#B59734', '#ECDD77')
names(ClusterCol) <- c("LRP_1","LRP_2","LRP3","LRP4","Martinotti_1","Martinotti_2","Martinotti_3","Non-Martinotti_1","Non-Martinotti_2","Non-Martinotti_3","Non-Martinotti_4","Non-Martinotti_5")

# Generate UMAP plots for minor label annotation
dittoDimPlot(integrated_v2_LRP,"minor_label_transferAnchors_BaseAtlas", split.by = "sample",color.panel=ClusterCol)
dittoDimPlot(integrated_v2_LRP,"minor_label_transferAnchors_BaseAtlas", color.panel=ClusterCol, size=2.5)

### Quality checks: Generate barplots and violin plots for various features across clusters
dittoBarPlot(integrated_v2_LRP, var = "sample", group.by = "integrated_snn_res.0.7")+ labs(title = NULL)

# Generate violin plots for Hsp90b1, Hspa5, percent ribosomal, mitochondrial, and feature counts
VlnPlot(integrated_v2_LRP, "Hsp90b1", pt.size = 0, group.by = "integrated_snn_res.0.7", cols = dittoColors())
VlnPlot(integrated_v2_LRP, "Hspa5", pt.size = 0, group.by = "integrated_snn_res.0.7", cols = dittoColors())
VlnPlot(integrated_v2_LRP, "percent.ribo", pt.size = 0, group.by = "integrated_snn_res.0.7")
VlnPlot(integrated_v2_LRP, "percent.mt", pt.size = 0, group.by = "integrated_snn_res.0.7")
VlnPlot(integrated_v2_LRP, "nFeature_RNA", pt.size = 0, group.by = "integrated_snn_res.0.7")


# Subset data for Martinotti cells and repeat integration, PCA, UMAP, clustering, and visualization (same steps as for LRP cells)
setwd('/lustre/projects/lim_core/results/tgi/Atlas_v2/integration_filtered_atlas/')
integrated_v2_Martinotti <- subset(integrated_v2, subset=major_label_transferAnchors_BaseAtlas %in% "Martinotti")
DefaultAssay(integrated_v2_Martinotti)<-"RNA"
dir.create('Atlas_v2_Martinotti_subset')
setwd('Atlas_v2_Martinotti_subset/')

# Remove mitochondrial and ribosomal genes for Martinotti subset
mt_genes <- grep(pattern = "^mt-", x = rownames(integrated_v2_Martinotti@assays$RNA@counts), value = TRUE)
ribo_genes <- grep(pattern = "^Rp[sl]", x = rownames(integrated_v2_Martinotti@assays$RNA@counts), value = TRUE)
keep_features <- !grepl(paste0("^", mt_genes, collapse = "|"), rownames(integrated_v2_Martinotti@assays$RNA@counts))
integrated_v2_Martinotti <- subset(integrated_v2_Martinotti, features = rownames(integrated_v2_Martinotti@assays$RNA@counts)[keep_features])
keep_features <- !grepl(paste0("^", ribo_genes, collapse = "|"), rownames(integrated_v2_Martinotti@assays$RNA@counts))
integrated_v2_Martinotti <- subset(integrated_v2_Martinotti, features = rownames(integrated_v2_Martinotti@assays$RNA@counts)[keep_features])

integrated_all_list <- SplitObject(integrated_v2_Martinotti, split.by = "sample")

for(sample in names(integrated_all_list)){
  obj <- integrated_all_list[[sample]]
  obj <- NormalizeData(obj)
  obj <- FindVariableFeatures(obj,nfeatures = 2000)
  integrated_all_list[[sample]]<-obj
}

i.anchors <- FindIntegrationAnchors(object.list = integrated_all_list, dims = 1:30,reduction='cca', scale=T,k.anchor=5, k.filter=100, k.score=15, anchor.features=3000)
integrated_v2_Martinotti <- IntegrateData(anchorset = i.anchors, dims = 1:30, normalization.method ='LogNormalize', k.weight=100)
DefaultAssay(integrated_v2_Martinotti) <- 'integrated'

integrated_v2_Martinotti <- ScaleData(integrated_v2_Martinotti, vars.to.regress = c("nFeature_RNA",'percent.mt','ccDiff'))
integrated_v2_Martinotti <- RunPCA(integrated_v2_Martinotti,npcs = 100)

data.use.integrated <- PrepDR(object = integrated_v2_Martinotti, genes.use = VariableFeatures(object = integrated_v2_Martinotti), use.imputed = F, assay.type = "integrated")
path_data <- getwd()

if (!file.exists(paste0(path_data, "/optimal_nPCs_5_integrated_filtered_Martinotti.RDS"))) { # Check if the optimal nPCs file for nPCs.data.use5 exists

  # If the file doesn't exist, calculate and save the nPCs for the integrated dataset
  nPCs.data.use5 <- PCA_estimate_nPC(data.use.integrated, 
                                     whereto = paste0(path_data, "/optimal_nPCs_5_integrated_filtered_Martinotti.RDS"), 
                                     by.nPC = 5, from.nPC = 30, to.nPC = 70)
} else {
  # If the file exists, load the precomputed nPCs
  nPCs.data.use5 <- readRDS(paste0(path_data, "/optimal_nPCs_5_integrated_filtered_Martinotti.RDS"))
}


if (!file.exists(paste0(path_data, "/optimal_nPCs_integrated_filtered_Martinotti.RDS"))) { # Check if the second nPCs file for nPCs.data.use exists
  # If the file doesn't exist, calculate and save the nPCs for the integrated dataset
  nPCs.data.use <- PCA_estimate_nPC(data.use.integrated, 
                                    whereto = paste0(path_data, "/optimal_nPCs_integrated_filtered_Martinotti.RDS"), 
                                    by.nPC = 1, from.nPC = nPCs.data.use5 - 5, to.nPC = nPCs.data.use5 + 5) #26
} else {
  # If the file exists, load the precomputed nPCs
  nPCs.data.use <- readRDS(paste0(path_data, "/optimal_nPCs_integrated_filtered_Martinotti.RDS"))
}

integrated_v2_Martinotti <- RunUMAP(integrated_v2_Martinotti, dims=1:nPCs.data.use)
integrated_v2_Martinotti <- FindNeighbors(integrated_v2_Martinotti, dims=1:nPCs.data.use)
integrated_v2_Martinotti <- FindClusters(integrated_v2_Martinotti, resolution = c(0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.5))
#saveRDS(integrated_v2_Martinotti,file="Integrated_atlas_v2_FILTERED_Martinotti_cluster.RDS")

#Plot UMAP to see martinotti distribution over samples
DimPlot(integrated_v2_Martinotti,group.by="sample",  reduction = 'umap')

# Generate UMAP plots for different cluster resolutions
pdf(paste0('umap_plot_atlas_v2_Martinotti_clusters.pdf'))
for(j in c(0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.5)){
  Idents(integrated_v2_Martinotti) <- paste0("integrated_snn_res.",j)
  print(DimPlot(integrated_v2_Martinotti, reduction = "umap",raster=FALSE, label=TRUE)  + ggtitle(paste0("Resolution: ",j)))
}  

# Generate a UMAP plot for nFeature_RNA expression
dittoDimPlot(integrated_v2_Martinotti, "nFeature_RNA", reduction.use = "umap", min.color = "lightgrey", max.color = "blue")

# Set the factor levels for minor labels
integrated_v2_Martinotti$minor_label_transferAnchors_BaseAtlas <- factor(integrated_v2_Martinotti$minor_label_transferAnchors_BaseAtlas, levels = c("LRP_1","LRP_2","LRP3","LRP4","Martinotti_1","Martinotti_2","Martinotti_3","Non-Martinotti_1","Non-Martinotti_2","Non-Martinotti_3","Non-Martinotti_4","Non-Martinotti_5"))

# Set cluster colors for visualization
ClusterCol <- c('#1E2D59', '#0C85B9', '#7CD5E2', '#4BC3AB','#61257A', '#C69BDF', '#D43DA7', '#C73A3A', '#F78173', '#FFB870', '#B59734', '#ECDD77')
names(ClusterCol) <- c("LRP_1","LRP_2","LRP3","LRP4","Martinotti_1","Martinotti_2","Martinotti_3","Non-Martinotti_1","Non-Martinotti_2","Non-Martinotti_3","Non-Martinotti_4","Non-Martinotti_5")

# Generate UMAP plots for minor label annotation
dittoDimPlot(integrated_v2_Martinotti,"minor_label_transferAnchors_BaseAtlas", split.by = "sample",color.panel=ClusterCol)
dittoDimPlot(integrated_v2_Martinotti,"minor_label_transferAnchors_BaseAtlas", color.panel=ClusterCol, size=1.5)


### Quality checks
dittoBarPlot(integrated_v2_Martinotti, var = "sample", group.by = "integrated_snn_res.0.7")+ labs(title = NULL)
VlnPlot(integrated_v2_Martinotti, "Hsp90b1", pt.size = 0, group.by = "integrated_snn_res.0.7", cols = dittoColors())
VlnPlot(integrated_v2_Martinotti, "Hspa5", pt.size = 0, group.by = "integrated_snn_res.0.7", cols = dittoColors())
VlnPlot(integrated_v2_Martinotti, "percent.ribo", pt.size = 0, group.by = "integrated_snn_res.0.7")
VlnPlot(integrated_v2_Martinotti, "percent.mt", pt.size = 0, group.by = "integrated_snn_res.0.7")
VlnPlot(integrated_v2_Martinotti, "nFeature_RNA", pt.size = 0, group.by = "integrated_snn_res.0.7")


# Subset data for Non-Martinotti cells and repeat integration, PCA, UMAP, clustering, and visualization (same steps as for Martinotti cells)
setwd('/lustre/projects/lim_core/results/tgi/Atlas_v2/integration_filtered_atlas/')
integrated_v2_NonMartinotti <- subset(integrated_v2, subset=major_label_transferAnchors_BaseAtlas %in% "Non-Martinotti")
DefaultAssay(integrated_v2_NonMartinotti) <- "RNA"
dir.create('Atlas_v2_NonMartinotti_subset')
setwd('Atlas_v2_NonMartinotti_subset/')

#Remove mitochondrial and ribosomal genes
mt_genes <- grep(pattern = "^mt-", x = rownames(integrated_v2_NonMartinotti@assays$RNA@counts), value = TRUE)
ribo_genes <- grep(pattern = "^Rp[sl]", x = rownames(integrated_v2_NonMartinotti@assays$RNA@counts), value = TRUE)
keep_features <- !grepl(paste0("^", mt_genes, collapse = "|"), rownames(integrated_v2_NonMartinotti@assays$RNA@counts))
integrated_v2_NonMartinotti <- subset(integrated_v2_NonMartinotti, features = rownames(integrated_v2_NonMartinotti@assays$RNA@counts)[keep_features])
keep_features <- !grepl(paste0("^", ribo_genes, collapse = "|"), rownames(integrated_v2_NonMartinotti@assays$RNA@counts))
integrated_v2_NonMartinotti <- subset(integrated_v2_NonMartinotti, features = rownames(integrated_v2_NonMartinotti@assays$RNA@counts)[keep_features])


integrated_all_list <- SplitObject(integrated_v2_NonMartinotti, split.by = "sample")

for(sample in names(integrated_all_list)){
  obj <- integrated_all_list[[sample]]
  obj <- NormalizeData(obj)
  obj <- FindVariableFeatures(obj,nfeatures = 2000)
  integrated_all_list[[sample]] <- obj
}

i.anchors <- FindIntegrationAnchors(object.list = integrated_all_list, dims = 1:30, reduction='cca',scale=T, k.anchor=5, k.filter=100, k.score=15, anchor.features=3000)
integrated_v2_NonMartinotti <- IntegrateData(anchorset = i.anchors, dims = 1:30, normalization.method ='LogNormalize', k.weight=100)
DefaultAssay(integrated_v2_NonMartinotti) <- 'integrated'
integrated_v2_NonMartinotti <- ScaleData(integrated_v2_NonMartinotti, vars.to.regress = c("nFeature_RNA",'percent.mt','ccDiff'))
integrated_v2_NonMartinotti <- RunPCA(integrated_v2_NonMartinotti,npcs = 100)
data.use.integrated<- PrepDR(object = integrated_v2_NonMartinotti, genes.use = VariableFeatures(object = integrated_v2_NonMartinotti), use.imputed = F, assay.type = "integrated")
path_data <- getwd()

# Check if the optimal nPCs file for nPCs.data.use5 exists
if (!file.exists(paste0(path_data, "/optimal_nPCs_5_integrated_filtered_NonMartinotti.RDS"))) {
  # If the file doesn't exist, calculate and save the nPCs for the integrated dataset
  nPCs.data.use5 <- PCA_estimate_nPC(data.use.integrated, 
                                     whereto = paste0(path_data, "/optimal_nPCs_5_integrated_filtered_NonMartinotti.RDS"), 
                                     by.nPC = 5, from.nPC = 30, to.nPC = 70)
} else {
  # If the file exists, load the precomputed nPCs
  nPCs.data.use5 <- readRDS(paste0(path_data, "/optimal_nPCs_5_integrated_filtered_NonMartinotti.RDS"))
}

# Check if the second nPCs file for nPCs.data.use exists
if (!file.exists(paste0(path_data, "/optimal_nPCs_integrated_filtered_NonMartinotti.RDS"))) {
  # If the file doesn't exist, calculate and save the nPCs for the integrated dataset
  nPCs.data.use <- PCA_estimate_nPC(data.use.integrated, 
                                    whereto = paste0(path_data, "/optimal_nPCs_integrated_filtered_NonMartinotti.RDS"), 
                                    by.nPC = 1, from.nPC = nPCs.data.use5 - 5, to.nPC = nPCs.data.use5 + 5)
} else {
  # If the file exists, load the precomputed nPCs
  nPCs.data.use <- readRDS(paste0(path_data, "/optimal_nPCs_integrated_filtered_NonMartinotti.RDS"))
}

integrated_v2_NonMartinotti <- RunUMAP(integrated_v2_NonMartinotti,dims=1:nPCs.data.use)
integrated_v2_NonMartinotti <- FindNeighbors(integrated_v2_NonMartinotti,dims=1:nPCs.data.use)
integrated_v2_NonMartinotti <- FindClusters(integrated_v2_NonMartinotti,resolution = c(0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.5))
#saveRDS(integrated_v2_NonMartinotti,file="Integrated_atlas_v2_FILTERED_NonMartinotti_cluster.RDS")


#Plot UMAP to see martinotti distribution over samples
DimPlot(integrated_v2_NonMartinotti,group.by="sample",  reduction = 'umap')
dev.off()

# Generate UMAP plots for different cluster resolutions
pdf(paste0('umap_plot_atlas_v2_NonMartinotti_clusters.pdf'))
for(j in c(0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.5)){
  Idents(integrated_v2_NonMartinotti) <- paste0("integrated_snn_res.",j)
  print(DimPlot(integrated_v2_NonMartinotti, reduction = "umap", raster=FALSE, label=TRUE)  + ggtitle(paste0("Resolution: ",j)))
}  
dev.off()

# Generate a UMAP plot for nFeature_RNA expression
dittoDimPlot(integrated_v2_NonMartinotti, "nFeature_RNA", reduction.use = "umap", min.color = "lightgrey", max.color = "blue")
dev.off()

# Set the factor levels for minor labels
integrated_v2_NonMartinotti$minor_label_transferAnchors_BaseAtlas <- factor(integrated_v2_NonMartinotti$minor_label_transferAnchors_BaseAtlas, levels = c("LRP_1","LRP_2","LRP3","LRP4","Martinotti_1","Martinotti_2","Martinotti_3","Non-Martinotti_1","Non-Martinotti_2","Non-Martinotti_3","Non-Martinotti_4","Non-Martinotti_5"))

# Set cluster colors for visualization
ClusterCol <- c('#1E2D59', '#0C85B9', '#7CD5E2', '#4BC3AB','#61257A', '#C69BDF', '#D43DA7', '#C73A3A', '#F78173', '#F69635', '#B59734', '#ECDD77')
names(ClusterCol) <- c("LRP_1","LRP_2","LRP3","LRP4","Martinotti_1","Martinotti_2","Martinotti_3","Non-Martinotti_1","Non-Martinotti_2","Non-Martinotti_3","Non-Martinotti_4","Non-Martinotti_5")

# Generate UMAP plots for minor label annotation
dittoDimPlot(integrated_v2_NonMartinotti, "minor_label_transferAnchors_BaseAtlas", split.by = "sample",color.panel=ClusterCol)
dittoDimPlot(integrated_v2_NonMartinotti,"minor_label_transferAnchors_BaseAtlas", color.panel=ClusterCol, size=1.5)

rm()
gc()
```


