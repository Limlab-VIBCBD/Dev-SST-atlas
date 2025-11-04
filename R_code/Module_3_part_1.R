## Installation needed
if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")

pacman::p_load(
  Seurat, dittoSeq, ggplot2, gtools, dplyr, R.utils, colorspace, tidyr,
  purrr, Matrix, patchwork, BiocManager
)
devtools::install_github("crazyhottommy/scclusteval")
BiocManager::install("scDblFinder")
BiocManager::install("MAST")
BiocManager::install("dittoSeq")

########################
# Load required packages
########################
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

###########
# Load data
###########
# Load previously Sst+ filtered cells from saved Seurat objects passing integrability tests
base_atlas <- readRDS("integrated_all_to_send.rds")
base_atlas$batch <- paste0("BaseAtlas_",base_atlas$orig.ident)
E16_Lim4 <- readRDS("lim4_E16_Sst_filtered.rds")
E16_Lim4$batch <- "E16_Lim4"
P1_Lim3 <- readRDS('lim3_P1_Sst_filtered.rds')
P1_Lim3$batch <- "P1_Lim3"
P1_Lim5 <- readRDS("lim5_P1_Sst_filtered.rds")
P1_Lim5$batch <- "P1_Lim5"
test_dataset_filtered_sst<-readRDS("P5_WT1_WT23_Sst_filtered.rds")
P5_WT1_Lim1 <- subset(test_dataset_filtered_sst, subset= orig.ident %in% "lim_P5_sorted_sept23_WT1")
P5_WT1_Lim1$batch <- "P5_WT1_Lim1"
P5_WT23_Lim2 <- subset(test_dataset_filtered_sst, subset= orig.ident %in% "lim_P5_sorted_sept23_WT23")
P5_WT23_Lim2$batch <- "P5_WT23_Lim2"
E18_Lippi<-readRDS("Lippi_lab_e18.5_Sst_filtered.rds")
E18_Lippi$batch <- "E18_Lippi"
P1_EMI014_Lim<-readRDS("lim_P1_EMI014_Sst_filtered.rds")
P1_EMI014_Lim$batch <- "P1_EMI014_Lim"
Wu_data <- readRDS("P2_P7_GSE272706_Sst_filtered.RDS")
list_Wu <-SplitObject(Wu_data, split.by = "orig.ident")
Wu_p2 <- list_Wu[["GSM8409508_P2MGE_Fezf2HET_BaxcHET"]]
Wu_p2$batch <- "P2_Wu"
E16_EMI018 <-readRDS("lim_E16.5_EMI018_Sst_filtered.rds")
E16_EMI018$batch <- 'E16_EMI018'
# Set the default assay to RNA for analysis
DefaultAssay(base_atlas) <- "RNA"
DefaultAssay(E16_Lim4) <- "RNA"
DefaultAssay(P1_Lim3) <- "RNA"
DefaultAssay(P1_Lim5) <- "RNA"
DefaultAssay(P5_WT1_Lim1) <- "RNA"
DefaultAssay(P5_WT23_Lim2) <- "RNA"
DefaultAssay(E18_Lippi) <- "RNA"
DefaultAssay(P1_EMI014_Lim) <- "RNA"
DefaultAssay(Wu_p2)<-"RNA"
DefaultAssay(E16_EMI018)<-"RNA"
# Split the base atlas into separate objects based on the sampleID 
integrated_all_list <- SplitObject(base_atlas, split.by = "batch")
# Add the cleaned dataset to the list of samples
integrated_all_list$E16_Lim4 <- E16_Lim4
integrated_all_list$P1_Lim3 <- P1_Lim3
integrated_all_list$P1_Lim5 <- P1_Lim5
integrated_all_list$P5_WT1_Lim1 <- P5_WT1_Lim1
integrated_all_list$P5_WT23_Lim2 <- P5_WT23_Lim2
integrated_all_list$E18_Lippi <- E18_Lippi
integrated_all_list$P1_EMI014_Lim <- P1_EMI014_Lim
integrated_all_list$Wu_p2 <- Wu_p2
integrated_all_list$E16_EMI018 <- E16_EMI018

#########################
# Doublets identification
#########################
# First, we annotate cells according to Atlas-v0 clusters using label transfer.  
# The predicted clusters will be used to identify potential doublets.
# Normalize the data, identify variable features, and scale the Atlas-v0 dataset.
base_atlas <- NormalizeData(base_atlas)
base_atlas <- FindVariableFeatures(base_atlas)
base_atlas <- ScaleData(base_atlas, vars.to.regress = c("nFeature_RNA",'percent.mt','ccDiff'))
base_atlas <- RunPCA(base_atlas,npcs = 50)
# Initialize empty vector to store minor labels
minor_label <- c()
# Perform label transfer
for (sample in c("E16_Lim4","P1_Lim3","P1_Lim5","P5_WT1_Lim1","P5_WT23_Lim2","E18_Lippi","P1_EMI014_Lim","Wu_p2","E16_EMI018")) {
  integrated_all_list[[sample]] <- NormalizeData(integrated_all_list[[sample]])
  integrated_all_list[[sample]] <- FindVariableFeatures(integrated_all_list[[sample]])
  integrated_all_list[[sample]] <- ScaleData(integrated_all_list[[sample]], vars.to.regress = c("nFeature_RNA",'percent.mt','ccDiff'))
  integrated_all_list[[sample]] <- RunPCA(integrated_all_list[[sample]],npcs = 50)
  # Find transfer anchors between the Atlas-v0 (reference) and the samples
  transfer.anchors <- FindTransferAnchors(reference = base_atlas, query = integrated_all_list[[sample]], dims = 1:40, reference.reduction = "pca", features=intersect(rownames(base_atlas), rownames(integrated_all_list[[sample]])))
  # Transfer labels based on the reference
  predictions <- TransferData(anchorset = transfer.anchors, refdata = base_atlas$cluster_label_trained_with_all, dims = 1:40)
  # Add the predicted minor labels to the sample object
  integrated_all_list[[sample]]$minor_label_transferAnchors_BaseAtlas <- predictions[colnames(integrated_all_list[[sample]]),'predicted.id']
  # Append the predicted labels to the overall list of minor labels
  minor_pred <- predictions$predicted.id
  names(minor_pred) <- rownames(predictions)
  minor_label <- c(minor_label,minor_pred)
}
for (sample in c("BaseAtlas_E16","BaseAtlas_lim_P5_fixed_sorted","BaseAtlas_P1","BaseAtlas_P5") ) {
  integrated_all_list[[sample]]$minor_label_transferAnchors_BaseAtlas <- integrated_all_list[[sample]]$cluster_label_trained_with_all
}
# Doublets Prediction Using Minor Labels
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
  integrated_all_list[[i]]$doublets <- doublets_anno[colnames(integrated_all_list[[i]])]
  # Subset the Seurat object to retain only the singlets
  integrated_all_list[[i]] <- subset(integrated_all_list[[i]], subset=doublets=="singlet")
}

################################
# CCA Integration and clustering
################################
# Normalize data and find variable features
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
integrated_v2 <- RunPCA(integrated_v2, npcs = 60)
ElbowPlot(integrated_v2, ndims = 60)
# Determine the optimal number of principal components (PCs) for downstream analysis
data.use.integrated <- PrepDR(object = integrated_v2, genes.use = VariableFeatures(object = integrated_v2), use.imputed = F, assay.type = "integrated")
path_data <- getwd()
nPCs.data.use5 <- PCA_estimate_nPC(data.use.integrated, 
                                     whereto = paste0(path_data, "/optimal_nPCs_5_integrated.RDS"), 
                                     k = 5, by.nPC = 5, from.nPC = 30, to.nPC = 50) #check ElbowPlot and adjust range acoordingly
nPCs.data.use <- PCA_estimate_nPC(data.use.integrated, 
                                    whereto = paste0(path_data, "/optimal_nPCs_integrated.RDS"), 
                                    k = 5, by.nPC = 1, from.nPC = nPCs.data.use5 - 5, to.nPC = nPCs.data.use5 + 5)
integrated_v2 <- RunUMAP(integrated_v2, dims=1:nPCs.data.use)
integrated_v2 <- FindNeighbors(integrated_v2, dims=1:nPCs.data.use)
integrated_v2 <- FindClusters(integrated_v2, resolution = c(0.4))

######################################
# PERFORM QC AND REMOVE STRESSED CELLS 
######################################
DefaultAssay(integrated_v2) <- 'RNA'
integrated_v2 <- ScaleData(integrated_v2)
# Extract UMAP embeddings (2D coordinates for visualization)
umap <- Embeddings(integrated_v2,reduction = "umap")
# Loop through a list of stress-related genes to visualize their expression on UMAP and export in pdf
for (gene in c("Hsp90b1","Hspa5","Mapk8")){
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
pdf('umap_plot_atlas_v2_clusters.pdf')
  Idents(integrated_v2) <- "integrated_snn_res.0.4"
  DimPlot(integrated_v2, reduction = "umap",raster=FALSE, label=TRUE) + ggtitle("Resolution: 0.4")
dev.off()
# Set the clustering identity to the specific resolution (0.4) and generate bar plots of samples per cluster
Idents(integrated_v2)<-'integrated_snn_res.0.4'
dittoBarPlot(integrated_v2, var = "batch", group.by = "integrated_snn_res.0.4")  + labs(title = NULL)
# Generate violin plots for expression of various genes by cluster (based on the integrated SNN resolution 0.4)
VlnPlot(integrated_v2, "Hsp90b1", pt.size = 0, group.by = "integrated_snn_res.0.4")
VlnPlot(integrated_v2, "Hspa5", pt.size = 0, group.by = "integrated_snn_res.0.4")
VlnPlot(integrated_v2, "Mapk8", pt.size = 0, group.by = "integrated_snn_res.0.4")
# Add a new feature for ribosomal RNA percentage (ribo genes)
integrated_v2[["percent.ribo"]] <- PercentageFeatureSet(integrated_v2, pattern = "^Rp[Sl]") 
# Violin plot for ribosomal RNA percentage by cluster
VlnPlot(integrated_v2, "percent.ribo", pt.size = 0, group.by = "integrated_snn_res.0.4")
# Violin plot for mitochondrial RNA percentage by cluster
VlnPlot(integrated_v2, "percent.mt", pt.size = 0, group.by = "integrated_snn_res.0.4")
# Violin plots for total RNA count and feature count by cluster
VlnPlot(integrated_v2, "nCount_RNA", pt.size = 0, group.by = "integrated_snn_res.0.4") + ylim(0,30000)
VlnPlot(integrated_v2, "nFeature_RNA", pt.size = 0, group.by = "integrated_snn_res.0.4")
# Filter cells - remove clusters with high expression of stress markers or mitochondrial genes, as well as clusters with low nFeatures or unbalanced sample composition.
# Specify clusters you want to remove
integrated_v2_filtered<-subset(integrated_v2, subset=integrated_snn_res.0.4 %in% c(8,14,15), invert = TRUE)

#################################
# CCA Integration after filtering
#################################
DefaultAssay(integrated_v2_filtered) <- 'RNA'
integrated_all_list <- SplitObject(integrated_v2_filtered, split.by = "batch")
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
integrated_v2_filtered <- ScaleData(integrated_v2_filtered, verbose = FALSE, vars.to.regress = c("nFeature_RNA",'percent.mt','ccDiff'))
integrated_v2_filtered <-RunPCA(integrated_v2_filtered, npcs = 60)
ElbowPlot(integrated_v2, ndims = 60)
# Determine the optimal number of principal components (PCs) for downstream analysis
data.use.integrated<- PrepDR(object = integrated_v2_filtered, genes.use = VariableFeatures(object = integrated_v2_filtered), use.imputed = F, assay.type = "integrated")
path_data <- getwd()
nPCs.data.use5 <- PCA_estimate_nPC(data.use.integrated, 
                                     whereto = paste0(path_data, "/optimal_nPCs_5_integrated_filtered.RDS"), 
                                     k = 5, by.nPC = 5, from.nPC = 40, to.nPC = 60) # Check Elbow plot and set the range for optimal  PC
nPCs.data.use <- PCA_estimate_nPC(data.use.integrated, 
                                    whereto = paste0(path_data, "/optimal_nPCs_integrated_filtered.RDS"), 
                                    k = 2, by.nPC = 1, from.nPC = nPCs.data.use5 - 5, to.nPC = nPCs.data.use5 + 5)
integrated_v2_filtered <- RunUMAP(integrated_v2_filtered, dims=1:nPCs.data.use)
integrated_v2_filtered <- FindNeighbors(integrated_v2_filtered, dims=1:nPCs.data.use)
integrated_v2_filtered <- FindClusters(integrated_v2_filtered, resolution = c(0.4,0.5,0.6,0.7))
# Here, another round of quality checking is performed

##################
# Cells annotation
##################
# We annotate cells to major classes using Atlas-v0 as reference
DefaultAssay(integrated_v2_filtered) <- 'RNA'
integrated_all_list <- SplitObject(integrated_v2_filtered, split.by = "batch")
base_atlas <- subset(integrated_v2_filtered, subset= batch %in% c("BaseAtlas_E16","BaseAtlas_lim_P5_fixed_sorted","BaseAtlas_P1","BaseAtlas_P5"))
DefaultAssay(base_atlas) <- 'RNA'
base_atlas <- NormalizeData(base_atlas)
base_atlas <- FindVariableFeatures(base_atlas)
base_atlas <- ScaleData(base_atlas, vars.to.regress = c("nFeature_RNA",'percent.mt','ccDiff'))
base_atlas<-RunPCA(base_atlas,npcs = 60)
# Initialize empty vectors to store predicted labels and scores
minor_label <- c()
major_label <- c()
pred_score_minor  <- c()
pred_score_major  <- c()
# Loop through each sample, normalize, identify variable features, scale data, and run PCA
for (sample in c("E16_Lim4","P1_Lim3","P1_Lim5","P5_WT1_Lim1","P5_WT23_Lim2","E18_Lippi","P1_EMI014_Lim","Wu_p2","E16_EMI018")) {
  # Normalize, find variable features, and scale the data for each sample
  integrated_all_list[[sample]] <- NormalizeData(integrated_all_list[[sample]])
  integrated_all_list[[sample]] <- FindVariableFeatures(integrated_all_list[[sample]])
  integrated_all_list[[sample]] <- ScaleData(integrated_all_list[[sample]], vars.to.regress = c("nFeature_RNA",'percent.mt','ccDiff'))
  integrated_all_list[[sample]] <- RunPCA(integrated_all_list[[sample]],npcs = 50)
  # Perform anchor finding for label transfer based on PCA reduction and intersected features
  transfer.anchors <- FindTransferAnchors(reference = base_atlas, query = integrated_all_list[[sample]], dims = 1:40,reference.reduction = "pca", features=intersect(rownames(base_atlas), rownames(integrated_all_list[[sample]])))
  # Transfer the predicted minor labels (clusters)
  predictions <- TransferData(anchorset = transfer.anchors, refdata = base_atlas$cluster_label_trained_with_all, dims = 1:40)
  minor_pred <- predictions$predicted.id
  names(minor_pred) <- rownames(predictions)
  minor_label <- c(minor_label,minor_pred)
  # Capture the prediction score for each minor label prediction
  pred_score_tmp <- predictions$prediction.score.max
  names(pred_score_tmp) <- rownames(predictions)
  pred_score_minor <- c(pred_score_minor,pred_score_tmp)
  # Transfer the predicted major labels (major subsets)
  predictions <- TransferData(anchorset = transfer.anchors, refdata = base_atlas$major_cluster_label_trained_with_all, dims = 1:40)
  major_pred <- predictions$predicted.id
  names(major_pred) <- rownames(predictions)
  major_label <- c(major_label,major_pred)
    # Capture the prediction score for each major label prediction
  pred_score_tmp <- predictions$prediction.score.max
  names(pred_score_tmp) <- rownames(predictions)
  pred_score_major <- c(pred_score_major,pred_score_tmp)
}
saveRDS(minor_label,"Minor_label_labelTranfer.RDS")
saveRDS(major_label,"Major_label_labelTranfer.RDS")
# Assign minor and major labels to the meta-data of the integrated filtered dataset
integrated_v2_filtered$minor_label_transferAnchors_BaseAtlas <- integrated_v2_filtered$cluster_label_trained_with_all
integrated_v2_filtered@meta.data[names(minor_label),'minor_label_transferAnchors_BaseAtlas'] <- minor_label
integrated_v2_filtered$major_label_transferAnchors_BaseAtlas <- integrated_v2_filtered$major_cluster_label_trained_with_all
integrated_v2_filtered@meta.data[names(major_label),'major_label_transferAnchors_BaseAtlas'] <- major_label
saveRDS(integrated_v2_filtered, "Integrated_atlas_V2_filtered_annotated.RDS")
# Generate UMAP plots for major labels
ClusterCol <- c('#35C1D5', '#4A2884',  '#E05F36')
DimPlot(integrated_v2_filtered, group.by="major_label_transferAnchors_BaseAtlas",  reduction = 'umap', cols=ClusterCol)

#######################################
# SAVE objects for clusters computation
#######################################
integrated_v2_filtered_LRP <- subset(integrated_v2_filtered, subset=major_label_transferAnchors_BaseAtlas == "LRP")
saveRDS(integrated_v2_filtered_LRP,'Integrated_atlas_v2_LRP.RDS')
integrated_v2_filtered_Martinotti <- subset(integrated_v2_filtered, subset=major_label_transferAnchors_BaseAtlas == "Martinotti")
saveRDS(integrated_v2_filtered_Martinotti,'Integrated_atlas_v2_Martinotti.RDS')
integrated_v2_filtered_NonMartinotti <- subset(integrated_v2_filtered, subset=major_label_transferAnchors_BaseAtlas == "Non-Martinotti")
saveRDS(integrated_v2_filtered_NonMartinotti,'Integrated_atlas_v2_NonMartinotti.RDS')

