## Installation needed:
if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(
  BiocManager, remotes, devtools, Seurat, dittoSeq, ggplot2, gtools, dplyr, R.utils, DoubletFinder,
  dismo, missMDA, hdf5r
)
#remotes::install_github("chris-mcginnis-ucsf/DoubletFinder") #install if needed
#BiocManager::install("dittoSeq")

##########################
# Load required packages
##########################

library(Seurat)
library(dittoSeq)
library(ggplot2)
library(gtools)
library(dplyr)
library(dittoSeq)
library(R.utils)
library(DoubletFinder)
#Source the external script 'Seurat_Utils' for further analysis (10-fold Singular Value Decomposition (SVD) cross validation to predict number of PCs)
source("Seurat_Utils.R") 

##########################
# Load data
##########################

#Loads scRNA-seq count matrices from CellRanger output. Ensure the folder containing your files is located in your working directory
#option 1 : Cellranger output files
#Dataset <- Read10X("filtered_feature_bc_matrix/") 
#Option 2: If your data is in an h5 file format, use the following line instead:
Dataset <- Read10X_h5("Module 1/Dataset_examples/h5 files/P1_DiBella.h5", use.names = TRUE, unique.features = TRUE)
#Create a Seurat object from the count matrices and specify the associated project: 
seurat_object <- CreateSeuratObject(counts = Dataset, project = "P1", min.cells = 10, min.features = 0)
#Add a metadata column to define the sample name:
seurat_object$sample <- "P1_Dibella" 
rm(Dataset)
gc()

###########################
# Perform QC and filtering
###########################

# Calculate the percentage of mitochondrial and ribosomal gene expression. These metrics indicate cell quality (e.g., high mitochondrial content suggests cell stress or death).
seurat_object[["percent.mt"]] <- PercentageFeatureSet(seurat_object, pattern = "^mt-")
seurat_object[["percent.ribo"]] <- PercentageFeatureSet(seurat_object, pattern = "^Rp[Sl]") # nolint: line_length_linter.
# Visualize the distributions of mitochondrial and ribosomal gene percentages, the number of expressed genes per cell, and the number of UMIs per cell
VlnPlot(seurat_object, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.ribo"), ncol = 2, pt.size = 0, group.by = "sample")
## Cell and gene filtering (thresholds are user dependent)
# Set the minimum number of UMIs per cell (calculated as the mean - 2x SD).
logUMI <- log1p(seurat_object[['nCount_RNA']])
lower.umi_WT1 <- exp(mean(logUMI[,1]) - (2*sd(logUMI[,1])))
# Remove cells based on UMI counts, percentage of mitochondrial genes and number of genes
Seurat_object_subset <- subset(seurat_object, subset = nFeature_RNA > 700 & percent.mt < 10 & nCount_RNA > lower.umi_WT1)
rm(list = c("logUMI", "seurat_object", "lower.umi_WT1"))
gc()


###########################
# Predict cell cycle
###########################

#Load the two external files containing gene names
s.genes <- readRDS("mouse_s.genes.rds")
g2m.genes <- readRDS("mouse_g2m.genes.rds")
#Assigns cells to cell cycle phases (S, G2M, or G1) based on marker genes
Seurat_object_subset <- CellCycleScoring(Seurat_object_subset, s.features = s.genes, g2m.features = g2m.genes)
ScoreDiffSeuratObject <- (Seurat_object_subset@meta.data[["S.Score"]] - Seurat_object_subset@meta.data[["G2M.Score"]])
Seurat_object_subset[['ccDiff']] <- ScoreDiffSeuratObject
rm(list = c("g2m.genes", "s.genes", "ScoreDiffSeuratObject"))


###########################
# Normalization and scaling
###########################

Seurat_object_subset <- NormalizeData(Seurat_object_subset)
Seurat_object_subset <- FindVariableFeatures(Seurat_object_subset, nfeatures = 2000)
Seurat_object_subset <- ScaleData(Seurat_object_subset,vars.to.regress = c("nFeature_RNA",'percent.mt','ccDiff'), model.use = 'linear', block.size=dim(Seurat_object_subset)[1])


#########
# PCA
#########

Seurat_object_subset <- RunPCA(Seurat_object_subset, npcs = 60)
# Determine the optimal number of principal components (PCs) for downstream analysis
ElbowPlot(Seurat_object_subset, ndims = 60)
# Adjusting nPC based on ndims on the ElbowPlot 
data.use.Seurat_object_subset <- PrepDR(object = Seurat_object_subset, genes.use = VariableFeatures(object = Seurat_object_subset), use.imputed = F, assay.type = "RNA")
path_data <- getwd()
# Adjust `from.nPC` and `to.nPC` after inspecting the ndims of the ElbowPlot
nPCs.Seurat_object_subset <- PCA_estimate_nPC(data.use.Seurat_object_subset, whereto = paste0(path_data,"/optimal_nPCs_WT1.RDS"), k = 5, by.nPC = 2, from.nPC = 30, to.nPC = 40) 
rm(path_data)
rm(data.use.Seurat_object_subset)
gc()

#################################
### Compute UMAP/tSNE and clusters
################################

# Perform tSNE and UMAP dimensionality reduction using the top PCs
Seurat_object_subset <- RunTSNE(Seurat_object_subset, dims = 1:nPCs.Seurat_object_subset)
Seurat_object_subset <- RunUMAP(Seurat_object_subset, dims = 1:nPCs.Seurat_object_subset)
# Find neighbors and clusters based on reduced dimensions
Seurat_object_subset <- FindNeighbors(Seurat_object_subset, dims = 1:nPCs.Seurat_object_subset)
Seurat_object_subset <- FindClusters(Seurat_object_subset, resolution = 0.8)

################################################
### Remove low quality clusters and contaminants
###############################################

# This visualizes the clusters identified by Seurat's clustering on UMAP
dittoDimPlot(Seurat_object_subset, "seurat_clusters", reduction.use = "umap", split.by = "sample", do.label = T)
#Create diagnostic plots to check quality control metrics on UMAP (Clipping outlier values to better visualization on the graph, the data itself it is not clipped)
clipped_data <- Seurat_object_subset
clipped_data$nCount_RNA <- pmin(pmax(Seurat_object_subset$nCount_RNA, 0), 10000) #remove higher values to rescale the color map (only for plotting on the umap)
clipped_data$nFeature_RNA <- pmin(pmax(Seurat_object_subset$nFeature_RNA, 2315), 7500) #remove higher values to rescale the color map (only for plotting on the umap)
#Plot these QC metrics on the scaled UMAP
dittoDimPlot(clipped_data, "nCount_RNA", reduction.use = "umap", min.color = "lightgrey", max.color = "blue", min = 0, max = 20000, split.by = "sample")
dittoDimPlot(clipped_data, "nFeature_RNA", reduction.use = "umap", min.color = "lightgrey", max.color = "blue", split.by = "sample")
dittoDimPlot(Seurat_object_subset, "percent.mt", reduction.use = "umap", min.color = "lightgrey", max.color = "blue", split.by = "sample")
VlnPlot(Seurat_object_subset, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.ribo"), ncol = 2, pt.size = 0, group.by = "seurat_clusters")
rm(clipped_data)
gc()

## User define genes to visual which clusters will be selected. We exluded (1) Stress clusters, (2) negative markers - pyramidal, striatal, (3) positve markers - in this case SST, Lhx6 etc
# Specificity plots (UMAP + Violin plot) - to extract cortical SST+ cells, we use Lhx6 and Sst 
# Stressed markers visualization (UMAP + Violin plot)
Stress <- c("Hspa5", "Mapk8", "Hsp90b1")
dittoDimPlot(Seurat_object_subset, Stress, reduction.use = "umap", min.color = "lightgrey", max.color = "blue")
VlnPlot(Seurat_object_subset, Stress, pt.size = 0, group.by = "seurat_clusters", ncol = 2)
# Cell-type contaminants plots, such as check for subcortical interneuron and pyramidal cells - these are to exclude
Contamin <- c("Meis2", "Nkx2-1")
dittoDimPlot(Seurat_object_subset, Contamin, reduction.use = "umap", min.color = "lightgrey", max.color = "blue", split.by = "sample")
VlnPlot(Seurat_object_subset, Contamin, pt.size = 0, group.by = "seurat_clusters", cols = dittoColors())
# Pyramidal cells markers (UMAP + Violin plot)
Pyramid <- c("Slc17a7", "Neurod6", "Tbr1")
dittoDimPlot(Seurat_object_subset, Pyramid, reduction.use = "umap", min.color = "lightgrey", max.color = "blue")
VlnPlot(Seurat_object_subset, Pyramid, pt.size = 0, group.by = "seurat_clusters", ncol = 2)
#Positive markers plots
Posmarkers <- c("Sst", "Lhx6", "Gad1", "Gad2")
dittoDimPlot(Seurat_object_subset, Posmarkers, reduction.use = "umap", min.color = "lightgrey", max.color = "blue")
VlnPlot(Seurat_object_subset, Posmarkers, pt.size = 0, group.by = "seurat_clusters", ncol = 2)
#use plot information to decide manually which clusters you want to keep 
# Remove contaminants 
# Specify clusters you want to remove
Seurat_object_cleaned <- subset(Seurat_object_subset, idents = c('8','12', '13', '14'), invert = T)
rm(Seurat_object_subset)
gc()

############################
### Recluster selected cells
############################

Seurat_object_cleaned <- NormalizeData(Seurat_object_cleaned)
Seurat_object_cleaned <- FindVariableFeatures(Seurat_object_cleaned,nfeatures = 2000)
Seurat_object_cleaned <- ScaleData(Seurat_object_cleaned,vars.to.regress = c("nFeature_RNA",'percent.mt','ccDiff'), model.use = 'linear', block.size=dim(Seurat_object_cleaned)[1])
Seurat_object_cleaned<-RunPCA(Seurat_object_cleaned, npcs = 60)
ElbowPlot(Seurat_object_subset, ndims = 60)
data.use.Seurat_object_cleaned <- PrepDR(object = Seurat_object_cleaned, genes.use = VariableFeatures(object = Seurat_object_cleaned), use.imputed = F, assay.type = "RNA")
path_data <- getwd()
nPCs.5_clean <- PCA_estimate_nPC(data.use.Seurat_object_cleaned, whereto = paste0(path_data,"/optimal_nPCs_5_clean.RDS"), by.nPC = 5, from.nPC = 30, to.nPC = 40) 
pdf('nPCs_estimate_clean.pdf')
#adjust from.nPC = XXX, to.nPC = XXX based on ElbowPlot ndims. You can also change the number of k
nPCs.clean <- PCA_estimate_nPC(data.use.Seurat_object_cleaned, whereto = paste0(path_data,"/optimal_nPCs_clean.RDS"), by.nPC = 1, from.nPC = nPCs.5_clean-5, to.nPC = nPCs.5_clean + 5) 
dev.off()
Seurat_object_cleaned <- RunUMAP(Seurat_object_cleaned, dims = 1:nPCs.clean)
Seurat_object_cleaned <- FindNeighbors(Seurat_object_cleaned, dims = 1:nPCs.clean)
Seurat_object_cleaned <- FindClusters(Seurat_object_cleaned, resolution = 0.8)
# plot clusters
dittoDimPlot(Seurat_object_cleaned, "seurat_clusters", reduction.use = "umap", do.label = T)
# Check again the expression of gene of interest
dittoDimPlot(Seurat_object_cleaned, "Sst", reduction.use = "umap", min.color = "lightgrey", max.color = "blue")
VlnPlot(Seurat_object_cleaned, "Sst", pt.size = 0, group.by = "seurat_clusters", cols = dittoColors())
# Quality control plots to asses different matrices
VlnPlot(Seurat_object_cleaned, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.ribo"), ncol = 2, pt.size = 0, group.by = "seurat_clusters")
# If contamination or unwanted clusters remain, this code can be re-run to clean the data further and perform reclustering.
#Save final seurat object
saveRDS(Seurat_object_cleaned, file = "Dataset_P1_Dibella.rds") 






