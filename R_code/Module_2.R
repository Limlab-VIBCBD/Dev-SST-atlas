
## Installation needed
if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(
  BiocManager, Seurat, dittoSeq, ggplot2, gtools, dplyr, R.utils, dbscan, EnvStats,
  gridExtra, gitcreds
)
#To install kBET you need a GitHub token. Generate a classic token via this link https://github.com/settings/tokens/new and insert you personal token in the following code:
#I also found this: 
#usethis::create_github_token()
if (!requireNamespace("kBET", quietly = TRUE)) {
  Sys.setenv(GITHUB_PAT = "your_token")
  devtools::install_github("theislab/kBET")
}
# donwnload and install Rtools from here https://cloud.r-project.org/bin/windows/Rtools/


########################
# Load required packages
########################
library(Seurat)
library(dittoSeq)
library(ggplot2)
library(gtools)
library(dplyr)
library(R.utils)
library(dbscan)
library(kBET)
library(EnvStats)
'%!in%' <- function(x,y)!('%in%'(x,y))
# Source the external script 'Seurat_Utils' for further analysis (10-fold Singular Value Decomposition (SVD) cross validation to predict number of PCs)
source("Seurat_Utils.R") 

###########
# Load data
###########

## Upload Base atlas and dataset to evaluate
# Set the working directory and open the dataset of interest
test_dataset_filtered_sst <- readRDS("P1_DS24_Clean.rds")
# Add a column to the metadata and specify the name of testing sample 
test_dataset_filtered_sst$sample <- "test_dataset"
# Set the working directory and open the reference Dataset
ref_atlas <- readRDS("Atlas_v0.rds")
ref_atlas$sample <- paste0("BaseAtlas_",ref_atlas$orig.ident)

DefaultAssay(ref_atlas) <- "RNA"
DefaultAssay(test_dataset_filtered_sst) <- "RNA"
# If your reference dataset is composed of multiple samples, split them to compare your Dataset of interest with each of the samples in the reference Dataset
ref_atlas_split <- SplitObject(ref_atlas, split.by = "sample")


###########################
# COMPUTE NUMBER OF ANCHORS
###########################

# subsample 1'000 cells of the test dataset
# subsamples are repeated to ensure a complete representation of the sample
cells_sub<-list()
for(i in 1:round(dim(test_dataset_filtered_sst)[2]/1000+0.5)){
  set.seed(345+i)
  cells_sub[[i]]<-sample(colnames(test_dataset_filtered_sst),size = 1000,replace = FALSE)
}
# compute number of anchors between each subsample of testing datset and the reference dataset
number_of_anchors_sub <- lapply(cells_sub, function(x){
  test_dataset <- subset(test_dataset_filtered_sst, cells=x)
  test_dataset<-SplitObject(test_dataset, split.by = "sample")
  integrated_all_list <- c(ref_atlas_split,test_dataset)
  for(sample in names(integrated_all_list)){
    obj<-integrated_all_list[[sample]]
    obj<-NormalizeData(obj)
    obj<-FindVariableFeatures(obj,nfeatures = 3000)
    integrated_all_list[[sample]]<-obj
  }
  i.anchors <- FindIntegrationAnchors(object.list = integrated_all_list, dims = 1:30,reduction='cca',scale=T,k.anchor=5,k.filter=100,k.score=15,anchor.features=3000, reference = c(5))
  number_of_anchors <- table(i.anchors@anchors$dataset1)[1:4]
  names(number_of_anchors) <- names(ref_atlas_split)
  return(number_of_anchors)
})
number_of_anchors_sub<-bind_rows(number_of_anchors_sub)
number_of_anchors <- colMeans(number_of_anchors_sub)
# Print number of anchors and write table
number_of_anchors
write.table(number_of_anchors, "Number_of_anchors_withBaseAtlas.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
#save results
saveRDS(number_of_anchors,"number_of_anchors_with_Reference atlas_P1_DS24.RDS")


#################################
# Perform CCA integration and PCA
#################################
test_dataset_filtered_sst_split<-SplitObject(test_dataset_filtered_sst, split.by = "sample")
integrated_all_list <- c(ref_atlas_split,test_dataset_filtered_sst_split)
for(sample in names(integrated_all_list)){
  obj<-integrated_all_list[[sample]]
  obj<-NormalizeData(obj)
  obj<-FindVariableFeatures(obj,nfeatures = 3000)
  integrated_all_list[[sample]]<-obj
}
# perform CCA integration of testing dataset and reference dataset
i.anchors <- FindIntegrationAnchors(object.list = integrated_all_list, dims = 1:30,reduction='cca',scale=T,k.anchor=5,k.filter=100,k.score=15,anchor.features=3000)
integrated_data <- IntegrateData(anchorset = i.anchors, dims = 1:30, normalization.method ='LogNormalize', k.weight=100)
DefaultAssay(integrated_data)<-'integrated'
integrated_data <- ScaleData(integrated_data, verbose = FALSE)
integrated_data <- RunPCA(integrated_data,npcs = 100)
# Determine the optimal number of principal components (PCs) for downstream analysis
ElbowPlot(Seurat_object_subset, ndims = 60)
# Adjusting nPC based on ndims on the ElbowPlot 
data.use.integrated<- PrepDR(object = integrated_data, genes.use = VariableFeatures(object = integrated_data), use.imputed = F, assay.type = "integrated")
path_data <- getwd()
# Adjust `from.nPC` and `to.nPC` after inspecting the ElbowPlot
nPCs.data.use <- PCA_estimate_nPC(data.use.integrated, whereto=paste0(path_data,"/optimal_nPCs_integrated.RDS"), by.nPC=1, from.nPC = 30,to.nPC = 40)
# save the most informative principal components
pca_integrated<-Embeddings(integrated_data, reduction = "pca")[,1:nPCs.data.use]


###################################
# Evaluate Neighborhood Composition
###################################
cells_test_dataset <- colnames(subset(integrated_data, subset=sample %in% 'test_dataset'))
cells_ref <- colnames(subset(integrated_data, subset=sample %!in% 'test_dataset'))
# Compute a k-nearest neighbors graph on integrated PCs (pca_integrated), for increasing number of k: seq(5,200,5)
# compute, for each cell in the testing dataset, the fraction of cells from the reference atlas among its k-nearest neighbors (kNN)
# report for each k, the fraction of cells in testing dataset having at least half of the global fraction of reference dataset cells among their kNN
percentage_nn <- c()
# Global fraction of cells from refence atlas
fraction_ref_in_integrated <- length(cells_ref)/(length(cells_ref)+length(cells_test_dataset))
for(k_score in seq(5,200,5)){
  knn_integrated <- dbscan::kNN(x = pca_integrated %>% as.matrix(), k = k_score, query = pca_integrated[cells_test_dataset,])
  knn_integrated.data <- data.frame(from = rep(rownames(knn_integrated$id), k_score),
                                    to = rownames(pca_integrated)[as.vector(knn_integrated$id)])
  knn_integrated.data <- knn_integrated.data[which(knn_integrated.data$to %in% cells_ref),]
  percentage_nn<-c(percentage_nn, sum(table(knn_integrated.data$from) > k_score*fraction_ref_in_integrated/2)/length(cells_test_dataset))
}
names(percentage_nn) <- seq(5,200,5)
# save outputs
saveRDS(percentage_nn,file = paste0("Percentages_of_cells_in_having_at_least_halfExpectedCells_from_BaseAtlas_in_kNN.RDS"))
write.table(percentage_nn["60"], paste0("Percentage_of_cells_in_having_at_least_halfExpectedCells_from_BaseAtlas_in_60NN.txt"),sep = "\t", quote = FALSE, col.names = FALSE, row.names=FALSE)
# print plot of results
pdf("Percentages_of_cells.pdf")
plot(seq(5,200,5),percentage_nn, pch=19, xlab="k",ylab="percentage of dataset", ylim=c(0,1))
abline(v=60, col="red2")
dev.off()


###################
# Perform kBET test
###################
batch <- integrated_data$sample
batch[which(batch %!in% c("test_dataset"))] <- "ref"
batch[which(batch %in% c("test_dataset"))] <- "new"
# subset 1'000 cells from each sample of reference dataset and testing dataset
set.seed(1234)
subset_ref_E16 <- sample(rownames(subset(integrated_data@meta.data, subset = sample == "BaseAtlas_E16")), size = 1000, replace=FALSE)
subset_ref_P1 <- sample(rownames(subset(integrated_data@meta.data, subset = sample == "BaseAtlas_P1")), size = 1000, replace=FALSE)
subset_ref_P5 <- sample(rownames(subset(integrated_data@meta.data, subset = sample == "BaseAtlas_P5")), size = 1000, replace=FALSE)
subset_ref_P5_fs <- sample(rownames(subset(integrated_data@meta.data, subset = sample == "BaseAtlas_lim_P5_fixed_sorted")), size = 1000, replace=FALSE)
cells_ref_subset <- c(subset_ref_E16,subset_ref_P1,subset_ref_P5, subset_ref_P5_fs)
cells_test_dataset_subset <- sample(cells_test_dataset, size = 1000, replace=FALSE)
pca_integrated_subset <- pca_integrated[c(cells_ref_subset,cells_test_dataset_subset),]
# Perform kBET analysis and plot results
pdf(paste0('kBET_results_boxplot.pdf'))
batch.estimate <- kBET(pca_integrated_subset, batch[rownames(pca_integrated_subset)], do.pca = FALSE, plot = TRUE , testSize = round(dim(pca_integrated_subset)[1] * 0.1), k0 = 20, n_repeat = 20)
dev.off()
# you can get stats for plotting from batch.estimate$stats
# save results
saveRDS(batch.estimate, file=paste0("kBET_batch_estimate.RDS"))
write.table(mean(batch.estimate$stats$kBET.observed), paste0("MeanRejectionRate_kBET_batch_estimate.txt"), sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)

