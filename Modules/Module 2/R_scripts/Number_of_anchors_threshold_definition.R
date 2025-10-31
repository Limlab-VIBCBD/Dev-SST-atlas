# This is the code used to estimate the threshold for the minimum number of anchors required. 
# Atlas-v0 is composed of four samples collected at different developmental stages, we computed anchors across all pairs of samples in Atlas-v0. 
# For each of the four samples in Atlas-v0, 1,000 cells were randomly subsampled, and anchors were identified with the remaining three samples. 
# The anchor values were compiled into a distribution, and the 25th percentile was chosen as the acceptance criterion.

## Installation needed
if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(
  BiocManager, Seurat,dplyr, R.utils
  gridExtra, gitcreds
)

########################
# Load required packages
########################
library(Seurat)
library(dplyr)
library(R.utils)


#########################################
# Compute anchors across Atlas-v0 samples
#########################################

#############
# E16 sample
#############
ref_atlas$sample <- paste0("BaseAtlas_",ref_atlas$orig.ident)
DefaultAssay(ref_atlas)<-"RNA"
ref_atlas_split<-SplitObject(ref_atlas, split.by = "sample")
# subsamples are repeated to ensure a complete representation of the sample
cells_sub<-list()
for(i in 1:4){
  cells_sub[[i]]<-sample(colnames(ref_atlas_split[["BaseAtlas_E16"]]),size = 1000,replace = FALSE)
}
number_of_anchors_sub_E16 <- lapply(cells_sub, function(x){
  ref_atlas_split<-SplitObject(ref_atlas, split.by = "sample")
  ref_atlas_split[["BaseAtlas_E16"]] <- subset(ref_atlas_split[["BaseAtlas_E16"]], cells=x)
  for(sample in names(ref_atlas_split)){
    obj<-ref_atlas_split[[sample]]
    obj<-NormalizeData(obj)
    obj<-FindVariableFeatures(obj,nfeatures = 3000)
    ref_atlas_split[[sample]]<-obj
  }
  i.anchors_BaseAtlas_E16 <- FindIntegrationAnchors(object.list = ref_atlas_split, dims = 1:30,reduction='cca',scale=T,k.anchor=5,k.filter=100,k.score=15,anchor.features=3000, reference = c(1))
  number_of_anchors <- table(i.anchors_BaseAtlas_E16@anchors$dataset1)[2:4]
  names(number_of_anchors) <- names(ref_atlas_split)[2:4]
  return(number_of_anchors)
})
number_of_anchors_sub_E16<-bind_rows(number_of_anchors_sub_E16)
number_of_anchors_sub_E16 <- colMeans(number_of_anchors_sub_E16)


############
# P1 sample
############
ref_atlas$sample <- paste0("BaseAtlas_",ref_atlas$orig.ident)
DefaultAssay(ref_atlas)<-"RNA"
ref_atlas_split<-SplitObject(ref_atlas, split.by = "sample")
# subsamples are repeated to ensure a complete representation of the sample
cells_sub<-list()
for(i in 1:3){
  cells_sub[[i]]<-sample(colnames(ref_atlas_split[["BaseAtlas_P1"]]),size = 1000,replace = FALSE)
}
number_of_anchors_sub_P1 <- lapply(cells_sub, function(x){
  ref_atlas_split<-SplitObject(ref_atlas, split.by = "sample")
  ref_atlas_split[["BaseAtlas_P1"]] <- subset(ref_atlas_split[["BaseAtlas_P1"]], cells=x)
  for(sample in names(ref_atlas_split)){
    obj<-ref_atlas_split[[sample]]
    obj<-NormalizeData(obj)
    obj<-FindVariableFeatures(obj,nfeatures = 3000)
    ref_atlas_split[[sample]]<-obj
  }
  i.anchors_BaseAtlas_P1 <- FindIntegrationAnchors(object.list = ref_atlas_split, dims = 1:30,reduction='cca',scale=T,k.anchor=5,k.filter=100,k.score=15,anchor.features=3000, reference = c(2))
  number_of_anchors <- table(i.anchors_BaseAtlas_P1@anchors$dataset1)[c(1,3,4)]
  names(number_of_anchors) <- names(ref_atlas_split)[c(1,3,4)]
  return(number_of_anchors)
})
number_of_anchors_sub_P1<-bind_rows(number_of_anchors_sub_P1)
number_of_anchors_sub_P1 <- colMeans(number_of_anchors_sub_P1)


############
# P5 sample
############
ref_atlas$sample <- paste0("BaseAtlas_",ref_atlas$orig.ident)
DefaultAssay(ref_atlas)<-"RNA"
ref_atlas_split<-SplitObject(ref_atlas, split.by = "sample")
# subsamples are repeated to ensure a complete representation of the sample
cells_sub<-list()
for(i in 1:3){
  cells_sub[[i]]<-sample(colnames(ref_atlas_split[["BaseAtlas_P5"]]),size = 1000,replace = FALSE)
}
number_of_anchors_sub_P5 <- lapply(cells_sub, function(x){
  ref_atlas_split<-SplitObject(ref_atlas, split.by = "sample")
  ref_atlas_split[["BaseAtlas_P5"]] <- subset(ref_atlas_split[["BaseAtlas_P5"]], cells=x)
  for(sample in names(ref_atlas_split)){
    obj<-ref_atlas_split[[sample]]
    obj<-NormalizeData(obj)
    obj<-FindVariableFeatures(obj,nfeatures = 3000)
    ref_atlas_split[[sample]]<-obj
  }
  i.anchors_BaseAtlas_P5 <- FindIntegrationAnchors(object.list = ref_atlas_split, dims = 1:30,reduction='cca',scale=T,k.anchor=5,k.filter=100,k.score=15,anchor.features=3000, reference = c(3))
  number_of_anchors <- table(i.anchors_BaseAtlas_P5@anchors$dataset1)[c(1,2,4)]
  names(number_of_anchors) <- names(ref_atlas_split)[c(1,2,4)]
  return(number_of_anchors)
})
number_of_anchors_sub_P5<-bind_rows(number_of_anchors_sub_P5)
number_of_anchors_sub_P5 <- colMeans(number_of_anchors_sub_P5)


########################
# P5 sample fixed sorted
########################
ref_atlas$sample <- paste0("BaseAtlas_",ref_atlas$orig.ident)
DefaultAssay(ref_atlas)<-"RNA"
ref_atlas_split<-SplitObject(ref_atlas, split.by = "sample")

cells_sub<-list()
for(i in 1:10){
  cells_sub[[i]]<-sample(colnames(ref_atlas_split[["BaseAtlas_lim_P5_fixed_sorted"]]),size = 1000,replace = FALSE)
}
# subsamples are repeated to ensure a complete representation of the sample
number_of_anchors_sub_P5_fixed_sorted <- lapply(cells_sub, function(x){
  ref_atlas_split<-SplitObject(ref_atlas, split.by = "sample")
  ref_atlas_split[["BaseAtlas_lim_P5_fixed_sorted"]] <- subset(ref_atlas_split[["BaseAtlas_lim_P5_fixed_sorted"]], cells=x)
  for(sample in names(ref_atlas_split)){
    obj<-ref_atlas_split[[sample]]
    obj<-NormalizeData(obj)
    obj<-FindVariableFeatures(obj,nfeatures = 3000)
    ref_atlas_split[[sample]]<-obj
  }
  i.anchors_BaseAtlas_P5_fixed_sorted <- FindIntegrationAnchors(object.list = ref_atlas_split, dims = 1:30,reduction='cca',scale=T,k.anchor=5,k.filter=100,k.score=15,anchor.features=3000, reference = c(4))
  number_of_anchors <- table(i.anchors_BaseAtlas_P5_fixed_sorted@anchors$dataset1)[c(1,2,3)]
  names(number_of_anchors) <- names(ref_atlas_split)[c(1,2,3)]
  return(number_of_anchors)
})
number_of_anchors_sub_P5_fixed_sorted<-bind_rows(number_of_anchors_sub_P5_fixed_sorted)
number_of_anchors_sub_P5_fixed_sorted <- colMeans(number_of_anchors_sub_P5_fixed_sorted)

number_of_anchors<-c(mean(number_of_anchors_sub_P5_fixed_sorted), mean(number_of_anchors_sub_P5), mean(number_of_anchors_sub_P1),mean(number_of_anchors_sub_E16))
# we extract the final threshold
quantile(number_of_anchors,0.25)

