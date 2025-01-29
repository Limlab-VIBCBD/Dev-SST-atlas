library(stringr)
library(Seurat)
library(dittoSeq)
library(ggplot2)
library(gtools)
library(dplyr)
library(dittoSeq)
library(R.utils)
library(colorspace)
library(scclusteval)
library(tidyr)
library(purrr)
library(parallel)
library(doMC)
source("Seurat_Utils.R")
source("Clustering_pipeline.r")

cr.dir<-getwd()
obj <- readRDS('Atlas_v2.RDS') # load seurat obj to compute clustering 
### create environment variable to save final clustering (after all iterations)
final.clusters<-rep("c", dim(obj)[2])
names(final.clusters) <- rownames(obj@meta.data)
final.clusters[names(clusters)] <- paste(final.clusters[names(clusters)],clusters,sep=".")
### launch pipeline
lapply(names(table(clusters)[table(clusters) >= 300]),function(cl){
  Iteration <<- Iteration + 1
  print(paste0("Iteration ",Iteration,"- Cluster ",cl))
  obj <- subset(obj, cells=names(clusters[clusters %in% cl]))
  setwd(cr.dir)
  dir.create(paste0("Cluster_",cl))
  setwd(paste0("Cluster_",cl))
  cr.dir<-paste0(cr.dir,"/Cluster_",cl)
  callIteration(obj, batch='sample', nPC=NULL, min.res=0.1, max.res=0.5, perc.sub=0.8, n_subsampling=20, jaccard_cutoff=0.75, percent_cutoff=0.74, minSize=149, DEscore.cutoff=60,cr.dir=cr.dir)
  gc()
})
#saveRDS(final.clusters,paste0(cr.dir,"/finalClusters.RDS")) #save resulting clusterings after all iterations end

