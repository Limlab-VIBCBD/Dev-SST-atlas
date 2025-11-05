## THIS IS AN ALTERNATIVE TO NEXTFLOW PIPELINE

########################
# Load required packages
########################
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
#########################################
# Load functions for iterative clustering
#########################################
source("Iterative_Clustering_pipeline_functions.R")

# set the working directory to the location where you want the output files to be saved.
setwd('output_dir/')
cr.dir<-getwd()
# load filtered Seurat object 
obj <- readRDS('Integrated_atlas_v2_NonMartinotti.RDS')
### create environment variable to save final clustering (after all iterations)
final.clusters<-rep("c", dim(obj)[2])
names(final.clusters) <- rownames(obj@meta.data)
### launch pipeline
# 'sample' — name of the column in the metadata of 'obj' containing sample IDs  
# nPC — number of PCs to use; if set to NULL, the optimal number is automatically estimated  
# min.res and max.res — range of resolutions to use when computing clusters  
# perc.sub — percentage of cells to subsample from the whole dataset for cluster stability evaluation  
# n_subsampling — number of subsampling iterations to perform for cluster stability evaluation  
# jaccard_cutoff — Jaccard index threshold used to define cluster stability  
# percent_cutoff — percentage of subsamplings that must meet the jaccard_cutoff to classify a cluster as stable  
# minSize — minimum cluster size  
# DEscore.cutoff — minimum DEscore value 
###################################
# RUN iterative clustering pipeline
###################################
callIteration(obj, 'batch', nPC=NULL, min.res=0.1, max.res=0.5, perc.sub=0.8, n_subsampling=20, jaccard_cutoff=0.75, percent_cutoff=0.749, minSize=149, DEscore.cutoff=60,cr.dir=cr.dir)
# save final clusters produced by the pipeline
saveRDS(final.clusters,paste0(cr.dir,"/finalClusters.RDS"))