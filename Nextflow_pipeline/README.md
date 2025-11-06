# Nextflow pipeline for iterative clustering 

This nextflow pipeline run the iterative clustering algorithm used to compute clusters for the DEV-Sst-Atlas. 

## Alghorithm
Iterative clustering is performed through the following steps:

+ Integration and clustering. Perform integration with CCA (samples with less than 20 cells are removed) and compute clusters (optimal number of PCs is estimated) for a limited set of resolutions (from 0.1 to 0.5). The tested resolution range is limited to prevent the formation of an excessively large number of clusters. Instead, we use an iterative approach, in which each cluster identified at one step is further subdivided in subsequent rounds of analysis.

+ Evaluate cluster stability. Subsample 80% of the dataset multiple times and recalculate clusters for each resolution. Compute jaccard index between clusters obtained on subsamplings and full dataset. Clusters with jaccard index higher than 0.75 in at least 75% of the subsamples are considered stable.

+ Select the resolution producing the most stable set of clusters. Choose the resolution that produces the clustering with the highest proportion of cells belonging to stable clusters. If multiple resolutions meet this criterion, select the one yielding the greater number of clusters. This selection is limited to clusterings where at least 70% of cells belong to stable clusters and at least 60% of clusters are stable. If no clustering meets these cutoffs, the iteration stops.

+ Check cluster size. Merge clusters containing fewer than 150 (or 300) cells with their nearest cluster. Distances between clusters are calculated based on centroids in PCA space using Euclidean distance.

+ Check cluster identity. Identify differentially expressed genes (DEGs) by comparing each cluster to all other cells (genes with avg_log2FC > 1 and p_adj < 0.01 are selected). Compute the DEscore (based on p_adj values of selected markers) and merge any cluster with a DEscore < 60 with its closest cluster. 

+ Replicate iteration. If only one cluster remains after completing all previous steps, the iteration stops. Otherwise, repeat the full procedure for all clusters until a single cluster is obtained.

## Requirements
To run the pipeline, you need to have Nextflow and [Conda](https://docs.conda.io/projects/conda/en/stable/user-guide/install/index.html) or Miniconda package manager installed on your system. \
The Nextflow pipeline will automatically create and activate a Conda environment with all required software packages and their dependecies. By default, Nextflow instructs Conda to save the required environments in the pipeline work/ directory.

## Input files
The pipeline requires as input a saved R Seurat object containing in the meta.data a column named ``batch`` with information about sample_id.\
Seurat object must be compatible with Seurat v4. If you have a v5 assy you can convert it using the following command:
```R
# convert a v5 assay to a v3 assay
input_object[["RNA3"]] <- as(object = input_object[["RNA"]], Class = "Assay")
```

## Config file
You can configure the input file path, adjust algorithm parameters, specify the workload manager used on your system (SLURM is set as the default), and modify the resource allocations for each process in the ``nextflow.config`` file.
We recommend increasing the memory allocation when working with datasets containing more than 50,000 cells.
## Run pipeline
To execute the pipeline, you can clone the repository and run the following command from within the pipeline directory: 
```nextflow
nextflow run clustering_pipeline.nf
```
If the pipeline stops for any reason, you can resume the run from the last successful process by adding the -resume option:
```nextflow
nextflow run clustering_pipeline.nf -resume
```
## Outputs
Once the pipeline has completed successfully, you will find the results in the specified output folder (by default, Nextflow creates a results/ folder in the pipeline directory).
This folder should contain a Seurat object named ``object_with_final_clusters.RDS``, which includes the final cluster assignments. These assignments are also saved in a separate file within the same folder.
You can plot the UMAP with clusters using the following code:
```R
library(Seurat)
obj<-readRDS('object_with_final_clusters.RDS')
DimPlot(obj, group.by="final_clusters')
```