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
source("Seurat_Utils.R")
'%!in%' <- function(x,y)!('%in%'(x,y))

### Compute integration of samples (batch) using cca and compute clusters
### clustering is computed for all resolution between min.res and max.res
### return processed seurat object
computeIntegrationAndClusters <- function(obj, batch, nPCs=NULL, min.res=0.1, max.res=1){
  DefaultAssay(obj) <- 'RNA'
  remove_samples<-names(table(obj@meta.data[,batch]))[table(obj@meta.data[,batch]) < 20]
  if(length(remove_samples > 0)){
    eval(parse(text=paste0("obj<-subset(obj, subset=",batch," %!in% remove_samples)")))
  }
  integrated_all_list <- SplitObject(obj, split.by = batch)
  integrated_all_list <- lapply(integrated_all_list, function(sample) {
    sample<-NormalizeData(sample, verbose = FALSE)
    sample<-FindVariableFeatures(sample,nfeatures = 2000, verbose = FALSE)
    return(sample)
  })
  i.anchors <- FindIntegrationAnchors(object.list = integrated_all_list, dims = 1:min(min(table(obj@meta.data[,batch]))-5,30),reduction='cca',scale=T,k.anchor=5,k.filter=100,k.score=15,anchor.features=3000,verbose = FALSE)
  obj <- IntegrateData(anchorset = i.anchors, dims = 1:min(min(table(obj@meta.data[,batch]))-5,30), normalization.method ='LogNormalize', k.weight=min(100,min(table(obj@meta.data[,batch]))),verbose = FALSE)
  DefaultAssay(obj)<-'integrated'
  obj<-ScaleData(obj,vars.to.regress = c("nFeature_RNA",'percent.mt','ccDiff'), verbose = FALSE)
  obj<-RunPCA(obj,npcs = 60, verbose = FALSE)
  if(is.null(nPCs)){
    data.use.integrated<- PrepDR(object = obj, genes.use = VariableFeatures(object = obj), use.imputed = F, assay.type = "integrated")
    path_data <- getwd()
    nPCs.data.use5 <- PCA_estimate_nPC(data.use.integrated, whereto=paste0(path_data,"/optimal_nPCs_5_integrated.RDS"), by.nPC=5,from.nPC = 10,to.nPC = 55) 
    nPCs <- PCA_estimate_nPC(data.use.integrated, whereto=paste0(path_data,"/optimal_nPCs_integrated.RDS"), by.nPC=1,from.nPC = nPCs.data.use5-5,to.nPC = nPCs.data.use5+5) 
  }
  obj<-RunUMAP(obj,dims=1:nPCs, verbose = FALSE)
  obj<-FindNeighbors(obj,dims=1:nPCs, verbose = FALSE)
  resolutions_selected <- seq(min.res, max.res, 0.1)
  obj<-FindClusters(obj,resolution = resolutions_selected, verbose = FALSE)
  saveRDS(obj, "obj_with_clustering.RDS")
  return(obj)
}

### check if clustering match at different resolutions (we don't want to run cluster stability analysis on same set of clusters more than once)
### returns resolutions on which to perform cluster stability analysis
findResolutions <- function(obj, min.res=0.1, max.res=1){
  resolutions_selected <- seq(min.res, max.res, 0.1)
  resolutions <- paste0("integrated_snn_res.",resolutions_selected)
  n_tmp <- length(resolutions)
  resolutions <- resolutions[resolutions %in% colnames(obj@meta.data)]
  if(length(resolutions) == 0){
    stop('Error: no requested resolution found in Seurat object')
  }
  n_clusters <- unlist(apply(obj@meta.data[,resolutions], 2, function(x) length(unique(x))))
  dup<-n_clusters[duplicated(n_clusters)]
  if(length(dup)==0){
    resolution_for_stability <- resolutions
  }else{
    resolution_for_stability<-names(n_clusters[!(n_clusters %in% dup)])
    resolution_for_stability<-c(resolution_for_stability, unlist(lapply(split(n_clusters[(n_clusters %in% dup)], as.factor(n_clusters[(n_clusters %in% dup)])), function(x){
      res<-c()
      for(k in 1:(length(x)-1)){
        n_match <- 0
        for (j in (k+1):length(x)){
          perc_match <- c()
          for (i in 0:(x[k]-1)){
            cells_in_cluster<-rownames(obj@meta.data[which(obj@meta.data[,names(x)[k]]==i),])
            perc_match <- c(perc_match, max(table(obj@meta.data[cells_in_cluster,names(x)[j]])/length(cells_in_cluster)))
          }
          if(min(perc_match) < 0.98){
            # no match
            n_match <- n_match +1
          }
        }
        if(n_match == (length(x)-k)){
          res<-c(res,names(x)[k])
        }
      }
      res<-c(res,names(x)[length(x)])
      return(res)
    })))
  }
  return(resolution_for_stability)
}


### Sample seurat obj extracting perc.sub of all cells and compute new clustering. 
### Subsampling is performed n_subsampling times
### Return tibble with clusters for all subsampling and all resolutions
performSubsamplingAndClustering <- function(obj, resolutions, nPCs, perc.sub=0.8, n_subsampling=20){
  subsample_idents<-Reduce( "bind_rows", lapply(resolutions, function(x){
    res<-str_remove(x,"integrated_snn_res.")
    out<-mclapply(c(1:n_subsampling), function(y){
      rand_test <- RandomSubsetData(obj, rate=perc.sub)
      rand_test <- ScaleData(rand_test,vars.to.regress = c("nFeature_RNA",'percent.mt','ccDiff'), verbose = FALSE)
      rand_test <- RunPCA(rand_test, verbose = FALSE,npcs = 60)
      rand_test <- FindNeighbors(rand_test, reduction = "pca", dims = 1:nPCs, verbose = FALSE)
      eval(parse(text=paste0("rand_test <- FindClusters(rand_test, resolution = ",res,", verbose = FALSE)")))
      recluster_ident<-rand_test@meta.data[,x]
      names(recluster_ident) <- rownames(rand_test@meta.data)
      original_ident<-obj@meta.data[rownames(rand_test@meta.data),x]
      names(original_ident) <- rownames(rand_test@meta.data)
      return(list(recluster_ident,original_ident))
    },mc.cores = 11)
    recluster_ident <-list()
    original_ident <-list()
    for(i in 1:n_subsampling){
      recluster_ident[[i]]=out[[i]][[1]]
      original_ident[[i]]=out[[i]][[2]]
    }
    data<-tibble(resolution=res, recluster_ident=recluster_ident, original_ident=original_ident,round=as.character(c(0:(n_subsampling-1))))
    return(data)
  }))
  return(subsample_idents)
}

### Return tibble with clusters for the whole dataset for all resolutions
prepareFullData <- function(obj, resolutions){
  original_ident_full<-lapply(resolutions, function(x){
    original_ident_full <- obj@meta.data[,x]
    names(original_ident_full) <- rownames(obj@meta.data)
    return(original_ident_full)
  })
  res<-unlist(lapply(resolutions, function(x){
    res<-str_remove(x,"integrated_snn_res.")
    return(res)
  }))
  fullsample_idents <- tibble(resolution=res,original_ident_full=original_ident_full)
  return(fullsample_idents)
}

#### Find most stable clustering.
### Stable clusters are defined as clusters that have jaccard higher than jaccard_cutoff in at least percent_cutoff of all subsampling
### Percentage of cell in stable cluster is computed.
### Sets of clusters with percentage of stable cells higher than 70% are selected. The clustering with highest percentage of stable cells and highest number of clusters is selected.
findStableClustering <- function(subsample_idents, fullsample_idents, jaccard_cutoff=0.8, percent_cutoff=0.8){
  subsample_idents_list<- subsample_idents %>% group_by(resolution) %>%  nest()
  stable_clusters<- subsample_idents_list %>% mutate(stable_cluster = map(data, ~ AssignStableCluster(.x$original_ident,.x$recluster_ident,jaccard_cutoff = jaccard_cutoff,method = "jaccard_percent", percent_cutoff = percent_cutoff)))
  df <- dplyr::left_join(stable_clusters, fullsample_idents) %>% dplyr::ungroup() %>% dplyr::mutate(total = map_dbl(stable_cluster, ~length(.x$stable_cluster))) %>% 
    dplyr::mutate(stable = map_dbl(stable_cluster, ~.x$number_of_stable_cluster)) %>% 
    dplyr::mutate(percentage_cluster = map2_dbl(stable, total, function(x, y) x/y)) %>%
    dplyr::mutate(percentage = map2_dbl(original_ident_full, stable_cluster, function(x, y) CalculatePercentCellInStable(x,y$stable_cluster))) %>% 
    dplyr::select(-data, -stable_cluster, -original_ident_full) %>% dplyr::mutate_if(is.character, function(x) as.factor(as.numeric(x))) %>% 
    tidyr::gather(total:stable,key = "category", value = "number")
  
  # check stability and choose best clustering
  df<-subset(df, subset=percentage > 0.7 & percentage_cluster >= 0.6 & category== "total")
  if (nrow(df) > 0){
    df<-subset(df, subset=percentage == max(df$percentage))
    df<-subset(df, subset=number == max(df$number))
    stable_resolution <- df$resolution[1]
  }else{
    stable_resolution <- NULL
  }
  return(stable_resolution)
}

### PassClusters is a logical vector with length equal to number of clusters 
### clusters with FALSE values in PassClusters vector are merged with their closest clusters in PCA space 
### centroids are computed for each cluster and euclidean distances between centroids are computed.  
mergeClosestCluster <-function(obj, clusters, PassClusters, nPCs){
  pca <- Embeddings(obj, reduction = "pca")[,1:nPCs]
  centroids<-Reduce("rbind", lapply(unique(clusters), function(x){
    cells<-names(clusters[which(clusters %in% x)])
    centroid<-colMeans(pca[cells,])
    return(centroid)
  }))
  rownames(centroids) <- unique(clusters)
  for (i in names(PassClusters[which(PassClusters == FALSE)])){
    if(i %in% clusters){
      euclidean_dist<-as.matrix(dist(centroids, method="euclidean", diag=TRUE, upper=TRUE))
      rownames(euclidean_dist)<-rownames(centroids)
      colnames(euclidean_dist)<-rownames(centroids)
      diag(euclidean_dist) <- Inf
      closest_cluster<-colnames(euclidean_dist)[euclidean_dist[i,] == min(euclidean_dist[i,])]
      cells_1<-names(clusters[which(clusters %in% i)])
      cells_2<-names(clusters[which(clusters %in% closest_cluster)])
      clusters[c(cells_1,cells_2)] <- paste(i,closest_cluster,sep="_")
      n_clusters <- length(unique(clusters))
      if(n_clusters > 1){
        centroids<-Reduce("rbind", lapply(unique(clusters), function(x){
          cells<-names(clusters[which(clusters %in% x)])
          centroid<-colMeans(pca[cells,])
          return(centroid)
        }))
        rownames(centroids) <- unique(clusters)
      }
    }
  }
  return(clusters)
}

renameClusters <- function(clusters){
  name <- 0
  for (i in names(table(clusters)[order(table(clusters), decreasing = TRUE)])){
    clusters[which(clusters == i)] <- as.character(name)
    name <- name + 1
  }
  return(clusters)
}

### Compute DEscore for all clusters
compute_DEscore <-function(obj, clusters){
  obj$clusters <- clusters[rownames(obj@meta.data)]
  Idents(obj) <- 'clusters'
  DEscore <-c()
  for (cluster in unique(clusters)){
    #markers<-FindMarkers(obj, ident.1 = cluster, only.pos = TRUE, assay = 'RNA',logfc.threshold = log(1.5), pseudocounts.use=0.1, min.pct = 0.3,min.diff.pct = 0.1)
    #markers<-FindMarkers(obj, ident.1 = cluster, only.pos = TRUE, assay = 'RNA',logfc.threshold = log(1.5), verbose = FALSE)
    ribo_genes <- grep(pattern = "^Rp[sl]", x = rownames(obj@assays$RNA@counts), value = TRUE)
    markers<-FindMarkers(obj, ident.1 = cluster, only.pos = TRUE, assay = 'RNA',logfc.threshold = log2(2), verbose = FALSE)
    markers<-markers[which(markers$p_val_adj<=0.01),]
    markers<-markers[which(rownames(markers) %!in% ribo_genes),]
    if(nrow(markers)==0){
      DEscore <- c(DEscore,0)
    }else{
      tmp<- log10(markers$p_val_adj)*(-1)
      tmp[which(tmp>20)]<-20
      DEscore<-c(DEscore,sum(tmp))
    }
  }
  names(DEscore) <- unique(clusters)
  return(DEscore)
} 

### check if all clusters in provided resolution pass DEscore.cutoff and have at least minSize cells
### If a cluster do not pass cutoffs it is merged to its closest cluster (with mergeClosestCluster function).
### Testing is repeated on the resulting clusters until all clusters pass cutoffs or we end up with only one cluster.
testStableClusters <- function(obj, resolution, nPCs, minSize=149, DEscore.cutoff=60){
  clusters <- obj@meta.data[,resolution]
  clusters <- as.vector(clusters)
  names(clusters) <- rownames(obj@meta.data)
  n_clusters <- length(unique(clusters))
  if (n_clusters == 1){
    clusters<-renameClusters(clusters)
    return(clusters)
  }else{
    # check cluster dimension and merge clusters with few cells with closest cluster (and repeat)
    PassClusters<-table(clusters) > minSize
    while (sum(PassClusters) < n_clusters){
      clusters <- mergeClosestCluster(obj, clusters, PassClusters, nPCs)
      PassClusters<-table(clusters) > minSize
      n_clusters <- length(unique(clusters))
    }
    # compute DEscore for resulting clusters and merge clusters not passing DEscore cutoff with closest cluster (and repeat)
    if (n_clusters == 1){
      clusters<-renameClusters(clusters)
      return(clusters)
    }else{
      DEscore<-compute_DEscore(obj, clusters)
      PassClusters<-DEscore > DEscore.cutoff
      while (sum(PassClusters) < n_clusters & n_clusters > 1){
        clusters <- mergeClosestCluster(obj, clusters, PassClusters, nPCs)
        n_clusters <- length(unique(clusters))
        if(n_clusters > 1){
          DEscore<-compute_DEscore(obj, clusters)
          PassClusters<-DEscore > DEscore.cutoff
        }
      }
    }
    clusters<-renameClusters(clusters)
    return(clusters)
  }
}

Iteration<-0
### call iteratively all processes
callIteration <- function(obj, batch, nPC=NULL, min.res=0.1, max.res=1, perc.sub=0.8, n_subsampling=20, jaccard_cutoff=0.8, percent_cutoff=0.8, minSize=149, DEscore.cutoff=60,cr.dir){
  parent_dir<-getwd()
  obj <- computeIntegrationAndClusters(obj, batch, nPCs=nPC, min.res=min.res, max.res=max.res)
  if(is.null(nPC)){
    nPCs <- readRDS("optimal_nPCs_integrated.RDS")
  }else{
    nPCs <- nPC
  }
  resolutions <- findResolutions(obj, min.res=min.res, max.res=max.res)
  subsample_idents <- performSubsamplingAndClustering(obj, resolutions, nPCs, perc.sub=perc.sub, n_subsampling=n_subsampling)
  saveRDS(subsample_idents,"subsample_idents.RDS")
  fullsample_idents <- prepareFullData(obj, resolutions)
  saveRDS(fullsample_idents,"fullsample_idents.RDS")
  stable_resolution <- findStableClustering(subsample_idents, fullsample_idents, jaccard_cutoff=jaccard_cutoff, percent_cutoff=percent_cutoff)
  if(!(is.null(stable_resolution))){
    resolution <- paste0("integrated_snn_res.",stable_resolution)
    write.table(resolution, "Selected_stable_resolution.txt")
    clusters <- testStableClusters(obj, resolution,nPCs, minSize=minSize, DEscore.cutoff=DEscore.cutoff)
    saveRDS(clusters, "clusters.RDS")
    n_clusters <- length(unique(clusters))
    if(n_clusters > 1){
      final.clusters[names(clusters)] <<- paste(final.clusters[names(clusters)],clusters,sep=".")
      if (length(names(table(clusters)[table(clusters) >= 300])) != 0){
        lapply(names(table(clusters)[table(clusters) >= 300]),function(cl){
          Iteration <<- Iteration + 1
          print(paste0("Iteration ",Iteration,"- Cluster ",cl))
          obj <- subset(obj, cells=names(clusters[clusters %in% cl]))
          setwd(cr.dir)
          dir.create(paste0("Cluster_",cl))
          setwd(paste0("Cluster_",cl))
          cr.dir<-paste0(cr.dir,"/Cluster_",cl)
          callIteration(obj, batch=batch, nPC=nPC, min.res=min.res, max.res=max.res, perc.sub=perc.sub, n_subsampling=n_subsampling, jaccard_cutoff=jaccard_cutoff, percent_cutoff=percent_cutoff, minSize=minSize, DEscore.cutoff=DEscore.cutoff,cr.dir=cr.dir)
          gc()
        })
      }
    }
  }else{
    write.table(print("no stable resolution found"), "No_stable_resolution.txt", row.names = FALSE, col.names = FALSE)
  }
  print("Iteration ended")
}

