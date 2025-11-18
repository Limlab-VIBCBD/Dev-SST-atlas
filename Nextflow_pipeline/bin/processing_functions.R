library(Seurat)
library(dplyr)

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
  cell_names <- names(clusters)
  clusters <- paste("cl",clusters, sep="_")
  names(clusters) <- cell_names
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

