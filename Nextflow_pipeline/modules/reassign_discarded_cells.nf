
process reassign_discarded_cells {

    input:
      val(n)
      val outpath
      val tmppath


    script:
    """
    #!/usr/bin/env Rscript
    library(Seurat)
    library(dbscan)
    library(dplyr)
    obj<-readRDS("${tmppath}/object_with_clustering.RDS")
    final_clusters<-rep("c", dim(obj)[2])
    names(final_clusters) <- rownames(obj@meta.data)
    sel_cells<-readLines("${tmppath}/remove_cells.txt")
    sel_cells<-unique(sel_cells)
    nPCs<-readRDS("${tmppath}/optimal_nPCs.RDS")
    obj<-RunUMAP(obj,dims=1:nPCs, verbose = FALSE)
    setwd("${tmppath}")
    files<-dir()
    files<-files[grep("clusters",files)]
    files<-files[grep(".RDS",files)]
    max_it<-max(nchar(files))
    tmp<-files[nchar(files) == 12]
    clusters<-Reduce( "c", lapply(tmp, function(x){
      cluster<-readRDS(x)
      return(cluster)
    }))
    final_clusters[names(clusters)] <-paste(final_clusters[names(clusters)],clusters,sep=".")
    i=14
    while (i <= (max_it-1)){
      tmp<-files[nchar(files) == i]
      clusters<-Reduce( "c", lapply(tmp, function(x){
        cluster<-readRDS(x)
        return(cluster)
      }))
      final_clusters[names(clusters)] <-paste(final_clusters[names(clusters)],clusters,sep=".")
      i=i+1
    }
    pca_integrated <- Embeddings(obj, reduction = "pca")[,1:nPCs]
    knn_integrated <- dbscan::kNN(x = pca_integrated %>% as.matrix(), k = 10)
    knn_integrated.data <- data.frame(from = rep(rownames(knn_integrated\$id), 10), to = rownames(pca_integrated)[as.vector(knn_integrated\$id)])
    knn_integrated.data <- knn_integrated.data[which(knn_integrated.data\$from %in% sel_cells),]
    assignment <-unlist(lapply(sel_cells, function(x){
      tmp <- knn_integrated.data[knn_integrated.data\$from == x,"to"]
      tmp <- tmp[!(tmp %in% sel_cells)]
      cluster <- names(table(final_clusters[tmp])[table(final_clusters[tmp]) == max(table(final_clusters[tmp]))])
      if(length(cluster) > 1){
        cluster_dist <- unlist(lapply(cluster, function(cl){
          cl_names<-names(final_clusters[tmp][final_clusters[tmp] %in% cl])
          dist_value<-median(dist(pca_integrated[c(x,cl_names),], method = "euclidean",upper = TRUE))
          return(dist_value)
        }))
        names(cluster_dist) <- cluster
        cluster <- names(cluster_dist)[which(cluster_dist == min(cluster_dist))]
      }
    return(cluster)
    }))
    final_clusters[sel_cells] <- assignment
    obj@meta.data<-obj@meta.data[,!(1:ncol(obj@meta.data) %in% grep("integrated_snn",colnames(obj@meta.data)))]
    obj\$final_clusters <- final_clusters[colnames(obj)]
    saveRDS(final_clusters,"${outpath}/final_clusters.RDS")
    saveRDS(obj, "${outpath}/object_with_final_clusters.RDS")
    """
}