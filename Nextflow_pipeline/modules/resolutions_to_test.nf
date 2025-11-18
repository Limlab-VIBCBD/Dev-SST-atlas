
process resolutions_to_test {

  input:
    path object
    val min_resolution
    val max_resolution
    val tmppath

  output:
    path "resolutions.csv" , emit: resolutions



  script:
    """
    #!/usr/bin/env Rscript
    library(Seurat)
    min.res=as.numeric("${min_resolution}")
    max.res=as.numeric("${max_resolution}")
    obj<-readRDS("${object}")
    resolutions_selected <- seq(min.res, max.res, 0.1)
    resolutions <- paste0("integrated_snn_res.",resolutions_selected)
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
  resolution_for_stability <- gsub("integrated_snn_res.","",resolution_for_stability)
  write.table(resolution_for_stability, "resolutions.csv",quote=FALSE, col.names=FALSE, row.names=FALSE,sep=",")
 """
}