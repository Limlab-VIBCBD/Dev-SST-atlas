
process compute_final_clusters {

    input:
      path object
      path nPCs
      path resolution
      val minSize
      val DEscore_cutoff


    output:
      path "clusters.RDS", emit: clusters
      path "n_clusters.txt", emit: num_clusters

    script:
    """
    #!/usr/bin/env Rscript
    '%!in%' <- function(x,y)!('%in%'(x,y))
    source("${projectDir}/bin/processing_functions.R")
    obj<-readRDS("${object}")
    minSize=as.numeric("${minSize}")
    nPCs <- readRDS('${nPCs}')
    DEscore.cutoff=as.numeric("${DEscore_cutoff}")
    resolution <- read.table("${resolution}")
    resolution<-resolution[,1]
    if(resolution == "integrated_snn_res.0.0"){
        clusters <-rep("0", nrow(obj@meta.data))
        names(clusters) <- rownames(obj@meta.data)
    }else{
        clusters <- obj@meta.data[,resolution]
        clusters <- as.vector(clusters)
        names(clusters) <- rownames(obj@meta.data)
        n_clusters <- length(unique(clusters))
        if (n_clusters == 1){
            clusters<-renameClusters(clusters)
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
        }
    }
    n_clusters <- length(unique(clusters))
    saveRDS(clusters,'clusters.RDS')
    write.table(n_clusters,"n_clusters.txt",sep="\t", quote=FALSE, col.names=FALSE, row.names=FALSE)
 """
}
