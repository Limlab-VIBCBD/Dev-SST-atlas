process compute_nPCs_and_clusterings_iteration {

    publishDir "${params.tmp_path}", mode: 'copy', pattern: "object_with_clustering*"
    publishDir "${params.tmp_path}", mode: 'copy', pattern: "optimal_nPCs_*"

  input:
    tuple val(cl), path(error_files),path(object), path(data_use), path(nPC_5)
    val min_resolution
    val max_resolution



  output:
    tuple val(cl), path("object_with_clustering_${cl}.RDS"), path("optimal_nPCs_${cl}.RDS")


  script:
    """
    #!/usr/bin/env Rscript
    options(warn=-1)
    library(Seurat)
    '%!in%' <- function(x,y)!('%in%'(x,y))
    min.res=as.numeric("${min_resolution}")
    max.res=as.numeric("${max_resolution}")
    obj<-readRDS("${object}")
    nPC_5<-readRDS("${nPC_5}")
    PC <-seq(from = nPC_5-5, to = nPC_5+5, by = 1)
    files<-"${error_files}"
    eval(parse(text=paste0("files<-c(\\"",gsub(" ","\\",\\"",files),"\\")")))
    error<-Reduce( "rbind", lapply(files, function(x){
        error<-readRDS(x)
        return(error)
    }))
    errors<-colSums(error)
    nPC=PC[which(errors == min(errors))]
    saveRDS(nPC,'optimal_nPCs_${cl}.RDS')
    #obj<-RunUMAP(obj,dims=1:nPC, verbose = FALSE)
    obj<-FindNeighbors(obj,dims=1:nPC, verbose = FALSE)
    resolutions_selected <- seq(min.res, max.res, 0.1)
    obj<-FindClusters(obj,resolution = resolutions_selected, verbose = FALSE)
    if ('integrated' %!in% Assays(obj)){
      obj@meta.data <- obj@meta.data[,colnames(obj@meta.data) %!in% colnames(obj@meta.data)[grep("integrated_snn_res.",colnames(obj@meta.data))]]
      colnames(obj@meta.data)[grep("RNA_snn_res.",colnames(obj@meta.data))] <- gsub("RNA_snn_res.","integrated_snn_res.",colnames(obj@meta.data)[grep("RNA_snn_res.",colnames(obj@meta.data))])
    }  
    saveRDS(obj, "object_with_clustering_${cl}.RDS")
    q('no')
    """
}