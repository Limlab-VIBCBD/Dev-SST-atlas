
process compute_nPCs_and_clusterings {

  input:
    path object
    path nPC_5
    path error_files
    val min_resolution
    val max_resolution


  output:
    path "object_with_clustering.RDS" , emit: obj_with_clustering
    path "optimal_nPCs.RDS" , emit: nPCs


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
    saveRDS(nPC,'optimal_nPCs.RDS')
    #obj<-RunUMAP(obj,dims=1:nPC, verbose = FALSE)
    obj<-FindNeighbors(obj,dims=1:nPC, verbose = FALSE)
    resolutions_selected <- seq(min.res, max.res, 0.1)
    obj<-FindClusters(obj,resolution = resolutions_selected, verbose = FALSE)
    if ('integrated' %!in% Assays(obj)){
      obj@meta.data <- obj@meta.data[,colnames(obj@meta.data) %!in% colnames(obj@meta.data)[grep("integrated_snn_res.",colnames(obj@meta.data))]]
      colnames(obj@meta.data)[grep("RNA_snn_res.",colnames(obj@meta.data))] <- gsub("RNA_snn_res.","integrated_snn_res.",colnames(obj@meta.data)[grep("RNA_snn_res.",colnames(obj@meta.data))])
    }  
    saveRDS(obj, "object_with_clustering.RDS")
    """
}