
process identify_stable_clusters_iteration {
    
    publishDir "${params.tmp_path}", mode: 'copy', pattern: "Selected_stable_resolution*"
    publishDir "${params.tmp_path}", mode: 'copy', pattern: "*_idents_*"


    input:
    tuple val(cl), path(recluster_ident_files), path(object), path(nPCs)
    val jaccard_cutoff
    val percent_cutoff

    output:
    tuple path(object), path(nPCs), val(cl), path("Selected_stable_resolution_${cl}.txt"), path("subsample_idents_${cl}.RDS"), path("fullsample_idents_${cl}.RDS")

        script:
    """
    #!/usr/bin/env Rscript
    library(dplyr)
    library(tidyr)
    library(purrr)
    library(scclusteval)
    files<-"${recluster_ident_files}"
    eval(parse(text=paste0("files<-c(\\"",gsub(" ","\\",\\"",files),"\\")")))
    obj<-readRDS("${object}")
    cl<-as.character("${cl}")
    jaccard_cutoff=as.numeric("${jaccard_cutoff}")
    percent_cutoff=as.numeric("${percent_cutoff}")
    res <- unique(unlist(lapply(files, function(x) strsplit(x,"_")[[1]][[3]])))
    n_subsampling <- max(as.numeric(unique(unlist(lapply(files, function(x) gsub(".RDS","",strsplit(x,"_")[[1]][[4]]))))))
    subsample_idents<-Reduce( "bind_rows", lapply(res, function(x){
        recluster_ident <-list()
        original_ident <-list()
        for(i in 1:n_subsampling){
            out<-readRDS(paste0("recluster_ident_",x,'_',i,"_",cl,".RDS"))
            recluster_ident[[i]]=out[[1]]
            original_ident[[i]]=out[[2]]
        }
        data<-tibble(resolution=x, recluster_ident=recluster_ident, original_ident=original_ident,round=as.character(c(0:(n_subsampling-1))))
    }))
    saveRDS(subsample_idents,'subsample_idents_${cl}.RDS')
    resolutions<-paste("integrated_snn_res.",res,sep="")
    original_ident_full<-lapply(resolutions, function(x){
      original_ident_full <- obj@meta.data[,x]
      names(original_ident_full) <- rownames(obj@meta.data)
      return(original_ident_full)
    })
    fullsample_idents <- tibble(resolution=res,original_ident_full=original_ident_full)
    saveRDS(fullsample_idents,'fullsample_idents_${cl}.RDS')
    subsample_idents_list<- subsample_idents %>% group_by(resolution) %>%  nest()
    stable_clusters<- subsample_idents_list %>% mutate(stable_cluster = map(data, ~ AssignStableCluster(.x[["original_ident"]],.x[["recluster_ident"]],jaccard_cutoff = jaccard_cutoff,method = "jaccard_percent", percent_cutoff = percent_cutoff)))
    df <- dplyr::left_join(stable_clusters, fullsample_idents) %>% dplyr::ungroup() %>% dplyr::mutate(total = map_dbl(stable_cluster, ~length(.x[["stable_cluster"]]))) %>%
      dplyr::mutate(stable = map_dbl(stable_cluster, ~.x[["number_of_stable_cluster"]])) %>%
      dplyr::mutate(percentage_cluster = map2_dbl(stable, total, function(x, y) x/y)) %>%
      dplyr::mutate(percentage = map2_dbl(original_ident_full, stable_cluster, function(x, y) CalculatePercentCellInStable(x,y[["stable_cluster"]]))) %>%
      dplyr::select(-data, -stable_cluster, -original_ident_full) %>% dplyr::mutate_if(is.character, function(x) as.factor(as.numeric(x))) %>%
      tidyr::gather(total:stable,key = "category", value = "number")
    # check stability and choose best clustering
    df<-subset(df, subset=percentage > 0.7 & percentage_cluster >= 0.6 & category== "total")
    if (nrow(df) > 0){
      df<-subset(df, subset=percentage == max(df[[3]]))
      df<-subset(df, subset=number == max(df[5]))
      stable_resolution <- df[[1]][1]
    }else{
      stable_resolution <- NULL
    }
    if(!(is.null(stable_resolution))){
      resolution <- paste0("integrated_snn_res.",stable_resolution)
      write.table(resolution, "Selected_stable_resolution_${cl}.txt", row.names = FALSE, col.names = FALSE)
    }else{
      write.table(print("integrated_snn_res.0.0"), "Selected_stable_resolution_${cl}.txt", row.names = FALSE, col.names = FALSE)
    }
 """
}