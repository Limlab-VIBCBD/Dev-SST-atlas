
process subsampling_and_clustering {

    input:
      path object
      path nPCs
      tuple val(resolution), val(n_subsampling)
      val perc_sub
      val reproducibility

    output:
      path "recluster_ident_${resolution}_${n_subsampling}.RDS"

    script:
    """
    #!/usr/bin/env Rscript
      library(Seurat)
      library(scclusteval)
      obj<-readRDS("${object}")
      nPCs <- readRDS("${nPCs}")
      res <- "${resolution}"
      n_subsampling <- ${n_subsampling}
      x<-paste0("integrated_snn_res.",res)
      if (as.character("${reproducibility}")){
        rand_test <- RandomSubsetData(obj, rate=${perc_sub},random.subset.seed=as.numeric(res)*100+as.numeric(n_subsampling))
      }else{
        rand_test <- RandomSubsetData(obj, rate=${perc_sub})
      }
      rand_test <- ScaleData(rand_test,vars.to.regress = c("nFeature_RNA",'percent.mt','ccDiff')[which(c("nFeature_RNA",'percent.mt','ccDiff') %in% colnames(obj@meta.data))], verbose = FALSE)
      rand_test <- RunPCA(rand_test, verbose = FALSE,npcs = 60)
      rand_test <- FindNeighbors(rand_test, reduction = "pca", dims = 1:nPCs, verbose = FALSE)
      eval(parse(text=paste0("rand_test <- FindClusters(rand_test, resolution = ",res,", verbose = FALSE)")))
      recluster_ident<-rand_test@meta.data[,x]
      names(recluster_ident) <- rownames(rand_test@meta.data)
      original_ident<-obj@meta.data[rownames(rand_test@meta.data),x]
      names(original_ident) <- rownames(rand_test@meta.data)
      saveRDS(list(recluster_ident,original_ident),"recluster_ident_${resolution}_${n_subsampling}.RDS")
    """
  }