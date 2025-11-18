
process integration {

    input:
      path object
      val tmppath
      val reproducibility

    output:
      path "obj_with_pca.RDS", emit: obj_with_pca
      path "data_use.RDS", emit: data_use
      path "dgem.kfold.RDS", emit: dgem_kfold

    script:
      """
      #!/usr/bin/env Rscript
      library(Seurat)
      library(missMDA)
      options(future.globals.maxSize=1000000000)
      '%!in%' <- function(x,y)!('%in%'(x,y))
      obj <- readRDS("${object}")
      if ('batch' %!in% colnames(obj@meta.data)) stop("Error: seurat object does not contain batch column in meta.data")
      if ('RNA' %!in% Assays(obj)) stop("Error: seurat object does not contain RNA assay")
      DefaultAssay(obj) <- 'RNA'
      remove_samples<-names(table(obj@meta.data[,"batch"]))[table(obj@meta.data[,"batch"]) < 20]
      if(length(remove_samples) > 0){
        remove_cells <- colnames(subset(obj, subset=batch %in% remove_samples))
        obj<-subset(obj, subset=batch %!in% remove_samples)
      }
      if(length(table(obj\$batch)) > 1){
        integrated_all_list <- SplitObject(obj, split.by = "batch")
        integrated_all_list <- lapply(integrated_all_list, function(sample) {
          sample<-NormalizeData(sample, verbose = FALSE)
          sample<-FindVariableFeatures(sample,nfeatures = 2000, verbose = FALSE)
          return(sample)
        })
        i.anchors <- FindIntegrationAnchors(object.list = integrated_all_list, dims = 1:min(min(table(obj@meta.data[,"batch"]))-5,30),reduction='cca',scale=T,k.anchor=5,k.filter=100,k.score=15,anchor.features=3000,verbose = FALSE)
        obj <- IntegrateData(anchorset = i.anchors, dims = 1:min(min(table(obj@meta.data[,"batch"]))-5,30), normalization.method ='LogNormalize', k.weight=min(100,min(table(obj@meta.data[,"batch"]))),verbose = FALSE)
        DefaultAssay(obj)<-'integrated'
        obj<-ScaleData(obj,vars.to.regress = c("nFeature_RNA",'percent.mt','ccDiff')[which(c("nFeature_RNA",'percent.mt','ccDiff') %in% colnames(obj@meta.data))], verbose = FALSE)
        data.use <- GetAssayData(obj, assay = "integrated",layer = "scale.data")
      }else{
        obj<-NormalizeData(obj, verbose = FALSE)
        obj<-FindVariableFeatures(obj,nfeatures = 2000, verbose = FALSE)
        obj<-ScaleData(obj,vars.to.regress = c("nFeature_RNA",'percent.mt','ccDiff')[which(c("nFeature_RNA",'percent.mt','ccDiff') %in% colnames(obj@meta.data))], verbose = FALSE)
        data.use <- GetAssayData(obj, assay = "RNA",layer = "scale.data")
      }
      obj<-RunPCA(obj,npcs = 60, verbose = FALSE)
      saveRDS(obj,'obj_with_pca.RDS')
      genes.use = VariableFeatures(object = obj)
      genes.use <- unique(x = genes.use[genes.use %in% rownames(x = data.use)])
      genes.var <- apply(X = data.use[genes.use, ], MARGIN = 1, FUN = var)
      genes.use <- genes.use[genes.var > 0]
      genes.use <- genes.use[! is.na(x = genes.use)]
      data.use <- data.use[genes.use, ]
      saveRDS(data.use,'data_use.RDS')
      if (as.character("${reproducibility}")){
        set.seed(1234)
        dgem.kfold<-dismo::kfold(t(data.use), k=10)
      }else{
        dgem.kfold<-dismo::kfold(t(data.use), k=10)
      }  
      saveRDS(dgem.kfold,'dgem.kfold.RDS')
      q('no')
      """
  }