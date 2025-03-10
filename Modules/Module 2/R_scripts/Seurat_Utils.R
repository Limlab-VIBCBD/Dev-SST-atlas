#Initialize Functions
PrepDR <- function( # From Seurat
  object,
  genes.use = NULL,
  use.imputed = FALSE,
  assay.type="RNA"
) {
  
  if (length(VariableFeatures(object = object)) == 0 && is.null(x = genes.use)) {
    stop("Variable genes haven't been set. Run MeanVarPlot() or provide a vector
         of genes names in genes.use and retry.")
  }
  if (use.imputed) {
    data.use <- t(x = scale(x = t(x = object@imputed)))
  } else {
    data.use <- GetAssayData(object, assay = assay.type,slot = "scale.data")
  }
  genes.use <- if(is.null(genes.use)) VariableFeatures(object = object) else genes.use # Changed
  genes.use <- unique(x = genes.use[genes.use %in% rownames(x = data.use)])
  genes.var <- apply(X = data.use[genes.use, ], MARGIN = 1, FUN = var)
  genes.use <- genes.use[genes.var > 0]
  genes.use <- genes.use[! is.na(x = genes.use)]
  data.use <- data.use[genes.use, ]
  return(data.use)
  }

PCA_estimate_nPC<-function(data, whereto, k=10, from.nPC = 2, to.nPC=150, by.nPC=5, maxit=200, seed=617) {
  library(missMDA)
  PC <-seq(from = from.nPC, to = to.nPC, by = by.nPC)
  # Init the error matrices
  error1<-matrix(0, nrow = length(c(1:k)), ncol = length(PC))
  error2<-matrix(0, nrow = length(c(1:k)), ncol = length(PC))
  print(paste0(k,"-fold paritioning..."))
  # K-Fold Partitioning
  dgem.kfold<-dismo::kfold(t(data), k=k)
  # SVD-CV based on https://stats.stackexchange.com/questions/93845/how-to-perform-leave-one-out-cross-validation-for-pca-to-determine-the-number-of
  for(i in c(1:k)) {
    print(paste0("k:",i))
    X.train<-t(data[, dgem.kfold!=i])
    X.test<-t(data[, dgem.kfold==i])
    # Find a few approximate singular values and corresponding singular vectors of a matrix.
    print("Running SVD...")
    # Seurat uses IRLBA to do PCA : https://github.com/satijalab/seurat/blob/cec7cb95c73fd6d605723e9af9a1f96eda5635de/R/dimensional_reduction.R
    pca.results<-irlba::irlba(A = X.train, nv = to.nPC, maxit = maxit) # Otherwise, default maxit=100 do not converge
    gl<-pca.results$v
    for(j in 1:length(PC)) {
      print(paste0("Ndims:",PC[j]))
      P<-gl[,1:PC[j]]%*%t(gl[,1:PC[j]])
      # Naive method
      err1<-X.test %*% (diag(dim(P)[1]) - P)
      # Approximate method
      err2<-X.test %*% (diag(dim(P)[1]) - P + diag(diag(P)))
      error1[i,j]<-sum(err1^2)
      error2[i,j]<-sum(err2^2)
      rm(err1)
      rm(err2)
    }
  }
  errors1<-colSums(error1)
  errors2<-colSums(error2)
  nPC=PC[which(errors2 == min(errors2))]
  saveRDS(nPC,whereto)
  plot(PC,errors1)
  plot(PC,errors2)
  return(nPC)
}