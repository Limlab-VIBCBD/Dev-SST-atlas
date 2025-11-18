
process prepare_PC_1 {

  input:
    path data_use
    path nPC_5
    path dgem_kfold
    val k

  output:
    path "error_${k}.RDS"

  script:
    """
    #!/usr/bin/env Rscript
    data<-readRDS("${data_use}")
    dgem.kfold<-readRDS("${dgem_kfold}")
    nPC_5<-readRDS("${nPC_5}")
    PC <-seq(from = nPC_5-5, to = nPC_5+5, by = 1)
    error<-c()
    k<-as.numeric("${k}")
    X.train<-t(data[, dgem.kfold!=k])
    X.test<-t(data[, dgem.kfold==k])
    pca.results<-irlba::irlba(A = X.train, nv = nPC_5+5, maxit = 200)
    gl<-pca.results[['v']]
    for(j in 1:length(PC)) {
        P<-gl[,1:PC[j]]%*%t(gl[,1:PC[j]])
        # Approximate method
        err<-X.test %*% (diag(dim(P)[1]) - P + diag(diag(P)))
        error<-c(error,sum(err^2))
    }
    saveRDS(error,'error_${k}.RDS')
    """
}