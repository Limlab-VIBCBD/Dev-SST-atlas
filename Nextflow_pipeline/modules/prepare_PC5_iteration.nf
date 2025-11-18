process prepare_PC_5_iteration {

  input:
    tuple path(data_use), path(dgem_kfold),val(cl), val(k)

  output:
    tuple path("error_${cl}_${k}.RDS"), val(cl)

  script:
    """
    #!/usr/bin/env Rscript
    data<-readRDS("${data_use}")
    dgem.kfold<-readRDS("${dgem_kfold}")
    PC <-seq(from = 10, to = 55, by = 5)
    error<-c()
    k<-as.numeric("${k}")
    X.train<-t(data[, dgem.kfold!=k])
    X.test<-t(data[, dgem.kfold==k])
    pca.results<-irlba::irlba(A = X.train, nv = 55, maxit = 200)
    gl<-pca.results[['v']]
    for(j in 1:length(PC)) {
        P<-gl[,1:PC[j]]%*%t(gl[,1:PC[j]])
        # Approximate method
        err<-X.test %*% (diag(dim(P)[1]) - P + diag(diag(P)))
        error<-c(error,sum(err^2))
    }
   saveRDS(error,'error_${cl}_${k}.RDS')
   q('no')
   """
}