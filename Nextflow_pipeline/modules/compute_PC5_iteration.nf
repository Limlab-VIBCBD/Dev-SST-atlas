process compute_PC_5_iteration {

  input:
    tuple val(cl), path(error_files), path(data_use)
    val reproducibility

  output:
    tuple path(data_use), path("n_PCs_5_${cl}.RDS"), path("dgem.kfold_${cl}.RDS"), val(cl)

  script:
    """
    #!/usr/bin/env Rscript
    library(missMDA)
    library(dplyr)
    files<-"${error_files}"
    eval(parse(text=paste0("files<-c(\\"",gsub(" ","\\",\\"",files),"\\")")))
    error<-Reduce( "rbind", lapply(files, function(x){
        error<-readRDS(x)
        return(error)
    }))
    errors<-colSums(error)
    PC <-seq(from = 10, to = 55, by = 5)
    nPC=PC[which(errors == min(errors))]
    saveRDS(nPC,'n_PCs_5_${cl}.RDS')
    data<-readRDS("${data_use}")
    if (as.character("${reproducibility}")){
      set.seed(2345)
      dgem.kfold<-dismo::kfold(t(data), k=10)
    }else{
      dgem.kfold<-dismo::kfold(t(data), k=10)
    }
    saveRDS(dgem.kfold,'dgem.kfold_${cl}.RDS')
    """
}