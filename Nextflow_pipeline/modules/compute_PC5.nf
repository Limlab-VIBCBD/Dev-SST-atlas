
process compute_PC_5 {

  input:
    path data_use
    path error_files
    val reproducibility

  output:
    path "n_PCs_5.RDS" , emit: nPCs_5
    path "dgem.kfold.RDS", emit: dgem_kfold

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
    saveRDS(nPC,'n_PCs_5.RDS')
    data<-readRDS("${data_use}")
    if (as.character("${reproducibility}")){
      set.seed(2345)
      dgem.kfold<-dismo::kfold(t(data), k=10)
    }else{
      dgem.kfold<-dismo::kfold(t(data), k=10)
    }
    saveRDS(dgem.kfold,'dgem.kfold.RDS')
    """
}