
process split_object {

    input:
      tuple path(object), path(clusters), val(cl), val(cl_updated) 
      val(min_size)

    output:
        tuple path("object_${cl_updated}.RDS"), val(cl_updated), path("${cl_updated}_size.txt")

    script:
    """
    #!/usr/bin/env Rscript
    library(Seurat)
    obj<-readRDS("${object}")
    clusters<-readRDS("${clusters}")
    cl=${cl}-1
    min_size<-as.numeric(${min_size})
    sel_cells <- names(clusters[which(clusters == cl)])
    obj_sub<-subset(obj, cells= sel_cells)
    saveRDS(obj_sub,"object_${cl_updated}.RDS")
    write.table(length(sel_cells),"${cl_updated}_size.txt",sep="\t",quote=FALSE,col.names=FALSE,row.names=FALSE)
 """
}