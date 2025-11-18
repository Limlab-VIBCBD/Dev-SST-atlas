
process update_cluster_id {

    input:
      tuple path(object), path(clusters), val(cl)

    output:
      tuple path(object), path(clusters), val(cl), path("cl_updated.txt")

    script:
    """
    #!/usr/bin/env Rscript
    cl<-as.character("${cl}")
    p_cl <- as.character("${object}")
    p_cl<-sub(".*object_with_clustering",'',p_cl)
    p_cl<-gsub("\\\\.RDS",'',p_cl)
    p_cl<-gsub("_",'',p_cl)
    cl_updated<-paste0(p_cl,cl)
    write.table(cl_updated,"cl_updated.txt",sep="\t", quote=FALSE, col.names=FALSE, row.names=FALSE)
    """
}