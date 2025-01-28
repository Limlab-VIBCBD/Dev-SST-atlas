library(stringr)
library(Seurat)
library(dittoSeq)
library(ggplot2)
library(gtools)
library(dplyr)
library(dittoSeq)
library(R.utils)
library(colorspace)
library(scclusteval)
library(tidyr)
library(purrr)
library(parallel)
library(MetaNeighbor)
'%!in%' <- function(x,y)!('%in%'(x,y))

atlas_v2<-readRDS('Integrated_atlas_filtered_annotated.RDS')

### compute subsamplings
cell_subsampling <-list()
for (perc.sub in seq(0.05, 0.95, 0.05)){
  cells_sub_list<-lapply(1:10, function(x){
    set.seed(NULL)
    cells <- sample(colnames(atlas_v2), size=round(dim(atlas_v2)[2]*perc.sub))
    return(cells)
  })
  cell_subsampling[[as.character(perc.sub)]] <- cells_sub_list
}
#saveRDS(cell_subsampling,"CellsSubsampling.RDS")


for (perc.sub in seq(0.05, 0.95, 0.05)){
  n <- lapply(cell_subsampling[[as.character(perc.sub)]], function(x){
    obj<-subset(atlas_v2, cells=x)
    table(obj$Stage)
  })
}


### Number of DEGs between clusters
atlas_v2_LRP<-readRDS("obj_with_finalClusters.RDS")
Idents(atlas_v2_LRP) <- 'final_clusters'
DefaultAssay(atlas_v2_LRP) <-'RNA'
majorClassMarkers<-FindAllMarkers(atlas_v2_LRP,logfc.threshold = log2(1.5),min.pct = 0.05,only.pos = TRUE, test.use = 'MAST', latent.vars = 'sample')
majorClassMarkers<-majorClassMarkers[which(majorClassMarkers$p_val_adj<=0.01),]
majorClassMarkers$subset <-"LRP"
atlas_v2_NM<-readRDS("obj_with_finalClusters.RDS")
Idents(atlas_v2_NM) <- 'final_clusters'
DefaultAssay(atlas_v2_NM) <-'RNA'
markers<-FindAllMarkers(atlas_v2_NM,logfc.threshold = log2(1.5),min.pct = 0.05,only.pos = TRUE, test.use = 'MAST', latent.vars = 'sample')
markers<-markers[which(markers$p_val_adj<=0.01),]
markers$subset <- "Non-Martinotti"
majorClassMarkers<-rbind(majorClassMarkers,markers)
atlas_v2_M<-readRDS("obj_with_finalClusters.RDS")
Idents(atlas_v2_M) <- 'final_clusters'
DefaultAssay(atlas_v2_M) <-'RNA'
markers<-FindAllMarkers(atlas_v2_M,logfc.threshold = log2(1.5),min.pct = 0.05,only.pos = TRUE, test.use = 'MAST', latent.vars = 'sample')
markers<-markers[which(markers$p_val_adj<=0.01),]
markers$subset <- "Martinotti"
majorClassMarkers<-rbind(majorClassMarkers,markers)
n_genes_all <- length(majorClassMarkers$gene)
print(n_genes_all) #710
saveRDS(majorClassMarkers,"allDEGs_final_clusters_subsets_Atlas_v2.RDS")

DEGs_list <-list()
for (perc.sub in seq(0.05, 0.95, 0.05)){
  print(perc.sub)
  markers_sub <- lapply(cell_subsampling[[as.character(perc.sub)]], function(x){
    markers_sub<-Reduce("rbind",lapply(list(atlas_v2_LRP,atlas_v2_NM,atlas_v2_M), function(obj){
      atlas_v2_sub <-subset(obj, cells=x[x %in% colnames(obj)])
      Idents(atlas_v2_sub) <- 'final_clusters'
      DefaultAssay(atlas_v2_sub) <-'RNA'
      tmp<-FindAllMarkers(atlas_v2_sub, logfc.threshold = log2(1.5),min.pct = 0.05,only.pos = TRUE, test.use = 'MAST', latent.vars = 'sample',verbose = FALSE)
      tmp<-tmp[which(tmp$p_val_adj<=0.01),]
    }))
    return(markers_sub)
  })
  DEGs_list[[as.character(perc.sub)]] <- markers_sub
}
saveRDS(DEGs_list,"allDEGs_final_clusters_subsets.RDS")


###Plot
DEGs_number<-list()
for (perc.sub in seq(0.05, 0.95, 0.05)){
  DEGs_number[[as.character(perc.sub)]]<-unlist(lapply(DEGs_list[[as.character(perc.sub)]], function(x) {
    x <- subset(x,pct.1 > 0.2 & avg_log2FC>log2(1.8))
    dim(x)[1]
  }))
}
DEGs_number_all<-dim(subset(majorClassMarkers,pct.1 > 0.2 & avg_log2FC>log2(1.8)))[1]

avg<-rev(unlist(lapply(DEGs_number, mean))/DEGs_number_all)
x<-round(seq(0.05, 0.95, 0.05)*51205)
err <- unlist(lapply(DEGs_number, function(x) sd(x/DEGs_number_all)/sqrt(length(x))))
pdf('Number_ofDEGs_finalClustersMarkers_MAST.pdf', width = 8, height = 6)
par(las=2, mar=c(6, 5, 4, 3))
plot(seq(0.05, 0.95, 0.05), rev(unlist(lapply(DEGs_number, mean))/DEGs_number_all), pch=20, xaxt='n', xlab="", ylab="Fraction of conserved DEGs", ylim=c(0.48,1.02))
axis(1, at=seq(0.05, 0.95, 0.05), labels = rev(x), cex.lab=1.2, cex.axis=1.2)
lines(seq(0.05, 0.95, 0.05), avg)
title(xlab="Number of cells",mgp=c(4,1,0))
abline(v=0.75, col="red3", lty=2) #inflection point
arrows(seq(0.05, 0.95, 0.05), avg-rev(err), seq(0.05, 0.95, 0.05), avg+rev(err), length=0.05, angle=90, code=3)

DEGs_number<-list()
for (perc.sub in seq(0.05, 0.95, 0.05)){
  DEGs_number[[as.character(perc.sub)]]<-unlist(lapply(DEGs_list[[as.character(perc.sub)]], function(x) {
    x <- subset(x,pct.1 > 0.2 & avg_log2FC>log2(2))
    dim(x)[1]
  }))
}
DEGs_number_all<-dim(subset(majorClassMarkers,pct.1 > 0.2 & avg_log2FC>log2(2)))[1]
avg<-rev(unlist(lapply(DEGs_number, mean))/DEGs_number_all)
x<-round(seq(0.05, 0.95, 0.05)*51205)
err <- unlist(lapply(DEGs_number, function(x) sd(x/DEGs_number_all)/sqrt(length(x))))
points(seq(0.05, 0.95, 0.05), rev(unlist(lapply(DEGs_number, mean))/DEGs_number_all), pch=20, col="red4")
lines(seq(0.05, 0.95, 0.05), avg,col="red4")
arrows(seq(0.05, 0.95, 0.05), avg-rev(err), seq(0.05, 0.95, 0.05), avg+rev(err), length=0.05, angle=90, code=3, col = "red4")

DEGs_number<-list()
for (perc.sub in seq(0.05, 0.95, 0.05)){
  DEGs_number[[as.character(perc.sub)]]<-unlist(lapply(DEGs_list[[as.character(perc.sub)]], function(x) {
    x <- subset(x,pct.1 > 0.2 & avg_log2FC>log2(2.2))
    dim(x)[1]
  }))
}
DEGs_number_all<-dim(subset(majorClassMarkers,pct.1 > 0.2 & avg_log2FC>log2(2.2)))[1]
avg<-rev(unlist(lapply(DEGs_number, mean))/DEGs_number_all)
x<-round(seq(0.05, 0.95, 0.05)*51205)
err <- unlist(lapply(DEGs_number, function(x) sd(x/DEGs_number_all)/sqrt(length(x))))
points(seq(0.05, 0.95, 0.05), rev(unlist(lapply(DEGs_number, mean))/DEGs_number_all), pch=20, col="blue4")
lines(seq(0.05, 0.95, 0.05), avg,col="blue4")
arrows(seq(0.05, 0.95, 0.05), avg-rev(err), seq(0.05, 0.95, 0.05), avg+rev(err), length=0.05, angle=90, code=3, col = "blue4")

legend("bottom",x.intersp = c(0.5,0.5,0.5),pch = 20,lty=1, legend = c(expression(log[2](FC) > log[2](1.8)),expression(log[2](FC) > log[2](2.0)),expression(log[2](FC) > log[2](2.2))), col = c("black","red4","blue4"), horiz = TRUE,inset=c(0,1), xpd=TRUE, bty="n" )
dev.off()


## compute Robust Hausdorff distance 
library(BiocNeighbors)
rHD <- function(A, B, q=50, distance="Euclidean",
                BNPARAM=KmknnParam(distance=distance)){
  nn = queryKNN(A, B, k=1, get.index=FALSE, BNPARAM=BNPARAM)
  k = length(nn$distance) - q
  return(sort(nn$distance)[k])
}

nPCs<-readRDS('optimal_nPCs_integrated_filtered.RDS')
pca <- Embeddings(atlas_v2, reduction = "pca")[,1:nPCs]


HausdorffDistance <-list()
for (perc.sub in seq(0.05, 0.95, 0.05)){
  hd <- unlist(lapply(cell_subsampling[[as.character(perc.sub)]], function(x){
    pca_A <- pca[x,]
    pca_B <- pca[rownames(pca)[!(rownames(pca) %in% x)],]
    hd <- rHD(pca_A,pca_B)
    return(hd)
  })) 
  HausdorffDistance[[as.character(perc.sub)]] <- hd
}
saveRDS(HausdorffDistance,"HausdorffDistance_q50.RDS")

pdf('HausdorffDistance_of_subsampling_q50.pdf', width = 12)
plot(seq(0.05, 0.95, 0.05), rev(unlist(lapply(HausdorffDistance, mean))), pch=19, xaxt='n', xlab="Subsampling - % of cells", ylab="Robust Hausdorff Distance")
axis(1, at=seq(0.05, 0.95, 0.05), labels = rev(seq(0.05, 0.95, 0.05)), cex.lab=1.2, cex.axis=1.2)
dev.off()


x<-round(seq(0.05, 0.95, 0.05)*51205)
avg<-rev(unlist(lapply(HausdorffDistance, mean)))
err <- unlist(lapply(HausdorffDistance, function(x) sd(x)/sqrt(length(x))))
pdf('HausdorffDistance_of_subsampling_q50_v2.pdf', width = 8, height = 6)
par(las=2, mar=c(6, 5, 3, 3))
plot(seq(0.05, 0.95, 0.05), rev(unlist(lapply(HausdorffDistance, mean))), pch=20, xaxt='n', xlab="", ylab="Robust Hausdorff Distance")
axis(1, at=seq(0.05, 0.95, 0.05), labels = rev(x), cex.lab=1.2, cex.axis=1.2)
lines(seq(0.05, 0.95, 0.05), rev(unlist(lapply(HausdorffDistance, mean))))
title(xlab="Number of cells",mgp=c(4,1,0))
abline(v=0.8, col="red3",lty=2) #inflection point
arrows(seq(0.05, 0.95, 0.05), avg-rev(err), seq(0.05, 0.95, 0.05), avg+rev(err), length=0.05, angle=90, code=3)
dev.off()



########## Meta Neighbour on clusters
atlas_v2<-readRDS('Integrated_atlas_V2_filtered_annotated.RDS')
cell_subsampling <-readRDS('CellsSubsampling.RDS')
clusters<-readRDS("finalClusters_LRP.RDS")
clusters_LRP<-paste0(clusters,"_LRP")
names(clusters_LRP) <- names(clusters)
clusters<-readRDS("finalClusters_NM.RDS")
clusters_NM<-paste0(clusters,"_NM")
names(clusters_NM) <- names(clusters)
clusters<-readRDS("finalClusters_M.RDS")
clusters_M<-paste0(clusters,"_M")
names(clusters_M) <- names(clusters)
final_clusters <- c(clusters_LRP,clusters_NM,clusters_M)
atlas_v2$comparison_cluster_label <- final_clusters[colnames(atlas_v2)]
atlas_v2$experiment<-'Atlas_v2'
DefaultAssay(atlas_v2) <-'RNA'
Idents(atlas_v2) <-'comparison_cluster_label'
atlasSCE<-as.SingleCellExperiment(atlas_v2,assay='RNA')


number_conserved_clusters <-list()
for (perc.sub in seq(0.05, 0.95, 0.05)){
  print(perc.sub)
  n_clusters <- unlist(lapply(cell_subsampling[[as.character(perc.sub)]], function(x){
    data_sub<-subset(atlas_v2, cells=x[x %in% colnames(atlas_v2)])
    data_sub$experiment<-'subset'
    data_sub <- RenameCells(object = data_sub, add.cell.id = 'sub')
    subDataSCE<-as.SingleCellExperiment(data_sub,assay='RNA')
    var_genes = variableGenes(dat = subDataSCE, exp_labels = subDataSCE$experiment)
    trained_model = trainModel(var_genes = var_genes,dat = subDataSCE,study_id = subDataSCE$experiment,cell_type = subDataSCE$comparison_cluster_label)
    celltype_NV <- MetaNeighborUS(var_genes = var_genes,trained_model= trained_model,dat = atlasSCE,study_id = atlasSCE$experiment,cell_type = atlasSCE$comparison_cluster_label,fast_version = TRUE)
    n_clusters<-sum(diag(celltype_NV) > 0.9)
    return(n_clusters)
    gc()
  })) 
  number_conserved_clusters[[as.character(perc.sub)]] <- n_clusters
}
saveRDS(number_conserved_clusters,"n_clusters_conserved_finalClusters.RDS")

pdf('Number_of_conservedClusters_finalClustering.pdf', width = 12)
plot(seq(0.05, 0.95, 0.05), rev(unlist(lapply(number_conserved_clusters, mean)))/35, pch=19, xaxt='n', xlab="Subsampling - % of cells", ylab="Percentage of conserved clusters")
axis(1, at=seq(0.05, 0.95, 0.05), labels = rev(seq(0.05, 0.95, 0.05)), cex.lab=1.2, cex.axis=1.2)
dev.off()



conserved_clusters <-list()
for (perc.sub in seq(0.05, 0.95, 0.05)){
  print(perc.sub)
  n_clusters <- lapply(cell_subsampling[[as.character(perc.sub)]], function(x){
    data_sub<-subset(atlas_v2, cells=x[x %in% colnames(atlas_v2)])
    data_sub$experiment<-'subset'
    data_sub <- RenameCells(object = data_sub, add.cell.id = 'sub')
    subDataSCE<-as.SingleCellExperiment(data_sub,assay='RNA')
    var_genes = variableGenes(dat = subDataSCE, exp_labels = subDataSCE$experiment)
    trained_model = trainModel(var_genes = var_genes,dat = subDataSCE,study_id = subDataSCE$experiment,cell_type = subDataSCE$comparison_cluster_label)
    celltype_NV <- MetaNeighborUS(var_genes = var_genes,trained_model= trained_model,dat = atlasSCE,study_id = atlasSCE$experiment,cell_type = atlasSCE$comparison_cluster_label,fast_version = TRUE)
    n_clusters<-sum(diag(celltype_NV) > 0.9)
    return(celltype_NV)
    gc()
  }) 
  conserved_clusters[[as.character(perc.sub)]] <- n_clusters
}
saveRDS(conserved_clusters,"MN_clusters_conserved_finalClusters.RDS")


###Plot
cluster_number<-list()
for (perc.sub in seq(0.05, 0.95, 0.05)){
  cluster_number[[as.character(perc.sub)]]<-unlist(lapply(conserved_clusters[[as.character(perc.sub)]], function(x) {
    sum(diag(x) >0.9)
  }))
}
cluster_number_all<-35

avg<-rev(unlist(lapply(cluster_number, mean))/cluster_number_all)
x<-round(seq(0.05, 0.95, 0.05)*51205)
err <- unlist(lapply(cluster_number, function(x) sd(x/cluster_number_all)/sqrt(length(x))))
pdf('Number_ofcluster_finalClustersMarkers_MAST_v2.pdf', width = 8, height = 6)
par(las=2, mar=c(6, 5, 4, 3))
plot(seq(0.05, 0.95, 0.05), rev(unlist(lapply(cluster_number, mean))/cluster_number_all), pch=20, xaxt='n', xlab="", ylab="Fraction of conserved cluster", ylim=c(0.54,1.02))
axis(1, at=seq(0.05, 0.95, 0.05), labels = rev(x), cex.lab=1.2, cex.axis=1.2)
lines(seq(0.05, 0.95, 0.05), avg)
title(xlab="Number of cells",mgp=c(4,1,0))
abline(v=0.8, col="red3", lty=2) #inflection point
arrows(seq(0.05, 0.95, 0.05), avg-rev(err), seq(0.05, 0.95, 0.05), avg+rev(err), length=0.05, angle=90, code=3)

cluster_number<-list()
for (perc.sub in seq(0.05, 0.95, 0.05)){
  cluster_number[[as.character(perc.sub)]]<-unlist(lapply(conserved_clusters[[as.character(perc.sub)]], function(x) {
    sum(diag(x) >0.88)
  }))
}
avg<-rev(unlist(lapply(cluster_number, mean))/cluster_number_all)
x<-round(seq(0.05, 0.95, 0.05)*51205)
err <- unlist(lapply(cluster_number, function(x) sd(x/cluster_number_all)/sqrt(length(x))))
points(seq(0.05, 0.95, 0.05), rev(unlist(lapply(cluster_number, mean))/cluster_number_all), pch=20, col="red4")
lines(seq(0.05, 0.95, 0.05), avg,col="red4")
arrows(seq(0.05, 0.95, 0.05), avg-rev(err), seq(0.05, 0.95, 0.05), avg+rev(err), length=0.05, angle=90, code=3, col = "red4")

cluster_number<-list()
for (perc.sub in seq(0.05, 0.95, 0.05)){
  cluster_number[[as.character(perc.sub)]]<-unlist(lapply(conserved_clusters[[as.character(perc.sub)]], function(x) {
    sum(diag(x) >0.86)
  }))
}
avg<-rev(unlist(lapply(cluster_number, mean))/cluster_number_all)
x<-round(seq(0.05, 0.95, 0.05)*51205)
err <- unlist(lapply(cluster_number, function(x) sd(x/cluster_number_all)/sqrt(length(x))))
points(seq(0.05, 0.95, 0.05), rev(unlist(lapply(cluster_number, mean))/cluster_number_all), pch=20, col="blue4")
lines(seq(0.05, 0.95, 0.05), avg,col="blue4")
arrows(seq(0.05, 0.95, 0.05), avg-rev(err), seq(0.05, 0.95, 0.05), avg+rev(err), length=0.05, angle=90, code=3, col = "blue4")

legend("bottom",x.intersp = c(0.5,0.5,0.5),pch = 20,lty=1, legend = c(expression(AUROC > 0.9),expression(AUROC > 0.88),expression(AUROC > 0.86)), col = c("black","red4","blue4"), horiz = TRUE,inset=c(0,1), xpd=TRUE, bty="n" )
dev.off()


#### Difference on mean expression within clusters
atlas_v2<-readRDS('Integrated_atlas_V2_filtered_annotated.RDS')
DefaultAssay(atlas_v2) <- 'RNA'
Idents(atlas_v2) <- 'final_clusters'
DEGenes<-FindAllMarkers(atlas_v2, logfc.threshold = log2(2), base=exp(1), only.pos=T)
DEGenes<-DEGenes[which(DEGenes$p_val_adj<0.01),]

mean_expr <- AverageExpression(atlas_v2,features = DEGenes$gene,assays = 'RNA')

diff_mean <-list()
for (perc.sub in seq(0.05, 0.95, 0.05)){
  diff <- unlist(lapply(cell_subsampling[[as.character(perc.sub)]], function(x){
    data_sub<-subset(atlas_v2, cells=x[x %in% colnames(atlas_v2)])
    Idents(data_sub)  <- 'final_clusters'
    mean_expr_sub <- AverageExpression(data_sub,features = DEGenes$gene,assays = 'RNA')
    #cor<-WGCNA::cor(t(mean_expr$RNA[rownames(mean_expr_sub$RNA),colnames(mean_expr_sub$RNA)]),t(mean_expr_sub$RNA))
    diff<-sum(rowMeans(abs(mean_expr$RNA[rownames(mean_expr_sub$RNA),colnames(mean_expr_sub$RNA)] - mean_expr_sub$RNA)))
    return(diff)
  })) 
  diff_mean[[as.character(perc.sub)]] <- diff
}
saveRDS(diff_mean,"diff_mean_subsampling_DEGs.RDS")

x<-round(seq(0.05, 0.95, 0.05)*51205)
avg<-rev(unlist(lapply(diff_mean, mean)))
err <- unlist(lapply(diff_mean, function(x) sd(x)/sqrt(length(x))))
pdf('Differences_MeanExpression_of_finalClustersMarkers_MAST_v2.pdf', width = 8, height = 6)
par(las=2, mar=c(6, 5, 3, 3))
plot(seq(0.05, 0.95, 0.05), rev(unlist(lapply(diff_mean, mean))), pch=20, xaxt='n', xlab="", ylab="Absolute Difference (mean expression)")
axis(1, at=seq(0.05, 0.95, 0.05), labels = rev(x), cex.lab=1.2, cex.axis=1.2)
lines(seq(0.05, 0.95, 0.05), rev(unlist(lapply(diff_mean, mean))))
title(xlab="Number of cells",mgp=c(4,1,0))
abline(v=0.75, col="red3", lty=2) #inflection point
arrows(seq(0.05, 0.95, 0.05), avg-rev(err), seq(0.05, 0.95, 0.05), avg+rev(err), length=0.05, angle=90, code=3)
dev.off()




