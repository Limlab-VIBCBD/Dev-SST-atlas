# Run this after calculating clusters with the iterative clustering pipeline (nextflow)
# Check whether any pairs of clusters should be merged based on marker gene expression.
library(Seurat)
library(ggplot2)
library(dplyr)
library(ComplexHeatmap)
# load output from nextflow clustering pipeline
obj<-readRDS('object_with_final_clusters.RDS')
DefaultAssay(obj) <-'RNA'
Idents(obj) <- 'final_clusters'
# Compute markers
markers <- FindAllMarkers(obj,logfc.threshold = log2(1.5), min.pct = 0.2,only.pos = TRUE, test.use = 'MAST', latent.vars = 'batch')
ribo_genes <- grep(pattern = "^Rp[sl]", x = rownames(obj@assays$RNA@counts), value = TRUE)
markers<-markers[which(markers$gene %!in% ribo_genes),]
mt_genes <- grep(pattern = "^mt-", x = rownames(obj@assays$RNA@counts), value = TRUE)
markers<-markers[which(markers$gene %!in% mt_genes),]
markers<-markers[order(markers$avg_log2FC, decreasing = TRUE),]
markers<-markers[which(markers$p_val_adj < 0.01),]
markers<-as.data.frame(markers)
markers<-markers[order(markers$cluster),]
# Report the average expression of marker genes for each cluster and compute hierachical clustering on clusters
p <- DotPlot(object = obj, assay="RNA", features=unique(markers$gene),scale=TRUE) + ylab("Clusters") + scale_colour_gradient2(low = "blue2", mid = "gray90", high = "red2", midpoint=0) + theme(axis.text.x = element_text(angle = 90, hjust = 1))
df<- p$data
exp_mat<-df %>% select(-pct.exp, -avg.exp) %>%  pivot_wider(names_from = id, values_from = avg.exp.scaled) %>% as.data.frame() 
row.names(exp_mat) <- exp_mat$features.plot  
exp_mat <- exp_mat[,-1] %>% as.matrix()
percent_mat<-df %>% select(-avg.exp, -avg.exp.scaled) %>%  pivot_wider(names_from = id, values_from = pct.exp) %>% as.data.frame() 
row.names(percent_mat) <- percent_mat$features.plot  
percent_mat <- percent_mat[,-1] %>% as.matrix()
col_fun = circlize::colorRamp2(c(-1.5, 0, 2.5), c("blue2","gray90","red2"))
layer_fun = function(j, i, x, y, w, h, fill){
  grid.rect(x = x, y = y, width = w, height = h, 
            gp = gpar(col = NA, fill = NA))
  grid.circle(x=x,y=y,r= pindex(percent_mat, i, j)/100 * unit(1.5, "mm"),
              gp = gpar(fill = col_fun(pindex(exp_mat, i, j)), col = NA))}
lgd_list = list(
  Legend( labels = c(0,0.25,0.5,0.75,1), title = "pt",
          graphics = list(
            function(x, y, w, h) grid.circle(x = x, y = y, r = 0 * unit(1.5, "mm"),
                                             gp = gpar(fill = "black")),
            function(x, y, w, h) grid.circle(x = x, y = y, r = 0.25 * unit(1.5, "mm"),
                                             gp = gpar(fill = "black")),
            function(x, y, w, h) grid.circle(x = x, y = y, r = 0.5 * unit(1.5, "mm"),
                                             gp = gpar(fill = "black")),
            function(x, y, w, h) grid.circle(x = x, y = y, r = 0.75 * unit(1.5, "mm"),
                                             gp = gpar(fill = "black")),
            function(x, y, w, h) grid.circle(x = x, y = y, r = 1 * unit(1.5, "mm"),
                                             gp = gpar(fill = "black")))
  ))
set.seed(123)    
hp<- Heatmap(exp_mat,
             heatmap_legend_param=list(title="expression"),
             column_title = "clustered dotplot", 
             col=col_fun,
             rect_gp = gpar(type = "none"),
             layer_fun = layer_fun,
             row_names_gp = gpar(fontsize = 7),
             border = "black",
             cluster_rows = FALSE)
#cluster_order<-hp@column_names_param$labels[column_order(hp)]
print(hp)
hp = draw(hp)
cluster_order <- column_order(hp)
tmp<-markers %>% arrange(factor(cluster, levels = unique(markers$cluster)[cluster_order]))
ordered_genes <- unique(tmp$gene)
exp_mat<-exp_mat[ordered_genes,]
percent_mat<-percent_mat[ordered_genes,]
hp<- Heatmap(exp_mat,
             heatmap_legend_param=list(title="expression"),
             column_title = "clustered dotplot", 
             col=col_fun,
             rect_gp = gpar(type = "none"),
             layer_fun = layer_fun,
             row_names_gp = gpar(fontsize = 7),
             border = "black",
             cluster_rows = FALSE)

# Identify couples with cophenetic distance below the 5th percentile of the distribution of cophenetic distances across all cluster pairs
hp = draw(hp)
dend_dist<-cophenetic(column_dend(hp))
cutoff<-quantile(as.vector(dend_dist), 0.05)

dend_dist <- as.matrix(dend_dist)
dend_dist[lower.tri(dend_dist,diag=TRUE)] <- 0     
clusters_to_test<-lapply(rownames(dend_dist), function(x){
  clusters <- colnames(dend_dist)[dend_dist[x,] < cutoff & dend_dist[x,] != 0]
  return(clusters)
})
names(clusters_to_test) <- rownames(dend_dist)
clusters_to_test <- clusters_to_test[lapply(clusters_to_test,length)>0]

# Compute DEscore between identified couples
DefaultAssay(obj) <-'RNA'
Idents(obj) <- 'final_clusters'
ribo_genes <- grep(pattern = "^Rp[sl]", x = rownames(obj@assays$RNA@counts), value = TRUE)
clusters_to_merge<-lapply(names(clusters_to_test), function(x){
  to_merge<-c()
  for (i in 1:length(clusters_to_test[[x]])){
    markers<-FindMarkers(obj, ident.1 = x, ident.2 = clusters_to_test[[x]][i], only.pos = TRUE, assay = 'RNA',logfc.threshold = 1, verbose = FALSE)
    markers<-markers[which(markers$p_val_adj<=0.01),]
    markers<-markers[which(rownames(markers) %!in% ribo_genes),]
    markers2<-FindMarkers(obj, ident.1 = clusters_to_test[[x]][i], ident.2 = x, only.pos = TRUE, assay = 'RNA',logfc.threshold = 1, verbose = FALSE)
    markers2<-markers2[which(markers2$p_val_adj<=0.01),]
    markers2<-markers2[which(rownames(markers2) %!in% ribo_genes),]
    if(nrow(markers)==0){
      DEscore <- 0
    }else{
      tmp<- log10(markers$p_val_adj)*(-1)
      tmp[which(tmp>20)]<-20
      DEscore<-sum(tmp)
      print(DEscore)
    }
    if(nrow(markers2)==0){
      DEscore2 <- 0
    }else{
      tmp<- log10(markers2$p_val_adj)*(-1)
      tmp[which(tmp>20)]<-20
      DEscore2<-sum(tmp)
      print(DEscore2)
    }
    if(DEscore < 60 | DEscore2 < 60){
      to_merge<-c(to_merge,clusters_to_test[[x]][i])
    }
  }
  return(to_merge)
})
names(clusters_to_merge) <- names(clusters_to_test)
# clusters_to_merge contains information on clusters that should be merged

