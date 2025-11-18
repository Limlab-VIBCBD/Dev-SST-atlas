########################
# Load required packages
########################
library(Seurat)
library(anndata)
library(SeuratDisk)
library(SeuratWrappers)
library(readr)
library(rWikiPathways)
library(gplots)
library(grDevices)

#############
# WOT on LRPs
#############
# Create directory for LRPs
# dir.create('LRPs')
# setwd('LRPs')
# load LRP subset with clusters
atlas_v2_LRP<-subset(atlas_v2_LRP, subset=final_clusters_renamed %in% c("LRP2.1","LRP1.0"))  
atlas_v2_LRP<-subset(atlas_v2_LRP, subset=time_point %in% c("E16","P1","P5"))     
atlas_v2_LRP$time <- parse_number(atlas_v2_LRP$time_point)
atlas_v2_LRP$time[which(atlas_v2_LRP$time == 16)] <- -4
atlas_v2_LRP$time <- atlas_v2_LRP$time + 4 # set first time point to 0
# Prepare files for WOT analysis
counts<-atlas_v2_LRP@assays$RNA$data
tmp<-CreateSeuratObject(counts)
var_genes <- rownames(atlas_v2_LRP@assays$integrated@data)
counts <- counts[var_genes,]
tmp<-CreateSeuratObject(counts)
SaveH5Seurat(tmp, filename = "matrix_varGenes.h5Seurat", overwrite = TRUE)
Convert("matrix_varGenes.h5Seurat", dest = "h5ad")
cell_days <- data.frame(id=rownames(atlas_v2_LRP@meta.data),day=atlas_v2_LRP@meta.data$time)
write.table(cell_days,"cell_days.txt", sep="\t", quote=FALSE, col.names = TRUE, row.names = FALSE)
embedding_coord<- Embeddings(atlas_v2_LRP,reduction = "umap")
embedding_coord<-data.frame(id=rownames(embedding_coord),x=embedding_coord[,"UMAP_1"],y=embedding_coord[,"UMAP_2"])
write.table(embedding_coord,"embedding_coord.txt", sep="\t", quote=FALSE, col.names = TRUE, row.names = FALSE)
df<-data.frame(cl_renamed = atlas_v2_LRP$final_clusters_renamed, cl=atlas_v2_LRP$final_clusters, id=rownames(atlas_v2_LRP@meta.data))
writeGMT(df, "cell_sets.gmt")
batches <- data.frame(id=rownames(atlas_v2_LRP@meta.data),batch=atlas_v2_LRP@meta.data$batch)
write.table(batches,"batches.txt", sep="\t", quote=FALSE, col.names = TRUE, row.names = FALSE)

# Run WOT using command line interface:
## wot optimal_transport --matrix matrix_varGenes.h5ad --cell_days cell_days.txt --growth_iters 3 --lambda1 1 --lambda2 50 --epsilon 0.05 --verbose
## wot trajectory --tmap tmaps --cell_set cell_sets.gmt --day 9 --embedding embedding_coord.txt
## wot transition_table --tmap tmaps --cell_set cell_sets.gmt --start_time 0 --end_time 9
## wot fates --tmap tmaps --cell_set cell_sets.gmt --day 5 --out to_P1
## wot fates --tmap tmaps --cell_set cell_sets.gmt --day 9 --out to_P5
## wot optimal_transport_validation --matrix matrix_varGenes.h5ad --cell_days cell_days.txt --covariate batches.txt --cell_growth_rates tmaps_g.txt --cell_growth_rates_field g2 --covariate_field batch --verbose

# Import WOT analysis outputs
# Fates to P1 
# Probability of each cell transitioning into a given cluster identity at subsequent time points
fates_to_p1<- read.table('to_P1_3_fates.txt', header = TRUE)
rownames(fates_to_p1) <- fates_to_p1$id
fates_to_p1$cluster <- atlas_v2_LRP@meta.data[rownames(fates_to_p1),"final_clusters_renamed"]
fates_to_p1$age <- atlas_v2_LRP@meta.data[rownames(fates_to_p1),"time_point"]
plot_to_p1 <- reshape2::melt(fates_to_p1, value.name = "Probability")
colnames(plot_to_p1) <- c('id','cluster',"age","Fate","Probability")
plot_to_p1_e16<- subset(plot_to_p1, subset=age=="E16")
# for each cluster we calculate the median fate probabilities toward all possible cluster identities at the following time point
by_cluster <- split(plot_to_p1_e16,plot_to_p1_e16$cluster)
median_table <- Reduce("rbind",lapply(by_cluster, function(data){
  by_fate<-split(data,data$Fate)
  median<-unlist(lapply(by_fate, function(data_sub){
    median(data_sub$Probability)
  }))
  return(median)
}))
rownames(median_table) <- names(by_cluster)
median_table<-median_table[,rownames(median_table)]
saveRDS(median_table,'median_probability_table_E16_to_P1.RDS')
median_table <-round(median_table,digits = 2)
min_val <- min(median_table)
max_val <- max(median_table)
breaks <-  seq(min_val, max_val, length.out=1001);
cols <- colorRampPalette(colors = c('black','#2D718EFF','#3CBC75FF','gold'))(1000)
pdf(file="Median_probability_LRP_E16_to_P1.pdf",width=5,height = 5)
gplots::heatmap.2(median_table,
                  margins=c(5,5),
                  key = FALSE, # Explicitly set key to TRUE
                  keysize=1,
                  key.xlab="",
                  key.title="Probability (median) - from E16 to P5",
                  trace = "none",
                  density.info = "none",
                  col = cols,
                  breaks = breaks,
                  offsetRow=0.1,
                  offsetCol=0.1,
                  cexRow = 0.8,
                  cexCol = 0.8,
                  cellnote = median_table,
                  notecex = 0.9,
                  notecol = 'white',
                  Colv = F,
                  Rowv = F,
                  dendrogram = "none")
dev.off()
# Fates to P5 
# Probability of each cell transitioning into a given cluster identity at subsequent time points
fates_to_p5<- read.table('to_P5_fates.txt', header = TRUE)
rownames(fates_to_p5) <- fates_to_p5$id
fates_to_p5$cluster <- obj@meta.data[rownames(fates_to_p5),"final_clusters_renamed"]
fates_to_p5$age <- obj@meta.data[rownames(fates_to_p5),"time_point"]
plot_to_p5 <- reshape2::melt(fates_to_p5, value.name = "Probability")
colnames(plot_to_p5) <- c('id','cluster',"age","Fate","Probability")
# we extract P1 to P5 probabilities
plot_to_p5_p1<- subset(plot_to_p5, subset=age=="P1")
# for each cluster we calculate the median fate probabilities toward all possible cluster identities at the following time point
by_cluster <- split(plot_to_p5_p1,plot_to_p5_p1$cluster)
median_table <- Reduce("rbind",lapply(by_cluster, function(data){
  by_fate<-split(data,data$Fate)
  median<-unlist(lapply(by_fate, function(data_sub){
    median(data_sub$Probability)
  }))
  return(median)
}))
rownames(median_table) <- names(by_cluster)
median_table<-median_table[,rownames(median_table)]
saveRDS(median_table,'median_probability_table_P1_to_P5.RDS')
median_table <-round(median_table,digits = 2)
min_val <- min(median_table)
max_val <- max(median_table)
breaks <-  seq(min_val, max_val, length.out=1001);
cols <- colorRampPalette(colors = c('black','#2D718EFF','#3CBC75FF','gold'))(1000)
pdf(file="Median_probability_LRP_P1_to_P5.pdf",width=5,height = 5)
gplots::heatmap.2(median_table,
                  margins=c(5,5),
                  key = FALSE, # Explicitly set key to TRUE
                  keysize=1,
                  key.xlab="",
                  key.title="Probability (median) - from E16 to P5",
                  trace = "none",
                  density.info = "none",
                  col = cols,
                  breaks = breaks,
                  offsetRow=0.1,
                  offsetCol=0.1,
                  cexRow = 0.8,
                  cexCol = 0.8,
                  cellnote = median_table,
                  notecex = 0.9,
                  notecol = 'white',
                  Colv = F,
                  Rowv = F,
                  dendrogram = "none")
dev.off()


###################
# WOT on Martinotti
###################
# Create directory for Martinotti
# dir.create('Martinotti')
# setwd('Martinotti')
# load Martinotti subset with clusters
atlas_v2_Martinotti<-subset(atlas_v2_Martinotti, subset=time_point %in% c("E16","P1","P5"))     
atlas_v2_Martinotti$time <- parse_number(atlas_v2_Martinotti$time_point)
atlas_v2_Martinotti$time[which(atlas_v2_Martinotti$time == 16)] <- -4
atlas_v2_Martinotti$time <- atlas_v2_Martinotti$time + 4
# As cluster sizes is imbalanced, we computed the median cluster size across all clusters and randomly subsampled this number of cells from each cluster. This approach allowed us to generate a more balanced dataset for fate prediction. 
# The subsampling procedure was repeated three times, and WOT analysis was performed independently on each random subset. Subsamples are repeated to ensure a complete representation of the dataset.
for (sub in 1:3){
  dir.create(paste0('cell_subset_',sub))
  setwd(paste0('cell_subset_',sub))
  n_cell<-median(table(atlas_v2_Martinotti$final_clusters_renamed))
  select_cells <- unlist(lapply(unique(obj$final_clusters_renamed), function(cl){
      cells<-rownames(obj@meta.data[which(obj@meta.data$final_clusters_renamed == cl),])
      if (length(cells) >= n_cell){
          set.seed(345+sub)
          cells <- sample(cells,size = n_cell,replace = FALSE)
      }
      return(cells)
  }))
  obj<-subset(atlas_v2_Martinotti, cells = select_cells) 
  # Prepare files for WOT analysis
  counts<-obj@assays$RNA$data
  tmp<-CreateSeuratObject(counts)
  var_genes <- rownames(obj@assays$integrated@data)
  counts <- counts[var_genes,]
  tmp<-CreateSeuratObject(counts)
  SaveH5Seurat(tmp, filename = "matrix_varGenes.h5Seurat", overwrite = TRUE)
  Convert("matrix_varGenes.h5Seurat", dest = "h5ad")
  cell_days <- data.frame(id=rownames(obj@meta.data),day=obj@meta.data$time)
  write.table(cell_days,"cell_days.txt", sep="\t", quote=FALSE, col.names = TRUE, row.names = FALSE)
  embedding_coord<- Embeddings(obj,reduction = "umap")
  embedding_coord<-data.frame(id=rownames(embedding_coord),x=embedding_coord[,"UMAP_1"],y=embedding_coord[,"UMAP_2"])
  write.table(embedding_coord,"embedding_coord.txt", sep="\t", quote=FALSE, col.names = TRUE, row.names = FALSE)
  df<-data.frame(cl_renamed = obj$final_clusters_renamed, cl=obj$final_clusters, id=rownames(obj@meta.data))
  writeGMT(df, "cell_sets.gmt")
  batches <- data.frame(id=rownames(obj@meta.data),batch=obj@meta.data$batch)
  write.table(batches,"batches.txt", sep="\t", quote=FALSE, col.names = TRUE, row.names = FALSE)

  # Run WOT using command line interface:
  ## wot optimal_transport --matrix matrix_varGenes.h5ad --cell_days cell_days.txt --growth_iters 3 --lambda1 1 --lambda2 50 --epsilon 0.05 --verbose
  ## wot trajectory --tmap tmaps --cell_set cell_sets.gmt --day 9 --embedding embedding_coord.txt
  ## wot transition_table --tmap tmaps --cell_set cell_sets.gmt --start_time 0 --end_time 9
  ## wot fates --tmap tmaps --cell_set cell_sets.gmt --day 5 --out to_P1
  ## wot fates --tmap tmaps --cell_set cell_sets.gmt --day 9 --out to_P5
  ## wot optimal_transport_validation --matrix matrix_varGenes.h5ad --cell_days cell_days.txt --covariate batches.txt --cell_growth_rates tmaps_g.txt --cell_growth_rates_field g2 --covariate_field batch --verbose
  
  # Import WOT analysis outputs
  # Fates to P1 
  # Probability of each cell transitioning into a given cluster identity at subsequent time points
  fates_to_p1<- read.table('to_P1_3_fates.txt', header = TRUE)
  rownames(fates_to_p1) <- fates_to_p1$id
  fates_to_p1$cluster <- obj@meta.data[rownames(fates_to_p1),"final_clusters_renamed"]
  fates_to_p1$age <- obj@meta.data[rownames(fates_to_p1),"time_point"]
  plot_to_p1 <- reshape2::melt(fates_to_p1, value.name = "Probability")
  colnames(plot_to_p1) <- c('id','cluster',"age","Fate","Probability")
  plot_to_p1_e16<- subset(plot_to_p1, subset=age=="E16")
  # for each cluster we calculate the median fate probabilities toward all possible cluster identities at the following time point
  by_cluster <- split(plot_to_p1_e16,plot_to_p1_e16$cluster)
  median_table <- Reduce("rbind",lapply(by_cluster, function(data){
  by_fate<-split(data,data$Fate)
  median<-unlist(lapply(by_fate, function(data_sub){
      median(data_sub$Probability)
  }))
  return(median)
  }))
  rownames(median_table) <- names(by_cluster)
  median_table<-median_table[,rownames(median_table)]
  median_table <-round(median_table,digits = 2)
  saveRDS(median_table,'median_probability_table_E16_to_P1.RDS')
  # Fates to P5 
  # Probability of each cell transitioning into a given cluster identity at subsequent time points
  fates_to_p5<- read.table('to_P5_fates.txt', header = TRUE)
  rownames(fates_to_p5) <- fates_to_p5$id
  fates_to_p5$cluster <- obj@meta.data[rownames(fates_to_p5),"final_clusters_renamed"]
  fates_to_p5$age <- obj@meta.data[rownames(fates_to_p5),"time_point"]
  plot_to_p5 <- reshape2::melt(fates_to_p5, value.name = "Probability")
  colnames(plot_to_p5) <- c('id','cluster',"age","Fate","Probability")
  # we extract P1 to P5 probabilities
  plot_to_p5_p1<- subset(plot_to_p5, subset=age=="P1")
  # for each cluster we calculate the median fate probabilities toward all possible cluster identities at the following time point
  by_cluster <- split(plot_to_p5_p1,plot_to_p5_p1$cluster)
  median_table <- Reduce("rbind",lapply(by_cluster, function(data){
  by_fate<-split(data,data$Fate)
  median<-unlist(lapply(by_fate, function(data_sub){
      median(data_sub$Probability)
  }))
  return(median)
  }))
  rownames(median_table) <- names(by_cluster)
  median_table<-median_table[,rownames(median_table)]
  median_table <-round(median_table,digits = 2)
  saveRDS(median_table,'median_probability_table_P1_to_P5.RDS')
}
# Merge subsamplings 
# E16 to P1
median_1 <-readRDS('cell_subset_1/median_probability_table_E16_to_P1.RDS')
median_2 <-readRDS('cell_subset_2/median_probability_table_E16_to_P1.RDS')
median_3 <-readRDS('cell_subset_3/median_probability_table_E16_to_P1.RDS')
# Compute mean across subsamplings
clusters<-rownames(median_1)
median<-matrix(ncol = length(clusters), nrow = length(clusters))
rownames(median) <- colnames(median) <- clusters
for (i in clusters){
  for (j in clusters){
    median[i,j] <- mean(c(median_1[i,j], median_2[i,j], median_3[i,j]))
  }
}
saveRDS(median, 'MedianProbability_table_Martinotti_finalClusters_E16_to_P1.RDS')
median_table <-round(median,digits = 2)
min_val <- min(median_table)
max_val <- max(median_table)
breaks <-  seq(min_val, max_val, length.out=1001);
cols <- colorRampPalette(colors = c('black','#2D718EFF','#3CBC75FF','gold'))(1000)
pdf(file="Median_probability_Martinotti_final_clusters_E16_to_P1.pdf",width=12,height = 12)
gplots::heatmap.2(median_table,
                  margins=c(8,8),
                  key = FALSE, # Explicitly set key to TRUE
                  keysize=1,
                  key.xlab="",
                  key.title="Probability (median) - from E16 to P5",
                  trace = "none",
                  density.info = "none",
                  col = cols,
                  breaks = breaks,
                  offsetRow=0.1,
                  offsetCol=0.1,
                  cexRow = 0.8,
                  cexCol = 0.8,
                  cellnote = median_table,
                  notecex = 0.9,
                  notecol = 'white',
                  Colv = F,
                  Rowv = F,
                  dendrogram = "none")
dev.off()
# P1 to P5
median_1 <-readRDS('cell_subset_1/median_probability_table_P1_to_P5.RDS')
median_2 <-readRDS('cell_subset_2/median_probability_table_P1_to_P5.RDS')
median_3 <-readRDS('cell_subset_3/median_probability_table_P1_to_P5.RDS')
# Compute mean across subsamplings
clusters<-rownames(median_1)
median<-matrix(ncol = length(clusters), nrow = length(clusters))
rownames(median) <- colnames(median) <- clusters
for (i in clusters){
  for (j in clusters){
    median[i,j] <- mean(c(median_1[i,j], median_2[i,j], median_3[i,j]))
  }
}
saveRDS(median, 'MedianProbability_table_Martinotti_finalClusters_P1_to_P5.RDS')
median_table <-round(median,digits = 2)
min_val <- min(median_table)
max_val <- max(median_table)
#breaks <-  seq(0, 0.5, length.out=1001);
breaks <-  seq(min_val, max_val, length.out=1001);
cols <- colorRampPalette(colors = c('black','#2D718EFF','#3CBC75FF','gold'))(1000)
pdf(file="Median_probability_Martinotti_final_clusters_P1_to_P5.pdf",width=12,height = 12)
gplots::heatmap.2(median_table,
                  margins=c(8,8),
                  key = FALSE, # Explicitly set key to TRUE
                  keysize=1,
                  key.xlab="",
                  key.title="Probability (median) - from E16 to P5",
                  trace = "none",
                  density.info = "none",
                  col = cols,
                  breaks = breaks,
                  offsetRow=0.1,
                  offsetCol=0.1,
                  cexRow = 0.8,
                  cexCol = 0.8,
                  cellnote = median_table,
                  notecex = 0.9,
                  notecol = 'white',
                  Colv = F,
                  Rowv = F,
                  dendrogram = "none")
dev.off()


######################
# WOT on NonMartinotti
######################
# Create directory for NonMartinotti
# dir.create('NonMartinotti')
# setwd('NonMartinotti')
# load Non Martinotti subset with clusters
atlas_v2_NonMartinotti<-subset(atlas_v2_NonMartinotti, subset=time_point %in% c("E16","P1","P5"))     
atlas_v2_NonMartinotti$time <- parse_number(atlas_v2_NonMartinotti$time_point)
atlas_v2_NonMartinotti$time[which(atlas_v2_NonMartinotti$time == 16)] <- -4
atlas_v2_NonMartinotti$time <- atlas_v2_NonMartinotti$time + 4
# As cluster sizes is imbalanced, we computed the median cluster size across all clusters and randomly subsampled this number of cells from each cluster. This approach allowed us to generate a more balanced dataset for fate prediction. 
# The subsampling procedure was repeated three times, and WOT analysis was performed independently on each random subset. Subsamples are repeated to ensure a complete representation of the dataset.
for (sub in 1:3){
  dir.create(paste0('cell_subset_',sub))
  setwd(paste0('cell_subset_',sub))
  n_cell<-median(table(atlas_v2_NonMartinotti$final_clusters_renamed))
  select_cells <- unlist(lapply(unique(obj$final_clusters_renamed), function(cl){
      cells<-rownames(obj@meta.data[which(obj@meta.data$final_clusters_renamed == cl),])
      if (length(cells) >= n_cell){
          set.seed(345+sub)
          cells <- sample(cells,size = n_cell,replace = FALSE)
      }
      return(cells)
  }))
  obj<-subset(atlas_v2_NonMartinotti, cells = select_cells) 
  # Prepare files for WOT analysis
  counts<-obj@assays$RNA$data
  tmp<-CreateSeuratObject(counts)
  var_genes <- rownames(obj@assays$integrated@data)
  counts <- counts[var_genes,]
  tmp<-CreateSeuratObject(counts)
  SaveH5Seurat(tmp, filename = "matrix_varGenes.h5Seurat", overwrite = TRUE)
  Convert("matrix_varGenes.h5Seurat", dest = "h5ad")
  cell_days <- data.frame(id=rownames(obj@meta.data),day=obj@meta.data$time)
  write.table(cell_days,"cell_days.txt", sep="\t", quote=FALSE, col.names = TRUE, row.names = FALSE)
  embedding_coord<- Embeddings(obj,reduction = "umap")
  embedding_coord<-data.frame(id=rownames(embedding_coord),x=embedding_coord[,"UMAP_1"],y=embedding_coord[,"UMAP_2"])
  write.table(embedding_coord,"embedding_coord.txt", sep="\t", quote=FALSE, col.names = TRUE, row.names = FALSE)
  df<-data.frame(cl_renamed = obj$final_clusters_renamed, cl=obj$final_clusters, id=rownames(obj@meta.data))
  writeGMT(df, "cell_sets.gmt")
  batches <- data.frame(id=rownames(obj@meta.data),batch=obj@meta.data$batch)
  write.table(batches,"batches.txt", sep="\t", quote=FALSE, col.names = TRUE, row.names = FALSE)

  # Run WOT using command line interface:
  ## wot optimal_transport --matrix matrix_varGenes.h5ad --cell_days cell_days.txt --growth_iters 3 --lambda1 1 --lambda2 50 --epsilon 0.05 --verbose
  ## wot trajectory --tmap tmaps --cell_set cell_sets.gmt --day 9 --embedding embedding_coord.txt
  ## wot transition_table --tmap tmaps --cell_set cell_sets.gmt --start_time 0 --end_time 9
  ## wot fates --tmap tmaps --cell_set cell_sets.gmt --day 5 --out to_P1
  ## wot fates --tmap tmaps --cell_set cell_sets.gmt --day 9 --out to_P5
  ## wot optimal_transport_validation --matrix matrix_varGenes.h5ad --cell_days cell_days.txt --covariate batches.txt --cell_growth_rates tmaps_g.txt --cell_growth_rates_field g2 --covariate_field batch --verbose
  
  # Import WOT analysis outputs
  # Fates to P1 
  # Probability of each cell transitioning into a given cluster identity at subsequent time points
  fates_to_p1<- read.table('to_P1_3_fates.txt', header = TRUE)
  rownames(fates_to_p1) <- fates_to_p1$id
  fates_to_p1$cluster <- obj@meta.data[rownames(fates_to_p1),"final_clusters_renamed"]
  fates_to_p1$age <- obj@meta.data[rownames(fates_to_p1),"time_point"]
  plot_to_p1 <- reshape2::melt(fates_to_p1, value.name = "Probability")
  colnames(plot_to_p1) <- c('id','cluster',"age","Fate","Probability")
  plot_to_p1_e16<- subset(plot_to_p1, subset=age=="E16")
  # for each cluster we calculate the median fate probabilities toward all possible cluster identities at the following time point
  by_cluster <- split(plot_to_p1_e16,plot_to_p1_e16$cluster)
  median_table <- Reduce("rbind",lapply(by_cluster, function(data){
  by_fate<-split(data,data$Fate)
  median<-unlist(lapply(by_fate, function(data_sub){
      median(data_sub$Probability)
  }))
  return(median)
  }))
  rownames(median_table) <- names(by_cluster)
  median_table<-median_table[,rownames(median_table)]
  median_table <-round(median_table,digits = 2)
  saveRDS(median_table,'median_probability_table_E16_to_P1.RDS')
  # Fates to P5 
  # Probability of each cell transitioning into a given cluster identity at subsequent time points
  fates_to_p5<- read.table('to_P5_fates.txt', header = TRUE)
  rownames(fates_to_p5) <- fates_to_p5$id
  fates_to_p5$cluster <- obj@meta.data[rownames(fates_to_p5),"final_clusters_renamed"]
  fates_to_p5$age <- obj@meta.data[rownames(fates_to_p5),"time_point"]
  plot_to_p5 <- reshape2::melt(fates_to_p5, value.name = "Probability")
  colnames(plot_to_p5) <- c('id','cluster',"age","Fate","Probability")
  # we extract P1 to P5 probabilities
  plot_to_p5_p1<- subset(plot_to_p5, subset=age=="P1")
  # for each cluster we calculate the median fate probabilities toward all possible cluster identities at the following time point
  by_cluster <- split(plot_to_p5_p1,plot_to_p5_p1$cluster)
  median_table <- Reduce("rbind",lapply(by_cluster, function(data){
  by_fate<-split(data,data$Fate)
  median<-unlist(lapply(by_fate, function(data_sub){
      median(data_sub$Probability)
  }))
  return(median)
  }))
  rownames(median_table) <- names(by_cluster)
  median_table<-median_table[,rownames(median_table)]
  median_table <-round(median_table,digits = 2)
  saveRDS(median_table,'median_probability_table_P1_to_P5.RDS')
}
# Merge subsamplings 
# E16 to P1
median_1 <-readRDS('cell_subset_1/median_probability_table_E16_to_P1.RDS')
median_2 <-readRDS('cell_subset_2/median_probability_table_E16_to_P1.RDS')
median_3 <-readRDS('cell_subset_3/median_probability_table_E16_to_P1.RDS')
# Compute mean across subsamplings
clusters<-rownames(median_1)
median<-matrix(ncol = length(clusters), nrow = length(clusters))
rownames(median) <- colnames(median) <- clusters
for (i in clusters){
  for (j in clusters){
    median[i,j] <- mean(c(median_1[i,j], median_2[i,j], median_3[i,j]))
  }
}
saveRDS(median, 'MedianProbability_table_NonMartinotti_finalClusters_E16_to_P1.RDS')
median_table <-round(median,digits = 2)
min_val <- min(median_table)
max_val <- max(median_table)
breaks <-  seq(min_val, max_val, length.out=1001);
cols <- colorRampPalette(colors = c('black','#2D718EFF','#3CBC75FF','gold'))(1000)
pdf(file="Median_probability_NonMartinotti_final_clusters_E16_to_P1.pdf",width=12,height = 12)
gplots::heatmap.2(median_table,
                  margins=c(8,8),
                  key = FALSE, # Explicitly set key to TRUE
                  keysize=1,
                  key.xlab="",
                  key.title="Probability (median) - from E16 to P5",
                  trace = "none",
                  density.info = "none",
                  col = cols,
                  breaks = breaks,
                  offsetRow=0.1,
                  offsetCol=0.1,
                  cexRow = 0.8,
                  cexCol = 0.8,
                  cellnote = median_table,
                  notecex = 0.9,
                  notecol = 'white',
                  Colv = F,
                  Rowv = F,
                  dendrogram = "none")
dev.off()
# P1 to P5
median_1 <-readRDS('cell_subset_1/median_probability_table_P1_to_P5.RDS')
median_2 <-readRDS('cell_subset_2/median_probability_table_P1_to_P5.RDS')
median_3 <-readRDS('cell_subset_3/median_probability_table_P1_to_P5.RDS')
# Compute mean across subsamplings
clusters<-rownames(median_1)
median<-matrix(ncol = length(clusters), nrow = length(clusters))
rownames(median) <- colnames(median) <- clusters
for (i in clusters){
  for (j in clusters){
    median[i,j] <- mean(c(median_1[i,j], median_2[i,j], median_3[i,j]))
  }
}
saveRDS(median, 'MedianProbability_table_NonMartinotti_finalClusters_P1_to_P5.RDS')
median_table <-round(median,digits = 2)
min_val <- min(median_table)
max_val <- max(median_table)
breaks <-  seq(min_val, max_val, length.out=1001);
cols <- colorRampPalette(colors = c('black','#2D718EFF','#3CBC75FF','gold'))(1000)
pdf(file="Median_probability_NonMartinotti_final_clusters_P1_to_P5.pdf",width=12,height = 12)
gplots::heatmap.2(median_table,
                  margins=c(8,8),
                  key = FALSE, # Explicitly set key to TRUE
                  keysize=1,
                  key.xlab="",
                  key.title="Probability (median) - from E16 to P5",
                  trace = "none",
                  density.info = "none",
                  col = cols,
                  breaks = breaks,
                  offsetRow=0.1,
                  offsetCol=0.1,
                  cexRow = 0.8,
                  cexCol = 0.8,
                  cellnote = median_table,
                  notecex = 0.9,
                  notecol = 'white',
                  Colv = F,
                  Rowv = F,
                  dendrogram = "none")
dev.off()
