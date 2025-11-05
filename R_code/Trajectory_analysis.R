# Trajectory analysis is performed on LRP neurons using destiny R package
# First we identify gene modules that underlie time-dependent transcriptomic states using Antler https://juliendelile.github.io/Antler/index.html

########################
# Load required packages
########################
library(Antler)
library(Seurat)
library(SeuratObject)
library(reshape2)
library(gprofiler2)
library(readr)
library(destiny)
library(scales)
library(cowplot)
library(tidyverse)
library(gtools)
library(ggthemes)
library(gam)
library(Hmisc)
library(circlize)
library(ComplexHeatmap)
'%!in%' <- function(x,y)!('%in%'(x,y))

#####################
# Antler gene modules 
#####################
# Prepare files for Antler
antler_path<-getwd()
atlas_v2_LRP<-subset(atlas_v2, subset=major_label_transferAnchors_BaseAtlas == "LRP")
atlas_v2_LRP<-subset(atlas_v2_LRP, subset=final_clusters_renamed %in% c("LRP2.1","LRP1.0"))  
atlas_v2_LRP<-subset(atlas_v2_LRP, subset=time_point %in% c("E16","P1","P5"))    
atlas_v2_LRP$time <- parse_number(atlas_v2_LRP$time_point)
atlas_v2_LRP$time[which(atlas_v2_LRP$time == 16)] <- -4
atlas_v2_LRP$replicate_id<-1
atlas_v2_LRP$replicate_id[which(atlas_v2_LRP$batch == "BaseAtlas_lim_P5_fixed_sorted")] <- 2
atlas_v2_LRP$replicate_id[which(atlas_v2_LRP$batch == "E16_Lim4")] <- 3
atlas_v2_LRP$replicate_id[which(atlas_v2_LRP$batch == "E16_EMI018")] <- 4
atlas_v2_LRP$replicate_id[which(atlas_v2_LRP$batch == "E18_Lippi")] <- 5
atlas_v2_LRP$replicate_id[which(atlas_v2_LRP$batch == "P1_EMI014_Lim")] <- 6
atlas_v2_LRP$replicate_id[which(atlas_v2_LRP$batch == "P1_Lim3")] <- 7
atlas_v2_LRP$replicate_id[which(atlas_v2_LRP$batch == "P1_Lim5")] <- 8
atlas_v2_LRP$replicate_id[which(atlas_v2_LRP$batch %in%  c("P5_WT1_Lim1","P5_WT23_Lim2"))] <- 9
atlas_v2_LRP$treatment<-'None'
# Run Antler and extract gene modules
antler <- Antler$new(output_folder=antler_path)
antler$load_dataset(assayData = as.matrix(atlas_v2_LRP@assays[['RNA']]@counts), phenoData = AnnotatedDataFrame(as.data.frame(atlas_v2_LRP@meta.data)) )
antler$remove_outliers(
  min_genes = 700,
  min_cells = 10)
antler$normalize('CPM')
antler$gene_modules$identify(
  name                  = "unbiasedGMs",
  corr_t                = 0.3,  # the Spearman correlation treshold
  corr_min              = 5,    # min. number of genes a gene must correlate with
  mod_min_cell          = 10,   # min. number of cells expressing the module
  mod_consistency_thres = 0.3,  # ratio of expressed genes among "positive" cells
  process_plots         = TRUE, # plot optimal module number heuristics
  num_cores             = 5,   # number of cores to use
  display=FALSE
)
modList<-antler$gene_modules$get("unbiasedGMs")
names(modList)<-paste("Mod",as.character(1:length(antler$gene_modules$get("unbiasedGMs"))))
file<-melt(modList)
colnames(file)<-c('Gene.Symbol','Module')
# Save to directory for later reference
saveRDS(modList,paste0(antler_path,'/Antler_gene_modules.RDS'))
# Generate functional enrichment for each gene module using gprofileR
modList<-readRDS(paste0(antler_path,'/Antler_gene_modules.RDS'))
linkList<-c()
GOenrichment<-list()
for(mod in 1:length(modList)){
  g<-modList[[mod]]
  #if there are enough genes, try gprofiler
  if(length(g)>=10){
    #warning("Attempting to generate functional enrichment using gprofileR. Sometimes this fails due to bad connection / problems with the gprofiler online service, and may simply need to be rerun.")
    res<-gost(g,organism='mmusculus',as_short_link=F)$result
    res <- data.frame(lapply(res, as.character), stringsAsFactors=FALSE)
    GOenrichment[[paste0('Mod ',mod)]] <- res
    write.csv(res,file=paste0(antler_path,'/Module_',as.character(mod),'_Functional_Enrichment.csv'))
    linkList<-rbind( linkList,  c(paste0('Mod',mod), gost(g,organism='mmusculus',as_short_link=T)))
  }
}
write.csv(linkList,file=paste0(antler_path,'/gprofiler_links.csv'))
saveRDS(GOenrichment, paste0(antler_path,'/Antler_gene_modules_GOenrichment.RDS'))
# Select gene modules associated with terms including “neuron differentiation”, “axon development”, and “synapse assembly”.
modules_sel<-c()
for(mod in 1:length(modList)){
  if(sum(GOenrichment[[mod]]$term_name %in% c("neuron differentiation","synapse formation","axon guidance","axon growth")) > 0) modules_sel<-c(modules_sel,names(GOenrichment)[mod])
}
# Extract genes from selected gene modules
keygenes<-c()
for(i in modules_sel){
  keygenes<-unique(c(keygenes, modList[[i]]))
}

######################
# Trajectory analysis
######################
# We run trajectory analysis using gene modules extracted with Antler
# extract the corrected gene expression for the genes in the selected modules
dat<-atlas_v2_LRP@assays$integrated@data
reducedDat<-dat[keygenes[keygenes %in% rownames(dat)],]
#run diffmap
if(dim(reducedDat)[1]<100){
  dm <- destiny::DiffusionMap(t(as.matrix(reducedDat)))
  
}else{
  dm <- destiny::DiffusionMap(t(as.matrix(reducedDat)), n_pcs = 20)
}
dpt<-destiny::DPT(dm)
# Select root cell from E16 time point cells
# root cell is identified as the cell from E16 samples that has the highest number of early sample cells in its neighborhood (computed on diffusion components). 
e16_cells<-rownames(atlas_v2_LRP@meta.data[which(atlas_v2_LRP@meta.data$time_point == "E16"),])
#get k nearest neighbours by eigenvalue distances
ev<-dm@eigenvectors
distance<-dist(ev)
distance<-as.matrix(distance)
k = 50
score<-rep(NA,length(e16_cells))
names(score)<-e16_cells
#get k nearest neighbours for each tip excluded cells from the same sample, and sum all the time values. smallest score = earliest tip. 
for(cell in e16_cells){
  dist<-distance[cell,]
  neighbours<-names(sort(dist)[1:k])
  score[cell]<-sum(atlas_v2_LRP$time[neighbours])
}
fetchval<-which(score == min(score))
rootcell<-names(score)[fetchval]
rootval<-grep(rootcell, rownames(dm@eigenvectors))
# assign root cell
sel_cells<-rownames(dpt@branch[which(dpt@branch[,1] == 3),])
tips[which(tips %in% sel_cells)] <- rootcell
# Recompute diffusion pseudotime with new root cell
dpt2<-destiny::DPT(dm,tips = which(rownames(dm@eigenvectors) %in% tips))
# extract pseudotime values
pseudotime<-dpt2[[paste0('DPT',rootval)]]
names(pseudotime)<-rownames(dm@eigenvectors)
# save pseudotime and diffusion components into Seurat object
atlas_v2_LRP$pseudotime <- pseudotime[colnames(atlas_v2_LRP)]
atlas_v2_LRP[["DC"]] <- CreateDimReducObject(embeddings = eigenvectors(dm)[colnames(obj),1:2], key = "DC_", assay = DefaultAssay(atlas_v2_LRP))
##########################
# Pseudotime related genes
##########################
DefaultAssay(atlas_v2_LRP) <- 'RNA'
atlas_v2_LRP <- ScaleData(atlas_v2_LRP,vars.to.regress = c("nFeature_RNA",'percent.mt','ccDiff'), verbose = TRUE)
# define early and late cells using pseudotime values
early_cells<-names(pseudotime[which(pseudotime<0.1)])
late_cells<-names(pseudotime[which(pseudotime>0.2)])
atlas_v2_LRP$groups_pseudotime <- 0
atlas_v2_LRP$groups_pseudotime[early_cells] <- 1
atlas_v2_LRP$groups_pseudotime[late_cells] <- 2
# identfy genes that are more expressed in early and late cells
Idents(atlas_v2_LRP) <- 'groups_pseudotime'
early<-FindMarkers(atlas_v2_LRP, ident.1 = "1", only.pos = TRUE, logfc.threshold = log2(1.5))
early<-early[which(early$p_val_adj<=0.01),]
late<-FindMarkers(atlas_v2_LRP, ident.1 = "2", only.pos = TRUE, logfc.threshold = log2(1.5))
late<-late[which(late$p_val_adj<=0.01),]
# compute correlation of genes with pseudotime and run GAM 
t <- pseudotime
t <- t[order(t)] 
y<- atlas_v2_LRP@assays$RNA@scale.data[,names(t)]
y <- y[which(rowSums(y==0)!=ncol(y)),]
gam.pval <- apply(y,1,function(z){
  d <- data.frame(z=z, t=t)
  tmp <- mgcv::gam(z ~ lo(t), data=d)
  p <- summary(tmp)
  return(p$p.p[2])
}) # estimate gam pvalue
corr <- apply(y,1,function(z){
  cor(t,z, method = "pearson")
}) # compute correlation with pseudotime
p_adj<-p.adjust(gam.pval)
# select early and late genes that correlate wuth pseudotime
genes_pseudotime <- data.frame(gam.p=gam.pval,gam.padj=p_adj, corr=corr)
genes_pseudotime <- genes_pseudotime[which(genes_pseudotime$gam.padj < 0.01),]
genes_pseudotime <- genes_pseudotime[order(genes_pseudotime$gam.padj),]
genes_pseudotime_early<-genes_pseudotime[rownames(early)[rownames(early) %in% rownames(genes_pseudotime)],]
genes_pseudotime_early<-genes_pseudotime_early[which(genes_pseudotime_early$corr < -0.15),]
genes_pseudotime_early <- genes_pseudotime_early[order(genes_pseudotime_early$corr, decreasing = FALSE),]
genes_pseudotime_late<-genes_pseudotime[rownames(late)[rownames(late) %in% rownames(genes_pseudotime)],]
genes_pseudotime_late<-genes_pseudotime_late[which(genes_pseudotime_late$corr > 0.2),]
genes_pseudotime_late <- genes_pseudotime_late[order(genes_pseudotime_late$corr, decreasing = TRUE),]
# Replicate same analysis on clusters
# Identify genes that are differentially expressed between LRP.1 and LRP.2
Idents(atlas_v2_LRP) <- 'final_clusters_renamed'
cl2<-FindMarkers(atlas_v2_LRP, ident.1 = "LRP2.1", only.pos = TRUE, logfc.threshold = log2(1.5))
cl2<-cl2[which(cl2$p_val_adj<=0.01),]
cl1<-FindMarkers(atlas_v2_LRP, ident.1 = "LRP1.0", only.pos = TRUE, logfc.threshold = log2(1.5))
cl1<-cl1[which(cl1$p_val_adj<=0.01),]
# genes_pseudotime_late<-genes_pseudotime_late[rownames(genes_pseudotime_late) %!in% c(rownames(cl1),rownames(cl2)),]
# gam and correlation on branch 2 and LRP2 cells
cells<-unique(c(rownames(obj@meta.data[which(obj@meta.data$final_clusters_renamed == "LRP2.1"),]),rownames(obj@meta.data[which(obj@meta.data$branch == 2),])))
t <- pseudotime[cells]
t <- t[order(t)] 
y<- obj@assays$RNA@scale.data[,names(t)]
y <- y[which(rowSums(y==0)!=ncol(y)),]
gam.pval <- apply(y,1,function(z){
  d <- data.frame(z=z, t=t)
  tmp <- mgcv::gam(z ~ lo(t), data=d)
  p <- summary(tmp)
  return(p$p.p[2])
})
corr <- apply(y,1,function(z){
  cor(t,z, method = "pearson")
}) # compute correlation with pseudotime
p_adj<-p.adjust(gam.pval)
sum(p_adj< 0.01)
genes_pseudotime <- data.frame(gam.p=gam.pval,gam.padj=p_adj, corr=corr)
genes_pseudotime <- genes_pseudotime[which(genes_pseudotime$gam.padj < 0.01),]
genes_pseudotime <- genes_pseudotime[order(genes_pseudotime$gam.padj),]
### for cluster2
genes_pseudotime2<-genes_pseudotime[rownames(cl2)[rownames(cl2) %in% rownames(genes_pseudotime)],]
# select genes that are more expressed in cluster LRP2 and whose expression correlates with pseudotime
genes_pseudotime2<-genes_pseudotime2[which(genes_pseudotime2$corr >0.2),]
ordered_LRP2<-names(t)
genes_pseudotime_cl2 <- genes_pseudotime2[order(genes_pseudotime2$corr, decreasing = TRUE),]
# gam and correlation on branch 1 and LRP1 cells
cells<-unique(c(rownames(obj@meta.data[which(obj@meta.data$final_clusters_renamed == "LRP1.0"),]),rownames(obj@meta.data[which(obj@meta.data$branch == 1),])))
t <- pseudotime[cells]
t <- t[order(t)] 
y<- obj@assays$RNA@scale.data[,names(t)]
y <- y[which(rowSums(y==0)!=ncol(y)),]
gam.pval <- apply(y,1,function(z){
  d <- data.frame(z=z, t=t)
  tmp <- mgcv::gam(z ~ lo(t), data=d)
  p <- summary(tmp)
  return(p$p.p[2])
})
corr <- apply(y,1,function(z){
  cor(t,z, method = "pearson")
}) # compute correlation with pseudotime
p_adj<-p.adjust(gam.pval)
genes_pseudotime <- data.frame(gam.p=gam.pval,gam.padj=p_adj, corr=corr)
genes_pseudotime <- genes_pseudotime[which(genes_pseudotime$gam.padj < 0.01),]
genes_pseudotime <- genes_pseudotime[order(genes_pseudotime$gam.padj),]
#for cluster1
genes_pseudotime1<-genes_pseudotime[rownames(cl1)[rownames(cl1) %in% rownames(genes_pseudotime)],]
# select genes that are more expressed in cluster LRP1 and whose expression correlates with pseudotime
genes_pseudotime1<-genes_pseudotime1[which(genes_pseudotime1$corr >0.15),]
ordered_LRP1<-names(t)
genes_pseudotime_cl1 <- genes_pseudotime1[order(genes_pseudotime1$corr, decreasing = TRUE),]
#saveRDS(genes_pseudotime_cl1,"genes_pseudotime_cl1.RDS")
# select early common genes
genes_pseudotime_early<-genes_pseudotime_early[rownames(genes_pseudotime_early) %!in% c(rownames(cl1),rownames(cl2)),]
#Plot heatmap
genes_pseudotime<-rbind(genes_pseudotime_cl1,genes_pseudotime_cl2,genes_pseudotime_early)
ordered_cells<-c(ordered_LRP1,ordered_LRP2)
data.use2 <- MinMax(data = obj@assays$RNA@scale.data[rownames(genes_pseudotime),ordered_LRP2], min = -3, max = 3)
data.use1 <- MinMax(data = obj@assays$RNA@scale.data[rownames(genes_pseudotime),ordered_LRP1], min = -3, max = 3)
gaps = factor(c(rep("LRP1.0", nrow(genes_pseudotime_cl1)),rep("LRP2.1",nrow(genes_pseudotime_cl2)),rep("early",nrow(genes_pseudotime_early))), levels=c("LRP1.0","LRP2.1","both","early"))
#col_fun = colorRamp2(c(0,0.18,0.36,0.54,0.72,0.9), alpha(c("darkorchid4","#554CA6","#4D89A6", "#7CC494",  "gold", "#F64444"),0.8))
col_fun = colorRamp2(c(0,0.07,0.13,0.24,0.4,0.6), alpha(c("darkorchid4","#554CA6","#4D89A6", "#7CC494",  "gold", "#F64444"),0.8))
pdf('Heatmap_pseudotime_genes.pdf', height = 10, width = 8)
Heatmap(as.matrix(data.use1), cluster_rows=F, cluster_columns=F, col = colorRamp2(c(-3,0,3), c("blue4","black","gold")), split=gaps, gap = unit(5, "mm"), show_row_names=T,row_names_side = "left",show_column_names=F,use_raster=F, name="Expression",
        top_annotation = HeatmapAnnotation(Pseudotime=anno_simple(pseudotime[ordered_LRP1], col=col_fun),show_annotation_name = F))  +
  Heatmap(as.matrix(data.use2), cluster_rows=F, cluster_columns=F, col = colorRamp2(c(-3,0,3), c("blue4","black","gold")), split=gaps, gap = unit(5, "mm"), show_row_names=T,row_names_side = "left",show_column_names=F,use_raster=F, name="Expression",
          top_annotation = HeatmapAnnotation(Pseudotime=anno_simple(pseudotime[ordered_LRP2], col=col_fun)))  
dev.off()  


