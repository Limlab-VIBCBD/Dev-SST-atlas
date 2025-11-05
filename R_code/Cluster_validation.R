# This is the code we used to validate cell identities, and in particular, the ability of a gene list to predict cluster identities.
# Five-fold cross-validation is performed, taking 20% of the data as a test set and the remaining 80% as the training set. 
# We train a random forest classifier with ntree trees on each pair of clusters in the training data, this classifier is then used to predict the identity of cells in the test set.
# The procedure is repeated for a total of nIteration iterations

# function to apply Random Forest validation and score the class membership. 
# A matrix is returned in which each row represents a cluster and each column indicates the fraction of cells predicted to belong to each cluster identity.
ComputeClassificationAccuracy <- function(seuratObject,
                            clustering,
                            minClusterRep=5,
                            ntree=100,
                            nIteration=10,
                            nCore=10,
                            ...){
  
  cells<-rownames(seuratObject[[clustering]])
  seuratObject <- subset(seuratObject, cells = cells[which(!is.na(seuratObject[[clustering]]))] )
  # take data from the RNA assay not the integrated assay
  DefaultAssay(seuratObject)<-'RNA'
  Input<-as.matrix(seuratObject@assays[['RNA']]@data)
  cate <- as.data.frame(as.factor(seuratObject[[clustering]][colnames(Input),]))
  rownames(cate)<-colnames(Input)
  feature <- t(as.matrix(Input))
  #size of test set
  ntest<-floor(dim(feature)[1]/5)
  #remainder (add on to the last group)
  rem<-dim(feature)[1]%%5
  cell_names<-rownames(feature)
  #Record possible combinations of classes
  classes<-unique(seuratObject[[clustering]])
  combos<-combn(classes[,1],m=2)
  #store predictions for each cell
  #predictions<-matrix(nrow=length(cell_names),ncol=10*dim(combos)[2])
  registerDoMC(nCore)

  out<-lapply(c(1:nIteration), function(iteration){
    #for reproducibility
    set.seed(iteration)
    #matrix holding individual cell predictions for this iteration
    predictions<-matrix(nrow=length(cell_names),ncol=dim(combos)[2])
    rownames(predictions)<-cell_names
    colnames(predictions)<-1:dim(combos)[2]
    print(paste0("starting iteration: ",as.character(iteration)))
    test_sets<-list()
    waitCounter<-0
    successfulTestSets<-FALSE
    while(successfulTestSets==FALSE){
      possibilities<-1:dim(feature)[1]
      #Generate 5 test sets, that are each a different fifth of the data
      for(i in 1:5){
        if(i==5){
          #not divisible by 5, so make up the remainder in the last round
          set<-base::sample(possibilities,size=ntest+rem,replace=FALSE)
        }else{
          set<-base::sample(possibilities,size=ntest,replace=FALSE)
        }
        test_sets[[i]]<-set        
        possibilities<-possibilities[is.na(match(possibilities,set))]
      }
      #Check: Does each test set contain at least 5 samples from each class in our clustering? 
      keeptry=TRUE
      for(i in 1:5){
        if(keeptry){
          
          counts<-table(seuratObject[[clustering]][test_sets[[i]],])
          
          if(sum(counts<=minClusterRep) > 0 ){
            test_sets<-list()
            waitCounter<-waitCounter+1
            keeptry=FALSE;
          }
          else if(waitCounter>10){
            print("class proportions too low! Generating test sets failed")
            break()
          }
        }
        
      }
      
      if(length(test_sets)==5){
        successfulTestSets=TRUE
      }
    }
    
    for(setnum in 1:5){
      test<-test_sets[[setnum]]
      test_feat<-feature[test,]
      test_cate<-cate[test,]

      training_feat<-feature[-test,]
      training_cate<-cate[-test,]
      
      preds<-list()
      for(i in 1:dim(combos)[2]){
        #print(i)
        comb<-combos[,i]
        ind<-c(which(!is.na(match(training_cate,comb[1]))),which(!is.na(match(training_cate,comb[2]))))
        training_feat_sub<-training_feat[ind,]
        training_cate_sub<-droplevels(training_cate[ind])
        #for subsampling
        minSize<-min(table(training_cate_sub))
        randf<-randomForest::randomForest(x=training_feat_sub, y=training_cate_sub, ntree=ntree,
                                          sampsize=rep(minSize,2)
                                          #classwt=table(training_cate_sub)
        ) 
        #Classify test data
        preds[[i]]<-predict(randf,test_feat,importance=T)
        #update objects
        #rfs[[length(rfs)+1]]<-list(randf,comb,setnum)
      }
      #save predictions to master table
      for(j in 1:dim(combos)[2]){
        dat<-preds[[j]]
        ind<-match(names(dat),rownames(predictions))
        predictions[ind,j]<-as.character(dat)
      }
    }
    return(predictions)
  })
  #End of 10 iterations. After this we have 10 predictions for each cell. 
  saveRDS(out,paste0("predictions_list_",clustering,".RDS"))
  
  #Merge all the predictions into one object
  predictionTable<-vector()
  for(i in 1:dim(combos)[2]){
    block<-vector()
    for(j in 1:length(out)){
      block<-cbind(block,out[[j]][,i])
    }
    predictionTable<-cbind(predictionTable,block)
  }
  #This has all the predictions together
  predictions<-predictionTable
  #This matrix will hold a class membership score for each cell, for each class
  membershipScore<-matrix(0,nrow=dim(predictions)[1],ncol=length(classes[,1]))
  rownames(membershipScore)<-rownames(predictions)
  colnames(membershipScore)<-sort(classes[,1])

  for(i in  1:dim(predictions)[1]){
    #grab the predictions just for this cell
    cell<-predictions[i,]
    #what are the possible class identities of the cell?
    class_set<-classes[,1]
    for(class in 1:dim(combos)[2]){
      #for each pairwise class comparison, we grab the predictions corresponding to that comparison. They are in blocks of nIteration
      chunk<-cell[((nIteration*class)-(nIteration-1)):(nIteration*class)]
      #record any class in this pairwise comparison which does not appear in the predictions. It has been dominated by the other possibility 
      dominated<-combos[which(is.na(match(combos[,class],chunk))),class]
      #remove any dominated classes from the set of possible identities
      if(length(dominated!=0)){
        class_set[match(dominated,class_set)]<-NA
      }
    }
    
    #consider the classes that remain
    class_set<-class_set[!is.na(class_set)]
    # evaluated only the predictions for the remaining classes
    cell<-cell[!is.na(match(cell,class_set))]
    #the membership score is the proportion of predicitons for each class
    scores<-table(cell)/length(cell)
    #record in matrix
    membershipScore[i,as.character(names(scores))]<-as.numeric(scores)
  }
  
  #We need to give a certainty category, based on the value of scores
  #We also need a final classification, taking highest valued score. 
  certainty<-vector(length=dim(membershipScore)[1])
  classification<-vector(length=dim(membershipScore)[1])
  for(i in 1:length(certainty)){
    #This holds all the scores for the cell. Most of them will be zero.
    cell<-membershipScore[i,]
    #seek the class with the highest membership score
    find<-grep(max(as.numeric(cell[1:length(classes[,1])])),as.numeric(cell[1:length(classes[,1])]))
    if(sum(as.numeric(as.character(cell[1:length(classes[,1])]))==1)==1){
      #If there is one class with a score of 1, we have identified a single dominant identity: this is a core cell
      certainty[i]<-'Core'
    }else if(sum(as.numeric(as.character(cell[1:length(classes[,1])])))==0  || length(find)>1){
      #if no identities have any positive score (all identities are dominated), or if the maximal score is found in more than one column (multiple classes have the same score) then we cannot discern identity.
      certainty[i]<-'Failure'
    }else{
      #The remaining cases are those with non-zero scores that are less than one, and one class has the highest score. These cells have multiple identity ties
      certainty[i]<-'Intermediate'
    }
    if(certainty[i]=='Failure'){
      classification[i]<-"U"
    }else{
      classification[i]<-find
    }
  }
  membershipScore<- as.data.frame(cbind(membershipScore,certainty,classification))
  out<-list(membershipScore,predictions)
  classification <- unlist(lapply(out[[1]]$classification, function(x){
    if(x != "U"){
      class<- colnames(out[[1]])[as.numeric(x)]
    }else{
      class <- "U"
    }
    return(class)
  }))
  out[[1]]$classification <- classification
  if(sum(out[[1]]$certainty == "Failure")/nrow(out[[1]]) > 0.1) {print("WARNING: Annotation failed for more than 10% of cells")}
  SM <- out[[1]][which(out[[1]]$certainty != "Failure"),]
  SM$classification <-factor(SM$classification, levels=unique(seuratObject@meta.data[,clustering]))
  accuracy_matrix <-vector()
  for (i in unique(seuratObject@meta.data[,clustering])){
    cells <-rownames(seuratObject@meta.data[which(seuratObject@meta.data[,clustering] %in% i),])
    cells <- cells[cells %in% rownames(SM)]
    accuracy_matrix <- rbind(accuracy_matrix,round(table(SM[cells,"classification"])/length(cells),3))
  }
  rownames(accuracy_matrix) <- unique(seuratObject@meta.data[,clustering])
  return(accuracy_matrix)
}


########################
# Load required packages
########################
library(Seurat)
library(SeuratWrappers)
library(SeuratDisk)
library(Seurat)
library(randomForest)
library(future)
library(doMC)
library(dplyr)
library(magrittr)
'%!in%' <- function(x,y)!('%in%'(x,y))
# set color palette 
breaks <-  seq(0, 0.5, length.out=101)
#cols <- colorRampPalette(colors = c('black','#EE9A3A'))(1000)
cols <- colorRampPalette(colors = c('black','#3CBC75FF','#F0B800'))(100)


#####################
# LRP subset atlas v2
#####################
atlas_v2_LRP<-subset(atlas_v2, subset=major_label_transferAnchors_BaseAtlas == "LRP")
# compute marker genes 
DefaultAssay(atlas_v2_LRP) <-'RNA'
Idents(atlas_v2_LRP) <- 'final_clusters_renamed'
markers <- FindAllMarkers(atlas_v2_LRP,logfc.threshold = log2(1.5), min.pct = 0.2,only.pos = TRUE, test.use = 'MAST', latent.vars = 'batch')
ribo_genes <- grep(pattern = "^Rp[sl]", x = rownames(atlas_v2_LRP@assays$RNA@counts), value = TRUE)
markers<-markers[which(markers$gene %!in% ribo_genes),]
mt_genes <- grep(pattern = "^mt-", x = rownames(atlas_v2_LRP@assays$RNA@counts), value = TRUE)
markers<-markers[which(markers$gene %!in% mt_genes),]
markers<-markers[order(markers$avg_log2FC, decreasing = TRUE),]
markers<-markers[which(markers$p_val_adj < 0.01),]
markers$logFC<-markers$avg_log2FC
# select top 10 markers (ranked by logFC) for each cluster 
markers<-markers %>% group_by(cluster) %>% top_n(n = 10)
markers<-as.data.frame(markers)
gene_list<-unique(markers$gene)
# subset count matrix (RNA assay) on the marker gene list
LRP_sub<-CreateSeuratObject(atlas_v2_LRP@assays$RNA@counts[gene_list,])
LRP_sub<-AddMetaData(LRP_sub, atlas_v2_LRP@meta.data[colnames(LRP_sub),])
# Compute the accuracy matrix on the subsetted count matrix to assess the ability of marker genes to predict cluster identities.
accuracy_matrix_v2<-ComputeClassificationAccuracy(LRP_sub, "final_clusters_renamed")
saveRDS(accuracy_matrix_v2, "accuracy_matrix_LRP_atlas_v2.RDS")
# plot matrix
accuracy_matrix_v2<-accuracy_matrix_v2[colnames(accuracy_matrix_v2)[order(colnames(accuracy_matrix_v2))],colnames(accuracy_matrix_v2)[order(colnames(accuracy_matrix_v2))]]
pdf('RF_accuracy_matrix_LRPs_atlas_v2.pdf')
gplots::heatmap.2(accuracy_matrix_v2, margins=c(7,7), key = FALSE, keysize=1, key.xlab="", key.title="Accuracy", trace = "none", 
                density.info = "none", col = cols, breaks = breaks, offsetRow=0.1, offsetCol=0.1, cexRow = 0.8, cexCol = 0.8,
                cellnote = accuracy_matrix_v2, notecex = 0.8, notecol = 'white', Colv = F, Rowv = F, dendrogram = "none")
dev.off()

#####################
# LRP subset atlas v1
#####################
atlas_v1_LRP<-subset(atlas_v1, subset=major_label_transferAnchors_BaseAtlas == "LRP")
# compute marker genes 
DefaultAssay(atlas_v1_LRP) <-'RNA'
Idents(atlas_v1_LRP) <- 'final_clusters_renamed'
markers <- FindAllMarkers(atlas_v1_LRP,logfc.threshold = log2(1.5), min.pct = 0.2,only.pos = TRUE, test.use = 'MAST', latent.vars = 'batch')
ribo_genes <- grep(pattern = "^Rp[sl]", x = rownames(atlas_v1_LRP@assays$RNA@counts), value = TRUE)
markers<-markers[which(markers$gene %!in% ribo_genes),]
mt_genes <- grep(pattern = "^mt-", x = rownames(atlas_v1_LRP@assays$RNA@counts), value = TRUE)
markers<-markers[which(markers$gene %!in% mt_genes),]
markers<-markers[order(markers$avg_log2FC, decreasing = TRUE),]
markers<-markers[which(markers$p_val_adj < 0.01),]
markers$logFC<-markers$avg_log2FC
# select top 10 markers (ranked by logFC) for each cluster 
markers<-markers %>% group_by(cluster) %>% top_n(n = 10)
markers<-as.data.frame(markers)
gene_list<-unique(markers$gene)
# subset count matrix (RNA assay) on the marker gene list
LRP_sub<-CreateSeuratObject(atlas_v1_LRP@assays$RNA@counts[gene_list,])
LRP_sub<-AddMetaData(LRP_sub, atlas_v1_LRP@meta.data[colnames(LRP_sub),])
# Compute the accuracy matrix on the subsetted count matrix to assess the ability of marker genes to predict cluster identities.
accuracy_matrix_v1<-ComputeClassificationAccuracy(LRP_sub, "final_clusters_renamed")
saveRDS(accuracy_matrix_v1, "accuracy_matrix_LRP_atlas_v1.RDS")
# plot matrix
accuracy_matrix_v1<-accuracy_matrix_v1[colnames(accuracy_matrix_v1)[order(colnames(accuracy_matrix_v1))],colnames(accuracy_matrix_v1)[order(colnames(accuracy_matrix_v1))]]
pdf('RF_accuracy_matrix_LRPs_atlas_v1.pdf')
gplots::heatmap.2(accuracy_matrix_v1, margins=c(7,7), key = FALSE, keysize=1, key.xlab="", key.title="Accuracy", trace = "none", 
                density.info = "none", col = cols, breaks = breaks, offsetRow=0.1, offsetCol=0.1, cexRow = 0.8, cexCol = 0.8,
                cellnote = accuracy_matrix_v1, notecex = 0.8, notecol = 'white', Colv = F, Rowv = F, dendrogram = "none")
dev.off()

# cluster matching - atlas v1 vs atlas v2
matches <- list("HIn1"=c("HIn1"), 
                "HIn2"=c("HIn2"), 
                "LRP1.0"=c("LRP1.0"), 
                "LRP2.1"=c("LRP2.1"), 
                "LRP2.2"=c("LRP2.2")
)
# compute delta on diagonal of accuracy matrices (atlas v2 - atlas v1)
delta_diag <- unlist(lapply(names(matches), function(cl){
  diag_2<-accuracy_matrix_v2[cl,cl]
  if(length(matches[[cl]]) > 1){
    #diag_1<-mean(diag(accuracy_matrix_v1[matches[[cl]],matches[[cl]]]))
    diag_1<-diag(accuracy_matrix_v1[matches[[cl]],matches[[cl]]])
  }else{
    diag_1<-accuracy_matrix_v1[matches[[cl]],matches[[cl]]]
  }
  delta_diag <- diag_2 - diag_1
  return(delta_diag)
}))
# compute delta on cross-talk of accuracy matrices (atlas v2 - atlas v1)
delta_crosstalk <- unlist(lapply(names(matches), function(cl){
  crosstalk_2<-sum(accuracy_matrix_v2[cl,!(colnames(accuracy_matrix_v2) %in% cl)])
  if(length(matches[[cl]]) > 1){
    crosstalk_1 <- unlist(lapply(matches[[cl]],function(x){
      crosstalk_1<-sum(accuracy_matrix_v1[x,!(colnames(accuracy_matrix_v1) %in% x)])
      return(crosstalk_1)
    }))
  }else{
    crosstalk_1<-sum(accuracy_matrix_v1[matches[[cl]],!(colnames(accuracy_matrix_v1) %in% matches[[cl]])])
  }
  delta_crosstalk <- crosstalk_2 - crosstalk_1
  return(delta_crosstalk)
}))
atlas_delta_summary_LRP <- data.frame(Family="LRPs",ClassIndex=c(0:(length(delta_crosstalk)-1)),diag_delta=delta_diag,ct_delta=delta_crosstalk)


############################
# Martinotti subset atlas v2
############################
atlas_v2_Martinotti<-subset(atlas_v2, subset=major_label_transferAnchors_BaseAtlas == "Martinotti")
# compute marker genes 
DefaultAssay(atlas_v2_Martinotti) <-'RNA'
Idents(atlas_v2_Martinotti) <- 'final_clusters_renamed'
markers <- FindAllMarkers(atlas_v2_Martinotti,logfc.threshold = log2(1.5), min.pct = 0.2,only.pos = TRUE, test.use = 'MAST', latent.vars = 'batch')
ribo_genes <- grep(pattern = "^Rp[sl]", x = rownames(atlas_v2_Martinotti@assays$RNA@counts), value = TRUE)
markers<-markers[which(markers$gene %!in% ribo_genes),]
mt_genes <- grep(pattern = "^mt-", x = rownames(atlas_v2_Martinotti@assays$RNA@counts), value = TRUE)
markers<-markers[which(markers$gene %!in% mt_genes),]
markers<-markers[order(markers$avg_log2FC, decreasing = TRUE),]
markers<-markers[which(markers$p_val_adj < 0.01),]
markers$logFC<-markers$avg_log2FC
# select top 10 markers (ranked by logFC) for each cluster 
markers<-markers %>% group_by(cluster) %>% top_n(n = 10)
markers<-as.data.frame(markers)
gene_list<-unique(markers$gene)
# subset count matrix (RNA assay) on the marker gene list
Martinotti_sub<-CreateSeuratObject(atlas_v2_Martinotti@assays$RNA@counts[gene_list,])
Martinotti_sub<-AddMetaData(Martinotti_sub, atlas_v2_Martinotti@meta.data[colnames(Martinotti_sub),])
# Compute the accuracy matrix on the subsetted count matrix to assess the ability of marker genes to predict cluster identities.
accuracy_matrix_v2<-ComputeClassificationAccuracy(Martinotti_sub, "final_clusters_renamed")
saveRDS(accuracy_matrix_v2, "accuracy_matrix_Martinotti_atlas_v2.RDS")
# plot matrix
accuracy_matrix_v2<-accuracy_matrix_v2[colnames(accuracy_matrix_v2)[order(colnames(accuracy_matrix_v2))],colnames(accuracy_matrix_v2)[order(colnames(accuracy_matrix_v2))]]
pdf('RF_accuracy_matrix_Martinotti_atlas_v2.pdf')
gplots::heatmap.2(accuracy_matrix_v2, margins=c(7,7), key = FALSE, keysize=1, key.xlab="", key.title="Accuracy", trace = "none", 
                density.info = "none", col = cols, breaks = breaks, offsetRow=0.1, offsetCol=0.1, cexRow = 0.8, cexCol = 0.8,
                cellnote = accuracy_matrix_v2, notecex = 0.8, notecol = 'white', Colv = F, Rowv = F, dendrogram = "none")
dev.off()

############################
# Martinotti subset atlas v1
############################
atlas_v1_Martinotti<-subset(atlas_v1, subset=major_label_transferAnchors_BaseAtlas == "Martinotti")
# compute marker genes 
DefaultAssay(atlas_v1_Martinotti) <-'RNA'
Idents(atlas_v1_Martinotti) <- 'final_clusters_renamed'
markers <- FindAllMarkers(atlas_v1_Martinotti,logfc.threshold = log2(1.5), min.pct = 0.2,only.pos = TRUE, test.use = 'MAST', latent.vars = 'batch')
ribo_genes <- grep(pattern = "^Rp[sl]", x = rownames(atlas_v1_Martinotti@assays$RNA@counts), value = TRUE)
markers<-markers[which(markers$gene %!in% ribo_genes),]
mt_genes <- grep(pattern = "^mt-", x = rownames(atlas_v1_Martinotti@assays$RNA@counts), value = TRUE)
markers<-markers[which(markers$gene %!in% mt_genes),]
markers<-markers[order(markers$avg_log2FC, decreasing = TRUE),]
markers<-markers[which(markers$p_val_adj < 0.01),]
markers$logFC<-markers$avg_log2FC
# select top 10 markers (ranked by logFC) for each cluster 
markers<-markers %>% group_by(cluster) %>% top_n(n = 10)
markers<-as.data.frame(markers)
gene_list<-unique(markers$gene)
# subset count matrix (RNA assay) on the marker gene list
Martinotti_sub<-CreateSeuratObject(atlas_v1_Martinotti@assays$RNA@counts[gene_list,])
Martinotti_sub<-AddMetaData(Martinotti_sub, atlas_v1_Martinotti@meta.data[colnames(Martinotti_sub),])
# Compute the accuracy matrix on the subsetted count matrix to assess the ability of marker genes to predict cluster identities.
accuracy_matrix_v1<-ComputeClassificationAccuracy(Martinotti_sub, "final_clusters_renamed")
saveRDS(accuracy_matrix_v1, "accuracy_matrix_Martinotti_atlas_v1.RDS")
# plot matrix
accuracy_matrix_v1<-accuracy_matrix_v1[colnames(accuracy_matrix_v1)[order(colnames(accuracy_matrix_v1))],colnames(accuracy_matrix_v1)[order(colnames(accuracy_matrix_v1))]]
pdf('RF_accuracy_matrix_Martinotti_atlas_v1.pdf')
gplots::heatmap.2(accuracy_matrix_v1, margins=c(7,7), key = FALSE, keysize=1, key.xlab="", key.title="Accuracy", trace = "none", 
                density.info = "none", col = cols, breaks = breaks, offsetRow=0.1, offsetCol=0.1, cexRow = 0.8, cexCol = 0.8,
                cellnote = accuracy_matrix_v1, notecex = 0.8, notecol = 'white', Colv = F, Rowv = F, dendrogram = "none")
dev.off()

# cluster matching - atlas v1 vs atlas v2
matches <- list("MC1.1"=c("MC1.1"), 
                 "MC1.10"=c("MC1.10"), 
                 "MC1.11"=c("MC1.4"), 
                 "MC1.2"=c("MC1.6"), 
                 "MC1.3"=c("MC1.2"), 
                 "MC1.4"=c("MC1.8","MC1.9"), 
                 "MC1.5"=c("MC1.3"), 
                 "MC1.6"=c("MC1.4"), 
                 "MC1.7"=c("MC1.7"), 
                 "MC1.8"=c("MC1.3"),
                 "MC1.9"=c("MC1.5"), 
                 "MC2.1"=c("MC2.1","MC2.2"), 
                 "MC2.2"=c("MC2.4"), 
                 "MC3.1"=c("MC3.1","MC3.2","MC2.3"), 
                 "MC3.2"=c("MC3.4"), 
                 "MC3.3"=c("MC3.3"), 
                 "MC3.4"=c("MC3.6"), 
                 "MC3.5"=c("MC3.5")
                 )
# compute delta on diagonal of accuracy matrices (atlas v2 - atlas v1)
delta_diag <- unlist(lapply(names(matches), function(cl){
    diag_2<-accuracy_matrix_v2[cl,cl]
    if(length(matches[[cl]]) > 1){
      #diag_1<-mean(diag(accuracy_matrix_v1[matches[[cl]],matches[[cl]]]))
      diag_1<-diag(accuracy_matrix_v1[matches[[cl]],matches[[cl]]])
    }else{
      diag_1<-accuracy_matrix_v1[matches[[cl]],matches[[cl]]]
    }
    delta_diag <- diag_2 - diag_1
    return(delta_diag)
}))
# compute delta on cross-talk of accuracy matrices (atlas v2 - atlas v1)
delta_crosstalk <- unlist(lapply(names(matches), function(cl){
  crosstalk_2<-sum(accuracy_matrix_v2[cl,!(colnames(accuracy_matrix_v2) %in% cl)])
  if(length(matches[[cl]]) > 1){
    crosstalk_1 <- unlist(lapply(matches[[cl]],function(x){
      crosstalk_1<-sum(accuracy_matrix_v1[x,!(colnames(accuracy_matrix_v1) %in% x)])
      return(crosstalk_1)
    }))
  }else{
    crosstalk_1<-sum(accuracy_matrix_v1[matches[[cl]],!(colnames(accuracy_matrix_v1) %in% matches[[cl]])])
  }
  delta_crosstalk <- crosstalk_2 - crosstalk_1
  return(delta_crosstalk)
}))
atlas_delta_summary_MC <- data.frame(Family="Martinotti",ClassIndex=c(0:(length(delta_crosstalk)-1)),diag_delta=delta_diag,ct_delta=delta_crosstalk)




################################
# Non-Martinotti subset atlas v2
################################
atlas_v2_NonMartinotti<-subset(atlas_v2, subset=major_label_transferAnchors_BaseAtlas == "Non-Martinotti")
# compute marker genes 
DefaultAssay(atlas_v2_NonMartinotti) <-'RNA'
Idents(atlas_v2_NonMartinotti) <- 'final_clusters_renamed'
markers <- FindAllMarkers(atlas_v2_NonMartinotti,logfc.threshold = log2(1.5), min.pct = 0.2,only.pos = TRUE, test.use = 'MAST', latent.vars = 'batch')
ribo_genes <- grep(pattern = "^Rp[sl]", x = rownames(atlas_v2_NonMartinotti@assays$RNA@counts), value = TRUE)
markers<-markers[which(markers$gene %!in% ribo_genes),]
mt_genes <- grep(pattern = "^mt-", x = rownames(atlas_v2_NonMartinotti@assays$RNA@counts), value = TRUE)
markers<-markers[which(markers$gene %!in% mt_genes),]
markers<-markers[order(markers$avg_log2FC, decreasing = TRUE),]
markers<-markers[which(markers$p_val_adj < 0.01),]
markers$logFC<-markers$avg_log2FC
# select top 10 markers (ranked by logFC) for each cluster 
markers<-markers %>% group_by(cluster) %>% top_n(n = 10)
markers<-as.data.frame(markers)
gene_list<-unique(markers$gene)
# subset count matrix (RNA assay) on the marker gene list
NonMartinotti_sub<-CreateSeuratObject(atlas_v2_NonMartinotti@assays$RNA@counts[gene_list,])
NonMartinotti_sub<-AddMetaData(NonMartinotti_sub, atlas_v2_NonMartinotti@meta.data[colnames(NonMartinotti_sub),])
# Compute the accuracy matrix on the subsetted count matrix to assess the ability of marker genes to predict cluster identities.
accuracy_matrix_v2<-ComputeClassificationAccuracy(NonMartinotti_sub, "final_clusters_renamed")
saveRDS(accuracy_matrix_v2, "accuracy_matrix_NonMartinotti_atlas_v2.RDS")
# plot matrix
accuracy_matrix_v2<-accuracy_matrix_v2[colnames(accuracy_matrix_v2)[order(colnames(accuracy_matrix_v2))],colnames(accuracy_matrix_v2)[order(colnames(accuracy_matrix_v2))]]
pdf('RF_accuracy_matrix_NonMartinotti_atlas_v2.pdf')
gplots::heatmap.2(accuracy_matrix_v2, margins=c(7,7), key = FALSE, keysize=1, key.xlab="", key.title="Accuracy", trace = "none", 
                density.info = "none", col = cols, breaks = breaks, offsetRow=0.1, offsetCol=0.1, cexRow = 0.8, cexCol = 0.8,
                cellnote = accuracy_matrix_v2, notecex = 0.8, notecol = 'white', Colv = F, Rowv = F, dendrogram = "none")
dev.off()

################################
# Non-Martinotti subset atlas v1
################################
atlas_v1_NonMartinotti<-subset(atlas_v1, subset=major_label_transferAnchors_BaseAtlas == "Non-Martinotti")
# compute marker genes 
DefaultAssay(atlas_v1_NonMartinotti) <-'RNA'
Idents(atlas_v1_NonMartinotti) <- 'final_clusters_renamed'
markers <- FindAllMarkers(atlas_v1_NonMartinotti,logfc.threshold = log2(1.5), min.pct = 0.2,only.pos = TRUE, test.use = 'MAST', latent.vars = 'batch')
ribo_genes <- grep(pattern = "^Rp[sl]", x = rownames(atlas_v1_NonMartinotti@assays$RNA@counts), value = TRUE)
markers<-markers[which(markers$gene %!in% ribo_genes),]
mt_genes <- grep(pattern = "^mt-", x = rownames(atlas_v1_NonMartinotti@assays$RNA@counts), value = TRUE)
markers<-markers[which(markers$gene %!in% mt_genes),]
markers<-markers[order(markers$avg_log2FC, decreasing = TRUE),]
markers<-markers[which(markers$p_val_adj < 0.01),]
markers$logFC<-markers$avg_log2FC
# select top 10 markers (ranked by logFC) for each cluster 
markers<-markers %>% group_by(cluster) %>% top_n(n = 10)
markers<-as.data.frame(markers)
gene_list<-unique(markers$gene)
# subset count matrix (RNA assay) on the marker gene list
NonMartinotti_sub<-CreateSeuratObject(atlas_v1_NonMartinotti@assays$RNA@counts[gene_list,])
NonMartinotti_sub<-AddMetaData(NonMartinotti_sub, atlas_v1_NonMartinotti@meta.data[colnames(NonMartinotti_sub),])
# Compute the accuracy matrix on the subsetted count matrix to assess the ability of marker genes to predict cluster identities.
accuracy_matrix_v1<-ComputeClassificationAccuracy(NonMartinotti_sub, "final_clusters_renamed")
saveRDS(accuracy_matrix_v1, "accuracy_matrix_NonMartinotti_atlas_v1.RDS")
# plot matrix
accuracy_matrix_v1<-accuracy_matrix_v1[colnames(accuracy_matrix_v1)[order(colnames(accuracy_matrix_v1))],colnames(accuracy_matrix_v1)[order(colnames(accuracy_matrix_v1))]]
pdf('RF_accuracy_matrix_NonMartinotti_atlas_v1.pdf')
gplots::heatmap.2(accuracy_matrix_v1, margins=c(7,7), key = FALSE, keysize=1, key.xlab="", key.title="Accuracy", trace = "none", 
                density.info = "none", col = cols, breaks = breaks, offsetRow=0.1, offsetCol=0.1, cexRow = 0.8, cexCol = 0.8,
                cellnote = accuracy_matrix_v1, notecex = 0.8, notecol = 'white', Colv = F, Rowv = F, dendrogram = "none")
dev.off()

# cluster matching - atlas v1 vs atlas v2
matches <- list("nMC1.0"=c("nMC1.1","nMC1.2","nMC1.3"), 
                "nMC2.0"=c("nMC2.1","nMC2.2","nMC2.3"), 
                "nMC3.0"=c("nMC3.0"), 
                "nMC4.0"=c("nMC4.0"), 
                "nMC5.1"=c("nMC5.0"), 
                "nMC5.2"=c("nMC5.0")
)
# compute delta on diagonal of accuracy matrices (atlas v2 - atlas v1)
delta_diag <- unlist(lapply(names(matches), function(cl){
  diag_2<-accuracy_matrix_v2[cl,cl]
  if(length(matches[[cl]]) > 1){
    #diag_1<-mean(diag(accuracy_matrix_v1[matches[[cl]],matches[[cl]]]))
    diag_1<-diag(accuracy_matrix_v1[matches[[cl]],matches[[cl]]])
  }else{
    diag_1<-accuracy_matrix_v1[matches[[cl]],matches[[cl]]]
  }
  delta_diag <- diag_2 - diag_1
  return(delta_diag)
}))
# compute delta on cross-talk of accuracy matrices (atlas v2 - atlas v1)
delta_crosstalk <- unlist(lapply(names(matches), function(cl){
  crosstalk_2<-sum(accuracy_matrix_v2[cl,!(colnames(accuracy_matrix_v2) %in% cl)])
  if(length(matches[[cl]]) > 1){
    crosstalk_1 <- unlist(lapply(matches[[cl]],function(x){
      crosstalk_1<-sum(accuracy_matrix_v1[x,!(colnames(accuracy_matrix_v1) %in% x)])
      return(crosstalk_1)
    }))
  }else{
    crosstalk_1<-sum(accuracy_matrix_v1[matches[[cl]],!(colnames(accuracy_matrix_v1) %in% matches[[cl]])])
  }
  delta_crosstalk <- crosstalk_2 - crosstalk_1
  return(delta_crosstalk)
}))
atlas_delta_summary_nMC <- data.frame(Family="Non-Martinotti",ClassIndex=c(0:(length(delta_crosstalk)-1)),diag_delta=delta_diag,ct_delta=delta_crosstalk)


######################
# atlas v2 vs atlas v1
######################
atlas_delta_summary <- rbind(atlas_delta_summary_LRP, atlas_delta_summary_MC, atlas_delta_summary_nMC)
# delta diagonal
means_diag <- df %>%
  group_by(Family) %>%
  summarise(mean_val = mean(diag_delta, na.rm = TRUE)) %>%
  ungroup()
pdf('Delta_diagonal_atlas_v1_vs_atlas_v2.pdf', width = 8, height = 6)
ggplot(df, aes(x = Family, y = diag_delta)) + ylim(-0.2,0.25) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
  stat_summary(fun = mean, geom = "bar", width = 0.6, fill = "#2C7BB6", alpha = 0.75) +
  geom_jitter(width = 0.12, size = 2.6, alpha = 0.9, color = "#313695") +
  geom_text(
    data = means_diag,
    aes(x = Family, y = mean_val, label = sprintf("%.3f", mean_val)),
    vjust = -0.6, fontface = "bold", color = "black"
  ) +
  labs(
    y = expression(Delta~"Diagonal Accuracy (v2 - v1)"),
    x = NULL,
    title = expression(Delta~"Diagonal Accuracy by Family")
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
    axis.text.x  = element_text(face = "bold", size = 13),
    panel.grid.major.x = element_blank()
  )
dev.off()
# delta cross-talk
means_ct <- df %>%
  group_by(Family) %>%
  summarise(mean_val = mean(ct_delta, na.rm = TRUE)) %>%
  mutate(vjust_lab = ifelse(mean_val >= 0, -0.6, 1.2)) %>%
  ungroup()

pdf('Delta_crosstalk_atlas_v1_vs_atlas_v2.pdf', width = 8, height = 6)
ggplot(df, aes(x = Family, y = ct_delta)) + ylim(-0.25,0.2) + 
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
  stat_summary(fun = mean, geom = "bar", width = 0.6, fill = "#D7191C", alpha = 0.75) +
  geom_jitter(width = 0.12, size = 2.6, alpha = 0.9, color = "#A50026") +
  geom_text(
    data = means_ct,
    aes(x = Family, y = mean_val, label = sprintf("%.3f", mean_val)),
    vjust = means_ct$vjust_lab,   # length matches rows in means_ct
    fontface = "bold", color = "black"
  ) +
  labs(
    y = expression(Delta~"Cross-talk (v2 - v1)"),
    x = NULL,
    title = expression(Delta~"Cross-talk by Family")
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
    axis.text.x  = element_text(face = "bold", size = 13),
    panel.grid.major.x = element_blank()
  )
dev.off()

###########################
# Spatial genes on atlas v2
###########################
# load genes from spatial data
gene_list<-readLines('gene_list_spatial.txt')
# subset count matrix (RNA assay) on the marker gene list
atlas_v2_sub<-CreateSeuratObject(atlas_v2@assays$RNA@counts[gene_list,])
atlas_v2_sub<-AddMetaData(atlas_v2_sub, atlas_v2@meta.data[colnames(atlas_v2_sub),])
# Compute the accuracy matrix on the subsetted count matrix to assess the ability of marker genes to predict cluster identities.
accuracy_matrix<-ComputeClassificationAccuracy(atlas_v2_sub, "final_clusters_renamed")
saveRDS(accuracy_matrix, "accuracy_matrix_NonMartinotti_atlas_v1.RDS")
# plot matrix
accuracy_matrix<-accuracy_matrix[colnames(accuracy_matrix)[order(colnames(accuracy_matrix))],colnames(accuracy_matrix)[order(colnames(accuracy_matrix))]]
pdf('RF_accuracy_matrix_atlas_v2_spatialGenes.pdf')
gplots::heatmap.2(accuracy_matrix, margins=c(7,7), key = FALSE, keysize=1, key.xlab="", key.title="Accuracy", trace = "none", 
                density.info = "none", col = cols, breaks = breaks, offsetRow=0.1, offsetCol=0.1, cexRow = 0.8, cexCol = 0.8,
                cellnote = accuracy_matrix, notecex = 0.8, notecol = 'white', Colv = F, Rowv = F, dendrogram = "none")
dev.off()
