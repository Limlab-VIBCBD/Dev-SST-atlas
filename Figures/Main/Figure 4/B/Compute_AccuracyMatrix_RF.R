
#function to apply RF validation and score the class membership. 
ComputeClassificationAccuracy <- function(seuratObject,
                            clustering,
                            minClusterRep=5,
                            ntree=100,
                            nGene=10,
                            nIteration=10,
                            nCore=10,
                            ...){
  
  cells<-rownames(seuratObject[[clustering]])
  seuratObject <- subset(seuratObject, cells = cells[which(!is.na(seuratObject[[clustering]]))] )
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


