########################
# Load required packages
########################
library(stringr)
library(Seurat)
library(SeuratWrappers)
library(SeuratDisk)
library(dittoSeq)
library(gtools)
library(dplyr)
library(R.utils)


####################
# MapMyCell to adult
####################
# Prepare files for MapMyCell
counts<-atlas_v2@assays$RNA$counts
obj<-CreateSeuratObject(counts)
SaveH5Seurat(obj, filename = "atlas_v2_raw_counts.h5Seurat", overwrite = TRUE)
Convert("atlas_v2_raw_counts.h5Seurat", dest = "h5ad")

# Assign cell types using MapMyCells
# These next steps are performed OUTSIDE of R in the MapMyCells web application.
# The steps to MapMyCells are as follows:
#1.  Go to (<https://knowledge.brain-map.org/mapmycells/process/>).\
#2.  (Optional) Log in to MapMyCells.\
#3.  Upload 'atlas_v2_raw_counts_MapMyCell.h5ad' to the site via the file system or drag and drop (Step 1).\
#4.  Choose "10x Whole Mouse Brain (CCN20230722)" as the "Reference Taxonomy" (Step 2).\
#5.  Choose the desired "Mapping Algorithm" (in this case "Hierarchical Mapping").\
#6.  Click "Start" and wait \~5 minutes. (Optional) You may have a panel on the left that says "Map Results" where you can also wait for your run to finish.\
#7.  When the mapping is complete, you will have an option to download the "tar" file with mapping results. If your browser is preventing popups, search for small folder icon to the right URL address bar to enable downloads.\
#8.  Unzip this file, which will contain three files: "validation_log.txt", "[NUMBER].json", and "[NUMBER].csv". [NUMBER].csv contains the mapping results, which you need. The validation log will give you information about the the run itself and the json file will give you extra information about the mapping (you can ignore both of these files if the run completes successfully).\
#9.  Copy [NUMBER].csv to your current working directory and rename it "mapmycell.csv").\
#10. You can now go back to R and continue the script below.

mapmycell<-read.csv('mapmycell.csv',comment.char="#")
mapmycell$cell_id <- gsub("\\.", "-", mapmycell$cell_id)
rownames(mapmycell) <- mapmycell$cell_id

### Add mapMyCell annotation to the atlas and filter high confidence cells
atlas_v2$class_name <- mapmycell[colnames(atlas_v2),"class_name"]
atlas_v2$subclass_name <- mapmycell[colnames(atlas_v2),"subclass_name"]
atlas_v2$subclass_bootstrapping_probability <- mapmycell[colnames(atlas_v2),"subclass_bootstrapping_probability"]
atlas_v2$supertype_name <- mapmycell[colnames(atlas_v2),"supertype_name"]
atlas_v2$supertype_bootstrapping_probability <- mapmycell[colnames(atlas_v2),"supertype_bootstrapping_probability"]
atlas_v2$cluster_name <- mapmycell[colnames(atlas_v2),"cluster_name"]

atlas_sub<-subset(atlas_v2, subset= class_name=="07 CTX-MGE GABA")
atlas_sub<-subset(atlas_sub, subset= subclass_name=="053 Sst Gaba")
atlas_sub<-subset(atlas_sub, subset=supertype_bootstrapping_probability > 0.9)

### one to one match clusters::supertypes
obj_sub<-subset(atlas_sub, subset=major_label_transferAnchors_BaseAtlas == "Martinotti")
#obj_sub<-subset(obj_sub, subset=major_label_transferAnchors_BaseAtlas == "Non-Martinotti") ## same analysis for both datasets
clusters_list <- names(table(obj_sub@meta.data$final_clusters_renamed))
matching_supertype <- lapply(clusters_list, function(cl){
  frac <- table(obj_sub@meta.data[which(obj_sub@meta.data$final_clusters_renamed %in% cl),"supertype_name"])/length(obj_sub@meta.data[which(obj_sub@meta.data$final_clusters_renamed %in% cl),"supertype_name"])
  matching_supertype <- names(frac)[which(frac > 0.65)]
  return(matching_supertype)
})
names(matching_supertype) <- clusters_list
matching_supertype <- matching_supertype[lapply(matching_supertype,length)>0]

supertypes_list <- names(table(obj_sub@meta.data$supertype_name))
matching_cluster <- lapply(supertypes_list, function(cl){
  frac <- table(obj_sub@meta.data[which(obj_sub@meta.data$supertype_name %in% cl),"final_clusters_renamed"])/length(obj_sub@meta.data[which(obj_sub@meta.data$supertype_name %in% cl),"final_clusters_renamed"])
  frac<-frac[which(frac > 0.3)]
  if(length(frac) > 0){
    matching_cluster <- names(frac)[which(frac == max(frac))]
  }
  return(matching_cluster)
})
names(matching_cluster) <- supertypes_list
matching_cluster <- matching_cluster[lapply(matching_cluster,length)>0]

matching_supertype<-unlist(matching_supertype)
matching_supertype <- matching_supertype[matching_supertype %in% names(unlist(matching_cluster))]
matching_supertype <- matching_supertype[names(matching_supertype) %in% unlist(matching_cluster)]


###################################
# Custom MapMyCell - celltypemapper
###################################
## Custom MapMyCell on a dowloaded dataset
## Download Yao 2023 dataset (CNN20230722) from these links: https://allen-brain-cell-atlas.s3.us-west-2.amazonaws.com/index.html#expression_matrices/WMB-10Xv2/20230630/, https://allen-brain-cell-atlas.s3.us-west-2.amazonaws.com/index.html#expression_matrices/WMB-10Xv3/20230630/ and for Multi https://allen-brain-cell-atlas.s3.us-west-2.amazonaws.com/index.html#expression_matrices/WMB-10XMulti/20230830/
# cl.df_CCN202307220.xlsx - Download annotation table here: https://portal.brain-map.org/explore/cell-type-references-and-algorithms 
# cell_metadata.csv - Download metadata here https://allen-brain-cell-atlas.s3.us-west-2.amazonaws.com/index.html#metadata/WMB-10X/20231215/

# Create seurat object and h5ad file for a subset of cells (if needed)
clusters_sub<-allen_csv[which(allen_csv$subclass_id %in% c(53,56)),'cl'] # in this example create object with cells from subclasses 053 and 056
metadata_sub<-subset(metadata, subset=cluster_alias %in% clusters_sub)
rownames(metadata_sub) <-metadata_sub$cell_label
# read raw counts
file_list<-names(table(metadata_sub$feature_matrix_label))
seurat_objects<-lapply(file_list,function(x){
  cells_sub <- metadata_sub[which(metadata_sub$feature_matrix_label == x),'cell_label']
  file<-paste0(x,"-raw.h5ad")
  tmp<-read_h5ad(file)
  counts<-tmp$X
  counts<-t(counts[cells_sub,])
  genes   <- h5read(file,"/var/gene_symbol")
  rownames(counts) <-genes
  obj <- CreateSeuratObject(counts=counts, min.cells = 0, min.genes = 0)
  obj <- AddMetaData(obj,metadata_sub[which(metadata_sub$feature_matrix_label == x),])
  rm(tmp)
  rm(counts)
  gc()
  return(obj)
})
names(seurat_objects) <- file_list
# Merge seurat objects
merged_obj <-Reduce(
  f = function(x, y) {merge(x, y, merge.data = FALSE)},
  x = seurat_objects # list of Seurat objects
)
# Add metadata info to the merged object
classes <- allen_csv[which(allen_csv$subclass_id %in% c(53,56)),'class_id_label']
names(classes) <- as.character(allen_csv[which(allen_csv$subclass_id %in% c(53,56)),'cl'])
merged_obj$class_name <- as.character(merged_obj$cluster_alias)
tmp<-classes[merged_obj$class_name]
names(tmp) <- colnames(merged_obj)
merged_obj$class_name <- tmp

subclasses <- allen_csv[which(allen_csv$subclass_id %in% c(53,56)),'subclass_id_label']
names(subclasses) <- as.character(allen_csv[which(allen_csv$subclass_id %in% c(53,56)),'cl'])
merged_obj$subclass_name <- as.character(merged_obj$cluster_alias)
tmp<-subclasses[merged_obj$subclass_name]
names(tmp) <- colnames(merged_obj)
merged_obj$subclass_name <- tmp

supertypes <- allen_csv[which(allen_csv$subclass_id %in% c(53,56)),'supertype_id_label']
names(supertypes) <- as.character(allen_csv[which(allen_csv$subclass_id %in% c(53,56)),'cl'])
merged_obj$supertype_name <- as.character(merged_obj$cluster_alias)
tmp<-supertypes[merged_obj$supertype_name]
names(tmp) <- colnames(merged_obj)
merged_obj$supertype_name <- tmp

clusters <- allen_csv[which(allen_csv$subclass_id %in% c(53,56)),'cluster_id_label']
names(clusters) <- as.character(allen_csv[which(allen_csv$subclass_id %in% c(53,56)),'cl'])
merged_obj$cluster_name <- as.character(merged_obj$cluster_alias)
tmp<-clusters[merged_obj$cluster_name]
names(tmp) <- colnames(merged_obj)
merged_obj$cluster_name <- tmp

# Add umap to seurat object
UMAP_coordinates<-merged_obj@meta.data[,c("x","y")]
colnames(UMAP_coordinates) <- c("UMAP_1","UMAP_2")
merged_obj[['umap']] <- CreateDimReducObject(embeddings = as.matrix(UMAP_coordinates), key = "UMAP_", global = T, assay = "RNA")

saveRDS(merged_obj,"Yao_2023_053_056_subclasses.RDS")
SaveH5Seurat(merged_obj, filename = "Yao_2023_053_056_subclasses.h5Seurat", overwrite = TRUE)
Convert("Yao_2023_053_056_subclasses.h5Seurat", dest = "h5ad") # create h5ad file that could be used as input by precompute_stats_scrattch function

# celltypemapper install is needed (https://github.com/AllenInstitute/cell_type_mapper)
# run celltypemapper from command line
# Create reference 
#python -m cell_type_mapper.cli.precompute_stats_scrattch --h5ad_path Yao_2023.h5ad --hierarchy '["class_name","subclass_name","supertype_name","cluster_name"]' --output_path Reference_markers_sub/precomputed_stats_Yao.h5
#python -m cell_type_mapper.cli.reference_markers --precomputed_path_list '["Reference_markers_sub/precomputed_stats_Yao.h5"]'  --output_dir Reference_markers_sub/ --tmp tmp_folder/ 
#python -m cell_type_mapper.cli.query_markers --reference_marker_path_list '["Reference_markers_sub/reference_markers.h5"]' --output_path Reference_markers_sub/reference_markers.json

# Map atlas to reference 
#python -m cell_type_mapper.cli.from_specified_markers --query_path atlas_v2_raw_counts.h5ad --extended_result_path MapMyCell/Reference/Reference_markers_sub/results.json --precomputed_stats.path MapMyCell/Reference/Reference_markers_sub/precomputed_stats_Yao.h5 --query_markers.serialized_lookup MapMyCell/Reference/Reference_markers_sub/reference_markers.json --type_assignment.normalization raw --csv_result_path mapmycell.csv


#########################
# celltypemapper for LRPs
#########################
# Create reference LRP atlas
LRP_subset <- subset(atlas_v2,subset=major_label_transferAnchors_BaseAtlas == "LRP")
LRP_subset <- subset(LRP_subset,subset=final_clusters_renamed %in% c("LRP2.1","LRP1.0"))
counts<-LRP_subset@assays$RNA$counts
LRP_subset$clusters<-LRP_subset$final_clusters_renamed
metadata<-as.data.frame(LRP_subset@meta.data[,c("clusters")])
colnames(metadata) <- 'clusters'
rownames(metadata) <- colnames(LRP_subset)
obj<-CreateSeuratObject(counts)
obj<-AddMetaData(obj, metadata)
SaveH5Seurat(obj, filename = "LRP_sub_MapMyCell.h5Seurat", overwrite = TRUE)
Convert("LRP_sub_MapMyCell.h5Seurat", dest = "h5ad")
# celltypemapper - command line: create reference 
#python -m cell_type_mapper.cli.precompute_stats_scrattch --h5ad_path LRP_sub_MapMyCell.h5ad --hierarchy '["clusters"]' --output_path Reference_markers_sub/precomputed_stats.h5
#python -m cell_type_mapper.cli.reference_markers --precomputed_path_list '["Reference_markers_sub/precomputed_stats.h5"]'  --output_dir Reference_markers_sub/ --tmp tmp_folder/ 
#python -m cell_type_mapper.cli.query_markers --reference_marker_path_list '["Reference_markers_sub/reference_markers.h5"]' --output_path Reference_markers_sub/reference_markers.json

# MAP TO LRP REFERENCE 
# Example on SST+ NPY+ human dataset
so <- readRDS('SeuratObj_NPY.SST.RDS') # downloaded from https://dev-ctx-meta-atlas.cells.ucsc.edu and further filtered
## for mapping to mouse data 
convert <- read.csv("mouse_human_marmoset_macaque_orthologs_20231113.csv") # can be downloaded here:https://github.com/AllenInstitute/GeneOrthology/blob/main/csv/mouse_human_marmoset_macaque_orthologs_20231113.csv
# Do the conversion by gene symbol 
convert_by_symbol <- convert[!(is.na(convert$human_Symbol)|is.na(convert$mouse_Symbol)),c("human_Symbol","mouse_Symbol")] # Remove NAs from conversion table
convert_by_symbol <- convert_by_symbol[is.element(convert_by_symbol$human_Symbol,rownames(so@assays$RNA@counts)),] # Remove genes not in data matrix
dataOut <- so@assays$RNA@counts[match(convert_by_symbol$human_Symbol,rownames(so@assays$RNA@counts)),] # Subset data to include only genes with mouse othologs
rownames(dataOut) <- convert_by_symbol$mouse_Symbol # Convert gene names to mouse
obj<-CreateSeuratObject(dataOut)
SaveH5Seurat(obj, filename = "NPY_SST_raw_counts_mouseOrt.h5Seurat", overwrite = TRUE)
Convert("NPY_SST_raw_counts_mouseOrt.h5Seurat", dest = "h5ad")
# celltypemapper - command line: map to reference 
#python -m cell_type_mapper.cli.from_specified_markers --query_path NPY_SST_raw_counts_mouseOrt.h5ad --extended_result_path results.json --precomputed_stats.path precomputed_stats.h5 --query_markers.serialized_lookup reference_markers.json --type_assignment.normalization log2CPM --csv_result_path results.csv
