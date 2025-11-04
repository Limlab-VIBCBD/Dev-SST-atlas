# Sample integration and clustering
Module 3 is subdivided into two major parts:

### QC and filtering, integration and annotation to major cell-types 
All samples selected in module 2 were integrated with CCA. After a first step of integration and clustering, computed with standard Seurat analysis, quality checks were performed and clusters eriched in stressed/low-quality cells were discarded. A new integration was performed on retained cells and UMAP was computed on integrated PCs. \
Cells were annotated using Label Trasfer on Atlas-v0 to three major classes: Long-range projecting (LRP) neurons, Martinotti interneurons and Non-Martinotti interneurons. \
You can find the code to reproduce this part [here]()

### Iterative clustering
Iterative clustering analysis was performed separately for each major cell type (LRP, Martinotti, and Non-Martinotti).\
To make this crucial step in the creation of the Dev-SST Atlas more reproducible and accessible to other users and applicable to different datasets, we developed a [**Nextflow pipeline**](). We recommend completing all preceding steps before running the iterative clustering to ensure the quality of the selected datasets and cells. However, the iterative clustering pipeline can also be executed independently, without running all prior steps.

