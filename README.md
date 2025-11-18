# Dev-SST-atlas
This repository contains the code used to build the SST+ interneuron atlas and reproduce all associated results. 

## Atlas creation
The code is organized into three major modules:

  [**Module 1**](R_code/Module_1.R): includes the code for processing a single cell RNA seq dataset, performing quality filterings and isolating a cell population of interest (Sst+ interneurons in the provided code).

  [**Module 2**](Modules/Module_2): contains the code to assess whether datasets from different sources can be integrated using Canonical Correlation Analysis (CCA). One of the datasets must be set as the reference (the Atlas-v0 is used as the reference in the provided code).

  [**Module 3**](Modules/Module_3): includes the code to integrate all samples selected in Module 2 using CCA and to perform iterative clustering. A [Nextflow pipeline](Nextflow_pipeline/) is avialable for running the whole iterative clustering process.

Downstream analyses performed on Dev-SST-v2 atlas include [Waddington-OT analysis](R_code/WaddingtonOT_analysis.R), [Clustering validation](R_code/Cluster_validation.R), [Trajectory analysis](R_code/Trajectory_analysis.R), [Mapping to adult counterparts](R_code/MapMyCell.R) and [Saturation analysis](R_code/Saturation_analysis.R). 

## Spatial data analysis
The code to reproduce the cell annotation of spatial data is provided [here](R_code/Merfish_annotation.R).
