# Dev-SST-atlas
This repository contains the code used to build the SST+ interneuron atlas and to visualize all associated results. 

## Atlas creation
The code is organized into three major modules:

  [**Module 1**](R_code/Module_1.R): includes the code for processing a single cell RNA seq dataset, performing quality filterings and isolating a cell population of interest (Sst+ interneurons in the reported code).

  [**Module 2**](Modules/Module%202): contains the code to assess whether datasets from different sources can be integrated using Canonical Correlation Analysis (CCA). One of the datasets must be set as the reference (the Atlas-v0 is used as the reference in the reported code).

  [**Module 3**](Modules/Module%203): includes the code to integrate all samples selected in Module 2 using CCA and to perform iterative clustering. A [Nextflow pipeline](Nextflow_pipeline/) is avialable for running the whole iterative clustering process.

Downstream analyses performed on Dev-SST-v2 atlas include [Waddington-OT analysis](), [Clustering validation](), [Trajectory analysis](), [Mapping to adult counterparts]() and [Saturation analysis]()

## Spatial data analysis

## Paper figures (Main)
### Figure 1
+ Figure 1B
+ Figure 1C & Figure 1D
+ Figure 1E
### Figure 2
### Figure 3
### Figure 4
### Figure 5
### Figure 6
### Figure 7


