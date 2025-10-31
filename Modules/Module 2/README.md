# Integrability tests
Only datasets containing at least 1'000 Sst⁺ cells were considered for integration into the Dev-SST atlas.
The integrability of each dataset with the reference atlas (Atlas-v0) was evaluated using three criteria:

+ Anchor Points — Number of Seurat integration anchors identified between the query dataset and Atlas-v0 (reference). To compute the anchors, 1'000 cells were randomly sampled from the query dataset multiple times to ensure results were independent of dataset size. The mean number of anchors across all subsamplings was calculated, and datasets with values above the predefined acceptance threshold were considered to have passed this criterion.

+ Neighborhood Composition — After integrating the query dataset with Atlas-v0 using Canonical Correlation Analysis (CCA), a k-nearest neighbor (kNN) graph was computed in the integrated PCA space using Euclidean distance, for increasing values of k. For each cell in the query dataset, we evaluated the fraction of Atlas-v0 cells among its k-nearest neighbors. Datasets in which at least 50% of the cells had at least half of the global fraction of reference atlas cells among their neighbors (k = 60) were considered to have passed this criterion.

+ kBET Rejection Rate — The kBET package was used to compute the rejection rate, assessing how well the query dataset and Atlas-v0 were mixed in the integrated PCA space. A rejection rate below 50% was set as the acceptance threshold.

Only datasets meeting all three criteria (anchor points, neighborhood composition, and kBET rejection rate) were included in the construction of the Dev-SST-v2 atlas.

You can find the code to reproduce the integrability tests [here](). \
The code used to compute the anchor point threshold is available [here]().