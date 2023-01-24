# ProjecTILs - Functions

* `load.reference.map`   Load or download the reference map for dataset projection. By the default it downloads the reference atlas of tumour-infiltrating lymphocytes (TILs).

* `read.sc.query`   Load a query expression matrix to be projected onto the reference atlas. Several formats (10x, hdf5, raw and log counts, etc.) are supported - see type parameter for details 

* `Run.ProjecTILs`	A wrapper for `make.projection` and `cellstate.predict` (see below). It returns the query embedded in the PCA and UMAP space of the reference, and predict cell state labels for the query cells.

* `ProjecTILs.classifier`	Similarly to `Run.ProjecTILs`, it perform reference projection and cell state prediction. However, the query embeddings are left intact and only cell state labels are returned. Useful if you wish to use ProjecTILs as a classifier and annotate your data in their native space.

* `make.projection`   Project a single-cell RNA-seq dataset onto a reference map of cellular states.

* `celltype.heatmap`   Generate a averaged gene expression heatmap from a Seurat object

* `plot.projection`   Plots the UMAP representation of the reference map, together with the projected coordinates of a query dataset.

* `cellstate.predict`   Use a nearest-neighbor algorithm to predict a feature (e.g. the cell state) of the query cells.

* `plot.statepred.composition`   Makes a barplot of the frequency of cell states in a query object.

* `plot.states.radar`   Makes a radar plot of the expression level of a specified set of genes.

* `find.discriminant.dimensions`   Searches PCA or ICA dimensions where the query set deviates the most from a control set or from the reference map.

* `plot.discriminant.3d`   Add an extra dimension to the reference map  to explore additional axes of variability in a query dataset compared to the reference map.

* `find.discriminant.genes` Performs differential expression analysis between a projected query and a control (either the reference map or a control sample), for
a selected reference subtype. Useful to detect whether specific cell states over/under-express genes between conditions or with respect to the reference.

* `make.reference`	Converts a Seurat object into a custom reference map for ProjecTILs.

* `recalculate.embeddings` After projection of query data into a reference, you may want to recalculate the low-dimensional embeddings accounting for the new data. The resulting object can be used as a new reference. 

* `merge.Seurat.embeddings` Given two Seurat objects, merge counts and data as well as dim reductions (PCA, UMAP, ICA, etc.)

* `compute_silhouette` Given a projected object and its reference, calculate silhouette coefficient for query cells with respect to reference cells with the same cell labels.

Find more information, syntax and examples using the R help function e.g. `?Run.ProjecTILs`

