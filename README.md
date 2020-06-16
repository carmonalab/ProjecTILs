# ProjecTILs - Projecting scRNA-seq data onto a reference map of T cell transcriptomic states

`ProjecTILs` is a computational method to project new data sets into a reference map of T cells, enabling their direct comparison in a stable, annotated system of coordinates. Because new cells are embedded in the same space of the reference, ProjecTILs enables the classification of new cells into annotated, discrete states, but also over a continuous space of intermediate states. By comparing multiple samples over the same map, and across alternative embeddings, the method allows exploring the effect of genetic perturbations (e.g. as the result of therapy) and identifying genetic programs significantly altered in the query compared to a control set or to the reference map.

Find the installation instructions for the package below, and a vignette detailing its functions at [Tutorial (html)](https://carmonalab.github.io/ProjecTILs/) and [Tutorial (repository)](https://gitlab.unil.ch/carmona/ProjecTILs.demo)

### Package Installation

To install `ProjecTILs` directly from its Git repository, run the following code from within R or RStudio:
```
if (!requireNamespace("remotes")) install.packages("remotes")
library(remotes)

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

if (!requireNamespace("Seurat", quietly = TRUE)) {
   BiocManager::install('multtest')
   install.packages("Seurat")
}

remotes::install_github("carmonalab/TILPRED")
remotes::install_github("carmonalab/ProjecTILs")
```

### Test the package

Load sample data and test your installation:
```
library(projecTILs)
data(query_example_seurat)
make.projection(query_example_seurat, skip.normalize=T)
```

On your first run, the `make.projection` call will download the reference TIL atlas used as default map. This may take some time depending on your connection, but it is only necessary on the first run.


### Data projection TUTORIAL

Find a step-by-step tutorial for `ProjecTILs` at: [ProjecTILs tutorial](https://carmonalab.github.io/ProjecTILs/)

To run the code of the tutorial on your machine, download the demo repository: [ProjecTILs tutorial repo](https://gitlab.unil.ch/carmona/ProjecTILs.demo) or obtain a Docker image (TODO) with all dependencies pre-installed.

### Documentation

See a description of the functions implemented in ProjecTILs at: [ProjecTILs functions](docs/functions.md)

### Atlases

Pre-computed atlases are available at:

* TIL atlas: [https://doi.org/10.6084/m9.figshare.12478571](https://doi.org/10.6084/m9.figshare.12478571)

* LCMV atlas: [https://doi.org/10.6084/m9.figshare.12489518](https://doi.org/10.6084/m9.figshare.12489518)
