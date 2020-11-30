# Projecting single-cell transcriptomics data onto a reference T cell atlas to interpret immune responses

<span title="Shuttling T cells into a reference transcriptomic space. A high-dimensional odyssey to interpret  immune responses" >
<p align="center">
  <img height="100" src="docs/projectils_logo_W_square.png">
</p>
</span>

`ProjecTILs` is a computational method to project scRNA-seq data into a reference map of T cells, enabling their direct comparison in a stable, annotated system of coordinates. Because new cells are embedded in the same space of the reference, ProjecTILs enables the classification of new cells into annotated, discrete states, but also over a continuous space of intermediate states. 
By comparing multiple samples over the same reference map, and across alternative embeddings, our method allows exploring the effect of cellular perturbations (e.g. as the result of therapy or genetic engineering) in terms of transcriptional states and altered genetic programs.

We have constructed two cross-study murine T cell reference atlases for ProjecTILs: the first describing tumor-infiltrating T lymphocytes (TILs), the second characterizing virus-specific CD8 T cells in acute and chronic infection. 

Find the installation instructions for the package below, and a vignette detailing its functions at [Tutorial (html)](https://carmonalab.github.io/ProjecTILs/tutorial.html) and [Tutorial (repository)](https://gitlab.unil.ch/carmona/ProjecTILs.demo)

For real-life applications, check out our list of [ProjecTILs Case Studies](https://carmonalab.github.io/ProjecTILs_CaseStudies/)

If you prefer to avoid installing R packages, you can run `ProjecTILs` in Docker. A ready-to-use Docker image with usage instructions is available on [DockerHub](https://hub.docker.com/repository/docker/mandrea1/projectils_demo)

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

if (!requireNamespace("TILPRED", quietly = TRUE)) {
  install.packages(c("doParallel","doRNG"))
  BiocManager::install(c("AUCell","SingleCellExperiment"))
  remotes::install_github("carmonalab/TILPRED")
}

remotes::install_github("carmonalab/ProjecTILs")


```

### Test the package

Load sample data and test your installation:
```
library(ProjecTILs)
data(query_example_seurat)
make.projection(query_example_seurat, skip.normalize=T)
```

On your first run, the `make.projection` call will download the reference TIL atlas used as default map. This may take some time depending on your connection, but it is only necessary on the first run.


### Data projection TUTORIAL

Find a step-by-step tutorial for `ProjecTILs` at: [ProjecTILs tutorial](https://carmonalab.github.io/ProjecTILs/tutorial.html)

To run the code of the tutorial on your machine, download the demo repository: [ProjecTILs tutorial repo](https://gitlab.unil.ch/carmona/ProjecTILs.demo) or obtain a [Docker image](https://hub.docker.com/repository/docker/mandrea1/projectils_demo) with all dependencies pre-installed.

For real-life applications, check out our list of [ProjecTILs Case Studies](https://carmonalab.github.io/ProjecTILs_CaseStudies/)

### Documentation

See a description of the functions implemented in ProjecTILs at: [ProjecTILs functions](docs/functions.md)

### Reference Atlases

Pre-computed reference atlases are available at:

* TIL atlas: [https://doi.org/10.6084/m9.figshare.12478571](https://doi.org/10.6084/m9.figshare.12478571) and interactive iSEE web app [http://TILatlas.unil.ch](http://TILatlas.unil.ch)

* viral infection CD8 T cell atlas: [https://doi.org/10.6084/m9.figshare.12489518](https://doi.org/10.6084/m9.figshare.12489518) and web app 
http://virustcellatlas.unil.ch/

If you wish to use your own **custom reference atlas**, follow this vignette to prepare it in a format that can be understood by ProjecTILs: [Building a custom reference atlas for ProjecTILs](https://carmonalab.github.io/ProjecTILs/build_ref_atlas.html)

### Troubleshooting 

* If a warning message prevented *remotes* from installing the package, try:
```Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS="true")```

* For analyzing datasets composed of multiple batches (e.g. different subjects, technologies), we recommend projecting each batch separately, by providing ProjecTILs a list of Seurat objects as input, e.g.:
```
data.seurat.list <- SplitObject(data.seurat, split.by = "batch")
query.projected.list <- make.projection(data.seurat.list)
```



### Citation

Projecting single-cell transcriptomics data onto a reference T cell atlas to interpret immune responses. Massimo Andreatta, Jesus Corria Osorio, Soren Muller, Rafael Cubas, George Coukos,  Santiago J Carmona. 2020 https://doi.org/10.1101/2020.06.23.166546

<p align="center">
  <img height="80" src="docs/projectils_logo_W_square.png">
</p>
