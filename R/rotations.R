#Rotations
run.umap.2 <- function(pca.obj, ndim=NULL, n.neighbors=15, n.components=2, min.dist=0.3, metric="cosine",seed=1234) {
  
  umap.config <- umap.defaults
  umap.config$n_neighbors = n.neighbors
  umap.config$min_dist = min.dist
  umap.config$metric = metric
  umap.config$n_components = n.components
  umap.config$random_state = seed
  umap.config$transform_state = seed
  
  if (is.null(ndim)) {
    ndim <- ncol(pca.obj$x)
  }
  
  ref.umap <- umap::umap(pca.obj$x[,1:ndim], config=umap.config)
  colnames(ref.umap$layout) <- c("UMAP_1","UMAP_2")
  return(ref.umap)
}

run.umap.uwot <- function(pca.obj, ndim=NULL, n.neighbors=15, n.components=2, min.dist=0.3, metric="cosine",seed=1234) {
  
  if (is.null(ndim)) {
    ndim <- ncol(pca.obj$x)
  }

  set.seed(seed)
  ref.umap <- uwot::umap(pca.obj$x[,1:ndim],
             metric=metric,
             min_dist=min.dist,
             n_neighbors = n.neighbors,
             ret_model=TRUE)
  
  colnames(ref.umap$embedding) <- c("UMAP_1","UMAP_2")
  
  return(ref.umap)
}

prcomp_seurat <- function(obj, assay=NULL, ndim=10, scale=TRUE) {
  
  if (is.null(assay)) {
    assay <- DefaultAssay(obj)
  }
  varfeat <- VariableFeatures(obj, assay=assay)
  mat <- GetAssayData(obj, assay=assay, slot="data")[varfeat,]
  refdata <- data.frame(t(as.matrix(mat)))
  
  refdata <- refdata[, sort(colnames(refdata))]
  ref.pca <- prcomp(refdata, rank. = ndim, scale. = scale, center = TRUE, retx=TRUE)
  
  #Save PCA rotation object
  obj@misc$pca_object <- ref.pca
  
  obj[["pca"]] <- CreateDimReducObject(embeddings=ref.pca$x, loadings=ref.pca$rotation, key = "PC_", assay = assay)
  return(obj)
}

apply.pca.obj.2 <- function(query, query.assay="RNA", pca.obj) {

  newdata <- data.frame(t(as.matrix(GetAssayData(query, assay=query.assay, slot="data"))))
  newdata <- newdata[ , order(names(newdata))]

  genes.use <-  sort(intersect(colnames(newdata), names(pca.obj$center)))

  newdata.var <- newdata[, genes.use]
  center.use <- pca.obj$center[genes.use]
  scale.use <- pca.obj$scale[genes.use]
  rotation.use <- pca.obj$rotation[genes.use,]

  npca <- scale(newdata.var, center.use, scale.use) %*% rotation.use

  return(npca)
}

apply.ica.obj <- function(query, query.assay="RNA", ica.obj) {

  newdata <- data.frame(t(as.matrix(GetAssayData(query, assay=query.assay, slot="data"))))
#  newdata <- data.frame(t(as.matrix(query@assays[[query.assay]]@data)))
  newdata <- newdata[ , order(names(newdata))]

  genes.use <-  sort(intersect(colnames(newdata), names(ica.obj$center)))

  newdata.var <- newdata[, genes.use]
  center.use <- ica.obj$center[genes.use]
  scale.use <- ica.obj$scale[genes.use]

  npca <- scale(newdata.var, center.use, scale.use) %*% ica.obj$K[genes.use,] %*% ica.obj$W
  colnames(npca) <- colnames(ica.obj$S)
  return(npca)
}

#dispatch to UMAP prediction method (complete of fast)
make.umap.predict <- function(ref.umap, fast.umap.predict=FALSE, ...) {
  
  if (fast.umap.predict) {
    nproj <- make.umap.predict.weighted.mean(ref.umap=ref.umap, ...)
  } else if (class(ref.umap) == "umap") {
    nproj <- make.umap.predict.2(ref.umap=ref.umap,
                                 method="umap", ...)
  } else if (!is.null(ref.umap$embedding)) {
    nproj <- make.umap.predict.2(ref.umap=ref.umap,
                                 method="uwot", ...)
  } else {
    warning("No UMAP-predict model available. Using fast.umap.predict approximation.")
    nproj <- make.umap.predict.weighted.mean(ref.umap=ref.umap, ...)
  }
  return(nproj)
}

#UMAP predict usign the umap package
make.umap.predict.2 <- function(ref.umap,
                                query,
                                query.assay="RNA",
                                pca.obj,
                                pca.query.emb=NULL,
                                method="uwot") {

  #if PCA query cell embeddings have been pre-calculated, read them from variable
  if (is.null(pca.query.emb)) {
    pca.query.emb <- apply.pca.obj.2(query=query, query.assay=query.assay, pca.obj=pca.obj)
  }  

  pca.dim <- dim(ref.umap$data)[2]
  
  if (method == "umap") {
    nproj.umap <- umap:::predict.umap(ref.umap, pca.query.emb[,1:pca.dim])
  } else if (method == "uwot") {
    nproj.umap <- uwot::umap_transform(pca.query.emb[,1:pca.dim], model = ref.umap)
  } else {
    stop("Unsupported UMAP method.")
  }
  return(nproj.umap)
}

#Fast projection mode: assign UMAP coordinates based on nearest neighbors in PCA space
make.umap.predict.weighted.mean <- function(ref.umap, query,
                                            query.assay="RNA",
                                            pca.obj, pca.query.emb=NULL,
                                            k=8) {
  
  if (is.null(pca.query.emb)) {
    pca.query.emb <- apply.pca.obj.2(query=query, query.assay=query.assay, pca.obj=pca.obj)
  }
  
  ref.space <- ref.umap$data
  pca.dim <- ncol(ref.space)
  query.space <- pca.query.emb[,1:pca.dim]
  
  nn.ranked <- Seurat:::NNHelper(data=ref.space, query=query.space, k = k, method = "rann")
  
  cellnames <- rownames(query.space)
  nproj.umap <- matrix(data = NA, nrow = length(cellnames), ncol = 2,
                       dimnames = list(cellnames, c("UMAP_1","UMAP_2")))
  
  for (cell in 1:length(cellnames)) {
    row <- exp(-nn.ranked@nn.dist[cell,])  #calculate exp(-dist) as weights for nearest neighbors
    weights = row/sum(row)
    nproj.umap[cell,] = weights %*% ref.umap$layout[nn.ranked@nn.idx[cell,],]  #assign UMAP coordinates of (weighted) neighbors
  }
  return(nproj.umap)

}

run.ica <- function(object, assay="integrated", ndim=50) {
  
  set.seed(1234)
  varfeat <- VariableFeatures(object, assay=assay)
  
  x <- scale(Matrix::t(GetAssayData(object, assay=assay, slot="data")[varfeat,]))
  set.seed(1234)
  ref.ica <- fastICA(x, n.comp=ndim, row.norm=T, maxit=1000, verbose=FALSE, tol=1e-13, method="R")
  
  ids <- paste0("ICA_", seq_len(ncol(ref.ica$K)))
  
  rownames(ref.ica$X) <- colnames(object)
  colnames(ref.ica$X) <- varfeat
  rownames(ref.ica$K) <- varfeat
  colnames(ref.ica$K) <- ids
  rownames(ref.ica$A) <- ids
  colnames(ref.ica$A) <- colnames(ref.ica$X)
  rownames(ref.ica$S) <- colnames(object)
  colnames(ref.ica$S) <- ids
  
  ref.ica$center <- attr(x,"scaled:center")
  ref.ica$scale <- attr(x,"scaled:scale")
  
  object[["ica"]] <- CreateDimReducObject(embeddings=ref.ica$S, loadings=t(ref.ica$A), key = "ICA_", assay = assay)
  object@misc$ica <- ref.ica
  return(object)
}

