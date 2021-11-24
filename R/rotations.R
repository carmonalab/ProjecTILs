#Rotations
apply.pca.obj.2 <- function(query, query.assay="RNA", pca.obj) {

  newdata <- data.frame(t(as.matrix(query@assays[[query.assay]]@data)))
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

  newdata <- data.frame(t(as.matrix(query@assays[[query.assay]]@data)))
  newdata <- newdata[ , order(names(newdata))]

  genes.use <-  sort(intersect(colnames(newdata), names(ica.obj$center)))

  newdata.var <- newdata[, genes.use]
  center.use <- ica.obj$center[genes.use]
  scale.use <- ica.obj$scale[genes.use]
  #rotation.use <- ica.obj$rotation[genes.use,]

  #ref.ica$X %*% ref.ica$K %*% ref.ica$W
  #npca <- scale(newdata.var, center.use, scale.use) %*% rotation.use
  npca <- scale(newdata.var, center.use, scale.use) %*% ica.obj$K[genes.use,] %*% ica.obj$W
  colnames(npca) <- colnames(ica.obj$S)
  return(npca)
}

#dispatch to UMAP prediction method (complete of fast)
make.umap.predict <- function(ref.umap, fast.mode=FALSE, ...) {
  
  if (!class(ref.umap) == "umap") {  #enforce fast.mode
     fast.mode <- TRUE
  }
  
  if (fast.mode) {
    nproj <- make.umap.predict.weighted.mean(ref.umap=ref.umap, ...)
  } else {
    nproj <- make.umap.predict.2(ref.umap=ref.umap, ...)
  }
  return(nproj)
}

#UMAP predict usign the umap package
make.umap.predict.2 <- function(ref.umap, query, query.assay="RNA", pca.obj, pca.query.emb=NULL) {

  if (is.null(pca.query.emb)) {
    pca.query.emb <- apply.pca.obj.2(query=query, query.assay=query.assay, pca.obj=pca.obj)
  }  #if PCA query cell embeddings have been pre-calculated, read them from variable

  pca.dim <- dim(ref.umap$data)[2]
  nproj.umap <- umap:::predict.umap(ref.umap, pca.query.emb[,1:pca.dim])

  return(nproj.umap)
}

#Fast projection mode: assign UMAP coordinates based on nearest neighbors in PCA space
make.umap.predict.weighted.mean <- function(ref.umap, query, query.assay="RNA", pca.obj, pca.query.emb=NULL, k=8) {
  
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
  require(fastICA)
  set.seed(1234)
  x <- scale(Matrix::t(object@assays[[assay]][object@assays[[assay]]@var.features,]))
  set.seed(1234)
  ref.ica <- fastICA(x, n.comp=ndim, row.norm=T, maxit=1000, verbose=F, tol=1e-13, method="R")
  
  colnames(ref.ica$X) <- ref@assays$integrated@var.features
  rownames(ref.ica$X) <- colnames(ref)
  rownames(ref.ica$K) <- ref@assays$integrated@var.features
  colnames(ref.ica$A) <- colnames(ref.ica$X)
  colnames(ref.ica$S) <- paste0("ICA_", seq_len(ncol(ref.ica$S)))
  
  ref.ica$center <- attr(x,"scaled:center")
  ref.ica$scale <- attr(x,"scaled:scale")
  
  object[["ica"]] <- CreateDimReducObject(embeddings=ref.ica$S, loadings=t(ref.ica$A), key = "ICA_", assay = assay)
  object@misc$ica <- ref.ica
  return(object)
}

