

filterCells <- function(query.object, filter=TRUE){
  sce <- as.SingleCellExperiment(query.object)
  sce.pred <- predictTilState(sce)
  query.object <- AddMetaData(query.object, metadata=sce.pred$predictedState, col.name = "TILPRED")
  query.object <- AddMetaData(query.object, metadata=sce.pred$cyclingScore, col.name = "cycling.score")

  if (filter) {
     cellRemove <- colnames(query.object)[!query.object$TILPRED %in% c("Non-Tcell","unknown")]
     print(paste(ncol(query.object)-length(cellRemove), "out of", ncol(query.object),
                 "(",round((ncol(query.object)-length(cellRemove))/ncol(query.object)*100),"% )",
                 "non-pure T cells removed.  Use filter.cells=FALSE to avoid pre-filtering (NOT RECOMMENDED)"))
     query.object <- subset(query.object, cells = cellRemove)
  }
  return(query.object)
}

randomSplit <- function(obj, n=2, seed=44, verbose=F) {
  set.seed(seed)
  lgt <- dim(obj)[2]
  ind <- sample.int(n, lgt, replace = T)
  cell.list <- split(colnames(obj), ind)
  seurat.list <- list()
  if (verbose==TRUE) {
     message(sprintf("Splitting object into %i random subsets", n))
  }
  for (h in 1:n) {
     seurat.list[[h]] <- subset(obj, cells= cell.list[[h]])
  }
  return(seurat.list)
}

#' Load Reference Atlas
#'
#' Load or download the reference map for dataset projection. By the default it downloads the reference atlas of tumour-infiltrating lymphocytes (TILs).
#'
#' @param ref Reference Atlas Seurat object (by default downloads the reference TIL atlas)
#' @examples
#' load.reference.map()
#' @export
load.reference.map <- function(ref="referenceTIL") {
  if(identical(ref,"referenceTIL")){
    print("Loading Default Reference Atlas...")
    refFileName <- paste0(getwd(),"/ref_TILAtlas_mouse_v1.rds")
    refUrl <- "https://ndownloader.figshare.com/files/23136746"
    if (file.exists(refFileName)){
      print(refFileName)
      tryCatch(ref <- readRDS(refFileName), error = function(e){ stop(paste("Reference object",refFileName,"is invalid"))})

    } else {
      print(paste0(refFileName," not found; downloading reference TIL map from the server..."))
      tryCatch(download.file(refUrl, refFileName), error = function(e){ stop("Sorry, it didn't work.")})
      tryCatch(ref <- readRDS(refFileName), error = function(e){ stop(paste("Reference object",refFileName,"is invalid"))})
    }
    tryCatch( print(paste0("Loaded Reference map ",ref@misc$projecTILs)),error = function(e){stop("Invalid Reference object")}   )

  } else {
    print("Loading Custom Reference Atlas...")
    tryCatch(ref <- readRDS(ref), error = function(e){ stop(paste("Reference object",ref,"is invalid"))})
    tryCatch(print(paste0("Loaded Custom Reference map ",ref@misc$projecTILs)),error = function(e){stop("Invalid Reference object")})
  }
  return(ref)
}

#' Project a query scRNA-seq dataset onto a reference atlas
#'
#' This function allows projecting a single-cell RNA-seq dataset onto a reference map of cellular states.
#' It can be useful to interpret a new dataset in the context of an annotated Atlas of known cell types. See `load.reference.map` to load
#' or download a reference atlas.
#'
#' @param query Seurat object with query data
#' @param ref Reference Atlas Seurat object - if NULL, downloads the default TIL reference atlas
#' @param filter.cells Logical value indicating if T cells should be pre-filtered using TILPRED. Default is TRUE. Only set to FALSE if the dataset has been previously checked for non-T cell contaminations.
#' @param query.assay Which assay slot to use for the query
#' @param direct.projection Logical. If true, apply PCA transformation directly without alignment
#' @param seurat.k.filter Integer. For alignment, how many neighbors (k) to use when picking anchors. Default is 200; try lower value in case of failure
#' @param skip.normalize Logical. By default, log-normalize the count data. If you have already normalized your data, you can skip normalization.
#' @return An augmented object \code{query} with projected UMAP coordinates on the reference map and cells classifications
#' @examples
#' data(query_example_seurat)
#' make.projection(query_example_seurat)
#' @export
make.projection <- function(query, ref=NULL, filter.cells=T, query.assay="auto", 
                            direct.projection=FALSE, seurat.k.filter=200, skip.normalize=FALSE) {

  if(is.null(ref)){
    print("Loading Default Reference Atlas...")
    refFileName <- paste0(getwd(),"/ref_TILAtlas_mouse_v1.rds")
    refUrl <- "https://ndownloader.figshare.com/files/23136746"
    if (file.exists(refFileName)){
      print(refFileName)
      tryCatch(ref <- readRDS(refFileName), error = function(e){ stop(paste("Reference object",refFileName,"is invalid"))})

    } else {
      print(paste0(refFileName," not found; I will try to download it and proceed, wish me luck..."))
      tryCatch(download.file(refUrl, refFileName), error = function(e){ stop("Sorry, it didn't work.")})
      tryCatch(ref <- readRDS(refFileName), error = function(e){ stop(paste("Reference object",refFileName,"is invalid"))})
    }
    tryCatch( print(paste0("Loaded Reference map ",ref@misc$projecTILs)),error = function(e){stop("Invalid Reference object")}   )

  }

  pca.dim=dim(ref@misc$umap_object$data)[2]  #use the number of PCs used to build the reference
  projected.list <- list()

  #Reference
  DefaultAssay(ref) <- "integrated"

  ref.var.features <- ref@assays$integrated@var.features

  if(!is.list(query)) {
     query.list <- list(query=query)
  } else {
    query.list <- query
  }
  rm(query)

  for (queryName in names(query.list)){

    #Query
    query <- query.list[[queryName]]
    
    retry.direct <- FALSE

    #Query dataset might be pre-integrated for batch effect correction. If 'integrated' assay exists, use it. Otherwise, use 'RNA'.
    if (!query.assay %in% c("RNA","integrated")){
      if ("integrated" %in% names(query@assays)) {
        query.assay = "integrated"
      } else {
        query.assay = "RNA"
      }
    }
    print(paste0("Using assay ",query.assay," for ",queryName))

    DefaultAssay(query) <- query.assay
    
    if(filter.cells){
      message("Pre-filtering of T cells (TILPRED classifier)...")
      query <- filterCells(query, filter=T)
    } else {
      query <- filterCells(query, filter=F)
    }

    #DefaultAssay(query) <- query.assay
    query <- RenameCells(query, add.cell.id = "Q")
    
    
    #Check if slots are there, and normalize data
    if (!skip.normalize) {
      if (dim(query@assays[[query.assay]]@counts)[1]==0) {
         stop("Counts slot not found in your Seurat object. If you already normalized your data, re-run with option skip.normalize=TRUE")
      }
      query <- NormalizeData(query)  
    } else {
      if (dim(query@assays[[query.assay]]@data)[1]==0) {
        stop("Data slot not found in your Seurat object. Please normalize the data")
      }  
    }
    
    genes4integration <- intersect(ref.var.features, row.names(query))

    if(length(genes4integration)/length(ref.var.features)<0.5){ stop("Too many genes missing. Check input object format") }
    #TODO implement ID mapping? e.g. from ENSEMBLID to symbol?

    if (length(genes4integration)/length(ref.var.features)<0.8) {
      print("Warning! more than 20% of variable genes not found in the query")
    }

    if (direct.projection) {
      projected <- query

      print("DIRECTLY projecting query onto Reference PCA space")
      query.pca.proj <-apply.pca.obj.2(query, pca.obj=ref@misc$pca_object, query.assay=query.assay)
      projected[["pca"]] <- CreateDimReducObject(embeddings = query.pca.proj, key = "PC_", assay = query.assay)

      print("DIRECTLY projecting query onto Reference UMAP space")
      query.umap.proj <- make.umap.predict.2(ref.umap=ref@misc$umap_obj, pca.query.emb = query.pca.proj)
      projected[["umap"]] <- CreateDimReducObject(embeddings = query.umap.proj, key = "UMAP_", assay = query.assay)

      DefaultAssay(projected) <- query.assay
    } else {
      tryCatch(    #Try to do alignment, if it fails (too few cells?) do direct projection
        expr = {

          print(paste0("Aligning ", queryName, " to reference map for batch-correction..."))
          
          if (dim(ref@assays$integrated@scale.data)[2]==0) {
              ref <- ScaleData(ref, do.center=FALSE, do.scale=FALSE, features = genes4integration)
          }
          query <- ScaleData(query, do.center=FALSE, do.scale=FALSE, features = genes4integration)
          
          ref <- RunPCA(ref, features = genes4integration,verbose = F)
          query <- RunPCA(query, features = genes4integration,verbose = F)

          #TODO optimize aligmment for speed? e.g. filter number of anchors STACAS
          proj.anchors <- FindIntegrationAnchors(object.list = c(ref, query), anchor.features = genes4integration,
                                                 dims = 1:pca.dim, k.filter = seurat.k.filter, scale = FALSE, assay=c("integrated",query.assay), reduction = "rpca")
          #Do integration
          all.genes <- intersect(row.names(ref), row.names(query))
          proj.integrated <- IntegrateData(anchorset = proj.anchors, dims = 1:pca.dim, features.to.integrate = all.genes,  preserve.order = T, verbose=F)

          #Subset query data from integrated space
          cells_query<- colnames(query)
          projected <- subset(proj.integrated, cells = cells_query)


          #Make PCA and UMAP projections
          cat("\nProjecting corrected query onto Reference PCA space\n")
          query.pca.proj <-apply.pca.obj.2(projected, pca.obj=ref@misc$pca_object, query.assay="integrated")
          projected[["pca"]] <- CreateDimReducObject(embeddings = query.pca.proj, key = "PC_", assay = "integrated")

          print("Projecting corrected query onto Reference UMAP space")
          query.umap.proj <- make.umap.predict.2(ref.umap=ref@misc$umap_obj, pca.query.emb=query.pca.proj)
          projected[["umap"]] <- CreateDimReducObject(embeddings = query.umap.proj, key = "UMAP_", assay = "integrated")

          DefaultAssay(projected) <- "integrated"
        },
        error = function(e) {
          message(paste("Alignment failed due to:", e, "\n"))
          message("Warning: alignment of query dataset failed - Trying direct projection...")
          retry.direct <<- TRUE
        }
      )
      if (retry.direct) {
        projected <- query

        print("DIRECTLY projecting query onto Reference PCA space")
        query.pca.proj <-apply.pca.obj.2(query, pca.obj=ref@misc$pca_object, query.assay=query.assay)
        projected[["pca"]] <- CreateDimReducObject(embeddings = query.pca.proj, key = "PC_", assay = query.assay)

        print("DIRECTLY projecting query onto Reference UMAP space")
        query.umap.proj <- make.umap.predict.2(ref.umap=ref@misc$umap_obj, pca.query.emb = query.pca.proj)
        projected[["umap"]] <- CreateDimReducObject(embeddings = query.umap.proj, key = "UMAP_", assay = query.assay)

        DefaultAssay(projected) <- query.assay
      }
    }

    projected@assays[[query.assay]]@var.features <- ref.var.features
    projected.list[[queryName]] <- projected
  }
  if(length(projected.list)==1)
     projected.list<-projected.list[[1]]

  return(projected.list)
}

#' Predict cell states of a projected dataset
#'
#' This function uses a nearest-neighbor algorithm to predict a feature (e.g. the cell state) of the query cells. Distances between
#' cells in the reference map and cells in the query are calculated in a reduced space (PCA or UMAP) and the feature is assigned to
#' query cells based on a consensus of its nearest neighbors in the reference object.
#'
#' @param ref Reference Atlas Seurat object
#' @param query Seurat object with query data
#' @param reduction The reduced space used to calculate pairwise distances. One of "pca" or "umap"
#' @param ndim How many dimensions in the reduced space to be used for distance calculations
#' @param k Number of neighbors to assign the cell type
#' @param labels.col The metadata field of the reference to annotate the clusters (default: functional.cluster)
#' @return The query object submitted as parameter, with two additional metadata slots for predicted state and its confidence score
#' @examples
#' cellstate.predict(ref, query_example.seurat)
#' @export
cellstate.predict = function(ref, query, reduction="pca", ndim=10, k=20, labels.col="functional.cluster") {
  require(Seurat)
  tdim <- dim(ref@reductions[[reduction]]@cell.embeddings)[2]
  if (ndim > tdim) {
     warning(sprintf("Number of dimensions ndim=%i is larger than the dimensions in reduction %s - Using only first %i dimensions",ndim,reduction,tdim))
     ndim = tdim
  }
  labels <- ref[[labels.col]][,1]

  ref.space <- ref@reductions[[reduction]]@cell.embeddings[,1:ndim]
  query.space <- query@reductions[[reduction]]@cell.embeddings[,1:ndim]

  pred.type <- rep("Unknown", dim(query.space)[1])
  pred.conf <- numeric(dim(query.space)[1])

  #Use NN-search wrapper in Seurat
  nn.method="rann"
  nn.ranked <- Seurat:::NNHelper(data=ref.space, query=query.space, k = k, method = nn.method)

  for (r in 1:dim(nn.ranked$nn.idx)[1]) {
    top.k <- nn.ranked$nn.idx[r,]
    scores <- sort(table(labels[top.k]), decreasing = T)/k
    pred.type[r] <- names(scores)[1]
    pred.conf[r] <- scores[1]
  }
  pred <- as.data.frame(cbind(row.names(query.space), pred.type, pred.conf), stringsAsFactors = FALSE)
  row.names(pred) <- row.names(query.space)
  colnames(pred) <- c("id","pred.state","confidence")

  pred <- transform(pred, confidence=as.numeric(confidence))

  query <- AddMetaData(query, metadata=pred$pred.state, col.name = labels.col)
  query <- AddMetaData(query, metadata=pred$confidence, col.name = paste0(labels.col,".conf"))
  message(paste0("Creating slots ",labels.col," and ",labels.col, ".conf in query object"))

  return(query)
}

# Function plot.projection
#' Show UMAP projection of query on reference map
#'
#' Plots the UMAP representation of the reference map, together with the projected coordinates of a query dataset.
#'
#' @param ref Reference Atlas Seurat object
#' @param query Seurat object with query data
#' @param labels.col The metadata field to annotate the clusters (default: functional.cluster)
#' @param cols Custom color palette for clusters
#' @return UMAP plot of reference map with projected query set in the same space
#' @examples
#' plot.projection(ref, query_example.seurat)
#' @export
plot.projection = function(ref, query=NULL, labels.col="functional.cluster", cols=NULL) {
  require(Seurat)
  require(ggplot2)

  labels <- ref[[labels.col]][,1]

  states_all <- levels(factor(labels))
  nstates <- length(states_all)
  
  if (!is.null(cols)) {  #custom palette
    if (nstates<=length(cols)) {
      stateColors_func <- cols[1:nstates]
    } else {  
      warning("Not enough colors provided. Making an automatic palette")
      stateColors_func <- rainbow(n=nstates)
    }  
  } else {   #default palette
    stateColors_func <- c("#edbe2a","#A58AFF","#53B400","#F8766D","#00B6EB","#d1cfcc","#FF0000","#87f6a5","#e812dd")
    
    if (nstates<=length(stateColors_func)) {
      stateColors_func <- stateColors_func[1:nstates]
    } else {   #make a new palette
      stateColors_func <- rainbow(n=nstates)
    }
  }
  names(stateColors_func) <- states_all
  cols_use <- stateColors_func[states_all]

  if (is.null(query)) {
    p <- DimPlot(ref, reduction="umap", label = F, group.by = labels.col, repel = T, cols=cols_use) +
      ggtitle ("Reference map") + theme(aspect.ratio=1)
  } else {
    p <- DimPlot(ref, reduction="umap", label = F, group.by = labels.col, repel = T, cols=cols_use) +
      geom_point(data.frame(query@reductions$umap@cell.embeddings), mapping=aes(x=UMAP_1,y=UMAP_2),alpha=0.6, size=1,shape=17, color="gray10") +
      geom_density_2d(data=data.frame(query@reductions$umap@cell.embeddings), mapping=aes(x=UMAP_1,y=UMAP_2),color="black",n=200,h=2) +
      ggtitle ("Projection of query on reference map") + theme(aspect.ratio=1)
  }
  return(p)
}

# Function plot.statepred.composition
#' Summarize the predicted cell states of an object
#'
#' Makes a barplot of the frequency of cell states in a query object.
#'
#' @param query Seurat object with query data
#' @param labels.col The metadata field used to annotate the clusters (default: functional.cluster)
#' @param metric One of `Count` or `Percent`. `Count` plots the absolute number of cells, `Percent` the fraction on the total number of cells.
#' @param cols Custom color palette for clusters
#' @return Barplot of predicted state composition
#' @examples
#' plot.statepred.composition(query_example.seurat)
#' @export
plot.statepred.composition = function(ref, query, labels.col="functional.cluster",cols=NULL, metric=c("Count","Percent")) {
  require(reshape2)
  require(ggplot2)
  
  metric <- tolower(metric[1])
  
  labels <- ref[[labels.col]][,1]

  states_all <- levels(factor(labels))
  nstates <- length(states_all)
  
  if (!is.null(cols)) {  #custom palette
    if (nstates<=length(cols)) {
      stateColors_func <- cols[1:nstates]
    } else {  
      warning("Not enough colors provided. Making an automatic palette")
      stateColors_func <- rainbow(n=nstates)
    }  
  } else {   #default palette
    stateColors_func <- c("#edbe2a","#A58AFF","#53B400","#F8766D","#00B6EB","#d1cfcc","#FF0000","#87f6a5","#e812dd")
    if (nstates<=length(stateColors_func)) {
      stateColors_func <- stateColors_func[1:nstates]
    } else {   #make a new palette
      stateColors_func <- rainbow(n=nstates)
    }
  }
  names(stateColors_func) <- states_all
  cols_use <- stateColors_func[states_all]

  tb <- table(factor(query[[labels.col]][,1], levels=states_all))
  
  if (metric=="percent") {  #normalize
    tb <- tb*100/sum(tb)
    tb.m <- melt(tb)
    colnames(tb.m) <- c("Cell_state","Perc_cells")
    p <- ggplot(tb.m, aes(x=Cell_state, y=Perc_cells, fill=Cell_state)) + geom_bar(stat="identity") +
      scale_fill_manual(values=cols_use) +
      theme(axis.text.x=element_blank(), legend.position="left")
  } else if (metric=="count") {
    tb.m <- melt(tb)
    colnames(tb.m) <- c("Cell_state","Ncells")
    p <- ggplot(tb.m, aes(x=Cell_state, y=Ncells, fill=Cell_state)) + geom_bar(stat="identity") +
      scale_fill_manual(values=cols_use) +
      theme(axis.text.x=element_blank(), legend.position="left")
  } else {
    stop("Unknown metric specified (Must be either Count or Percent")
  }
  
  return(p)
}

# Function plot.states.radar
#' Show expression level of key genes
#'
#' Makes a radar plot of the expression level of a set of genes. It can be useful to compare the gene expression profile of different cell
#' states in the reference atlas vs. a projected set.
#'
#' @param ref Seurat object with reference atlas
#' @param query Query data, either as a Seurat object or as a list of Seurat objects
#' @param labels.col The metadata field used to annotate the clusters (default: functional.cluster)
#' @param genes4radar Which genes to use for plotting (default: c("Foxp3","Cd4","Cd8a","Tcf7","Ccr7","Gzmb","Gzmk","Pdcd1","Havcr2","Tox,"Mki67")
#' @param min.cells Only display cell states with a minimum number of cells
#' @param cols Custom color palette for clusters
#' @param return Return the ggplot object instead of printing it to the default device
#' @return Radar plot of gene expression of key genes by cell subtype
#' @examples
#' plot.states.radar(ref)
#' @export
plot.states.radar = function(ref, query=NULL, labels.col="functional.cluster",
                             genes4radar=NULL, min.cells=10, cols=NULL, return=F) {
  require(ggplot2)
  require(gridExtra)

  if (is.null(genes4radar)) {
     genes4radar <- c("Foxp3","Cd4","Cd8a","Tcf7","Ccr7","Gzmb","Gzmk","Pdcd1","Havcr2","Tox","Mki67")
  }
  genes4radar <- sort(genes4radar)
  order <- match(genes4radar, row.names(ref@assays$RNA@data))
  rr <- ref@assays$RNA@data[order,]

  labels <- ref[[labels.col]][,1]
  states_all <- levels(factor(labels))
  nstates <- length(states_all)

  if (!is.null(query)) {
    if (!is.list(query)) {
       query <- list("Query" = query)
    }
    if (is.null(names(query))) {
      for (i in 1:length(query)) {
         names(query)[[i]] <- paste0("Query",i)
      }
    }
    labels.q <- list()
    qq <- list()

    for (i in 1:length(query)) {
       labels.q[[i]] <- query[[i]][[labels.col]][,1]
       order <- match(genes4radar, row.names(query[[i]]@assays$RNA@data))
       qq[[i]] <- query[[i]]@assays$RNA@data[order,]
    }
  }

  if (!is.null(cols)) {  #custom palette
    if (nstates<=length(cols)) {
      stateColors_func <- cols[1:nstates]
    } else {  
       warning("Not enough colors provided. Making an automatic palette")
       stateColors_func <- rainbow(n=nstates)
    }  
  } else {   #default palette
    stateColors_func <- c("#edbe2a","#A58AFF","#53B400","#F8766D","#00B6EB","#d1cfcc","#FF0000","#87f6a5","#e812dd")
    if (nstates<=length(stateColors_func)) {
      stateColors_func <- stateColors_func[1:nstates]
    } else {   #make a new palette
      stateColors_func <- rainbow(n=nstates)
    }
  }
  names(stateColors_func) <- states_all
  cols_use <- stateColors_func[states_all]

  #Get raw expression means, to normalize by gene
  m <- matrix(, nrow = length(states_all), ncol = length(genes4radar))
  rownames(m) <- states_all
  colnames(m) <- genes4radar
  for (i in 1:length(states_all)) {
    s <- states_all[i]
    m[i,] <- apply(rr[, labels == s], MARGIN=1, mean)
  }
  normfacs <- apply(m, MARGIN=2, function(x) {max(c(1,x))})

  pll <- list()
  for (j in 1:length(states_all)) {
    s <- states_all[j]
    fill_use <- c(rep(cols_use[j], length(query)+1))
    names(fill_use) <- c("Reference",names(query))

    col_use <- c("black","orange")
    if (length(query)>1) { col_use <- c(col_use, rainbow(length(query)-1)) }
    names(col_use) <-  c("Reference",names(query))

    this.mean <- apply(rr[, labels == s], MARGIN=1, mean)
    this.mean <- this.mean/normfacs

    this.df <- data.frame(t(rbind(names(this.mean), this.mean, "Reference")))
    colnames(this.df) <- c("Gene","Expression","Dataset")
    this.df$Expression <- as.numeric(as.character(this.df$Expression))

    i <- 1
    while (i <= length(query)) {
      m <- as.matrix(qq[[i]][, labels.q[[i]] == s])
      if (dim(m)[2] >= min.cells) {
        q.mean <- apply(m, MARGIN=1, mean)
        q.mean <- q.mean/normfacs
        q.df <- data.frame(t(rbind(names(q.mean), q.mean, names(query)[[i]])))
        colnames(q.df) <- c("Gene","Expression","Dataset")
        q.df$Expression <- as.numeric(as.character(q.df$Expression))
        this.df <- rbind(this.df, q.df)
      }
      i=i+1
    }

    ymin <- min(c(-0.1, min(this.df$Expression)))
    ymax <- max(c(1, max(this.df$Expression)))


    pll[[j]] <- ggplot(data=this.df,  aes(x=Gene, y=Expression, group= Dataset, colour=Dataset, fill=Dataset)) +
      geom_point(size=2) +
      geom_polygon(size = 1, alpha= 0.2) +
      ylim(ymin, ymax) + ggtitle(s)  +
      scale_x_discrete() +
      scale_fill_manual(values= fill_use) +
      scale_colour_manual(values= col_use) +
      theme_light()+
      coord_polar()

  }
  g <- do.call("arrangeGrob", c(pll, ncol=3, top=paste0("Radar plots for ", labels.col)))

  if(return){
    return(g)
  } else{
    return(plot(g))
  }

}

# Function find.discriminant.dimensions
#' Find discriminant dimensions
#'
#' Searches PCA or ICA dimensions where the query set deviates the most from a control set or from the reference map. It can
#' be useful to suggest novel cell states that escape from the main axes of diversity of the UMAP
#'
#' @param ref Seurat object with reference atlas
#' @param query Seurat object with query data
#' @param query.control Optionally, you can compare your query with a control sample, instead of the reference
#' @param query.assay The data slot to be used for enrichment analysis
#' @param state Perform discriminant analysis on this cell state. Can be either:
#' \itemize{
#'   \item{"largest" - Performs analysis on the cell state most represented in the query set(s)}
#'   \item{"all" - Performs analysis on the complete dataset, using all cells}
#'   \item{A specific cell state, one of the states in metadata field labels.col}
#' }
#' @param labels.col The metadata field used to annotate the clusters (default: functional.cluster)
#' @param reduction Which dimensionality reduction to use (either ICA or PCA)
#' @param test Which test to perform between the dataset distributions in each ICA/PCA dimension. One of `ks` (Kolmogorov-Smirnov) or `t.test` (T-test)
#' @param ndim How many dimensions to consider in the reduced ICA/PCA space
#' @param print.n The nunmber of top dimensions to return to STDOUT
#' @return A list of PCA or ICA dimensions, ranked by the enrichment of the query vs. control state (if provided), otherwise of the query vs. the reference map
#' @examples
#' find.discriminant.dimensions(ref, query=query.set)
#' find.discriminant.dimensions(ref, query=query.set, query.control=control.set)
#' @export
find.discriminant.dimensions <- function(ref, query, query.control=NULL, query.assay="RNA",
                                         state="largest", labels.col="functional.cluster",
                                         reduction="ICA", test=c("ks","t.test"), ndim=50, print.n=3) {

  reduction=tolower(reduction)
  test=test[1]

  if (is.null(ref)) {stop("Please provide the reference object (ref")}
  if (is.null(query)) {stop("Please provide a query object (query)")}
  #Determine cell state for analysis
  if (is.null(state) | state=="largest") {
    if (!is.null(query.control)) {
      ss <- table(rbind(query[[labels.col]], query.control[[labels.col]]))
    } else {
      ss <- table(query[[labels.col]])
    }
    state <- names(sort(ss, decreasing = T))[1]
    message(paste0("Performing discriminant analysis using the most abundant cell state in the query - ", state))
  } else if (state=="all") {
    message("Performing discriminant analysis using all cells")
  } else {
    if (!state %in% query[[labels.col]][,1]) {
      stop(sprintf("State %s not found in query metadata colum %s", state, labels.col))
    }
    if (!is.null(query.control) & !state %in% query.control[[labels.col]][,1]) {
      stop(sprintf("State %s not found in query.control metadata colum %s", state, labels.col))
    }
    message(paste0("Performing discriminant analysis with user-specified state - ", state))
  }
  #Subset query data on specific state
  if (state=="all") {
     ref.cells=seq(1, dim(ref)[2])
     query.cells=seq(1, dim(query)[2])
     if (!is.null(query.control)) {
        query.c.cells=seq(1, dim(query.control)[2])
     }
  } else {
     ref.cells=which(ref[[labels.col]]==state)
     query.cells=which(query[[labels.col]]==state)
     if (!is.null(query.control)) {
       query.c.cells=which(query.control[[labels.col]]==state)
     }
  }

  if (reduction=="ica") {
    if (is.null(ref@reductions[[reduction]])) {
      message("Reduction ICA not found. Calculating ICA for reference object")
      ref <- run.ica(ref, ndim=ndim)
    }
    ref_dimRed <- ref@misc$ica
    perturb_dimRed <- apply.ica.obj(query=query, query.assay=query.assay, ica.obj=ref_dimRed)
    if (!is.null(query.control)) {
      control_dimRed <- apply.ica.obj(query=query.control, query.assay=query.assay, ica.obj=ref_dimRed)
    }
  } else if (reduction=="pca") {
    ref_dimRed <- ref@misc$pca_obj
    perturb_dimRed <- apply.pca.obj.2(query=query, query.assay=query.assay, pca.obj=ref_dimRed)
    if (!is.null(query.control)) {
      control_dimRed <- apply.pca.obj.2(query=query.control, query.assay=query.assay, pca.obj=ref_dimRed)
    }
  } else {
    stop(paste0("Unrecognized reduction slot: ", reduction))
  }

  ndim <- min(ndim, length(colnames(perturb_dimRed)))
  w <- rep(1,ndim)
  names(w) <- colnames(perturb_dimRed)[1:ndim]
  stat <- rep(1,ndim)
  names(stat) <- colnames(perturb_dimRed)[1:ndim]
  pval <- rep(1,ndim)
  names(pval) <- colnames(perturb_dimRed)[1:ndim]

  if (!is.null(query.control)) {
     message("Query and control datasets was provided. Determining discriminant components of Query vs. Control...")
  } else {
     message("Single query dataset was provided. Determining discriminant components of Query vs. Reference...")
  }

  for (pc in 1:ndim) {
    d1 <- perturb_dimRed[query.cells,pc]
    if (!is.null(query.control)) {
       d2 <- control_dimRed[query.c.cells,pc]
    } else {
       d2 <- ref@reductions[[reduction]]@cell.embeddings[ref.cells, pc]
    }
    if (test=="ks") {
       this.test <- ks.test(d1, d2, alternative="two.sided")
    } else {
      this.test <- t.test(d1, d2)
    }

    w[pc] <- mean(d1) - mean(d2)
    stat[pc] <- this.test$statistic
    pval[pc] <- this.test$p.value * ndim   #multiple testing
  }
  w.sorted <- sort(abs(stat), decreasing = T)

  for (i in 1:print.n) {
    topPC <- names(w.sorted)[i]
    pc.index <- match(topPC, names(stat))
    feats <- ref@reductions[[reduction]]@feature.loadings[,pc.index]
    topgenes <- names(head(sort(abs(feats), decreasing = T), 10))
    topgenes.sign <- feats[topgenes]
    if (w[pc.index]>0) {
       topgenes.p <- topgenes.sign[topgenes.sign>=0]
       topgenes.n <- topgenes.sign[topgenes.sign<0]
    } else {
      topgenes.p <- topgenes.sign[topgenes.sign<0]
      topgenes.n <- topgenes.sign[topgenes.sign>=0]
    }

    pval2print <- ifelse(pval[pc.index]<0.001, sprintf("%.1e",pval[pc.index]), sprintf("%.4f",pval[pc.index]))
    message(sprintf("-----------------------\nTop %s component %i: %s Stat %.3f p.val %s",
                    reduction, i, topPC, stat[pc.index], pval2print))
    message("Driver genes for this component:")
    message(paste(c("Higher in query +++ ", names(topgenes.p)), collapse = " "))
    message(paste(c("Lower in query  --- ", names(topgenes.n)), collapse = " "))
  }

  return(w.sorted)
}


# Function plot.discriminant.3d
#' 3D plot of reference map with extra discriminant dimension
#'
#' Add an extra dimension to the reference map (it can be suggested by `find.discriminant.dimensions`), to explore additional axes of variability
#' in a query dataset compared to the reference map.
#'
#' @param ref Seurat object with reference atlas
#' @param query Seurat object with query data
#' @param query.control Optionally, you can compare your query with a control sample, instead of the reference
#' @param query.assay The data slot to be used for enrichment analysis
#' @param labels.col The metadata field used to annotate the clusters (default: functional.cluster)
#' @param reduction Which dimensionality reduction to use (either ICA or PCA)
#' @param query.state Only plot the query cells from this specific state
#' @param extra.dim The additional dimension to be added on the z-axis of the plot. Can be either:
#' \itemize{
#'   \item{An ICA or PCA dimension (e.g. ICA_10). See `find.discriminant.dimensions`}
#'   \item{"cycling.score" - The enrichment score for the cycling signature calculated by TILPRED}
#' }
#' @return A three dimensional plot with UMAP_1 and UMAP_2 on the x and y axis respectively, and the specified `extra.dim` on the z-axis.
#' @examples
#' plot.discriminant.3d(ref, query=query, extra.dim="ICA_19")
#' plot.discriminant.3d(ref, query=treated.set, query.control=control.set, extra.dim="ICA_2")
#' @export

plot.discriminant.3d <- function(ref, query, query.control=NULL, reduction="ICA", query.assay="RNA",
                             labels.col="functional.cluster", extra.dim="ICA_1", query.state=NULL) {
  require(plotly)

  reduction=tolower(reduction)
  message(paste0("Generating UMAP with 3rd dimension on ",extra.dim))

  ref.3d <- ref
  ref.3d$queryGroup <- "Reference"

  #Only show cells of a specific state for the query
  if(!is.null(query.state)) {
    query.cells=colnames(query)[which(query[[labels.col]]==query.state)]
    query <- subset(query, cells=query.cells)
    if (!is.null(query.control)) {
      query.c.cells=colnames(query.control)[which(query.control[[labels.col]]==query.state)]
      query.control <- subset(query.control, cells=query.c.cells)
    }
  }

  #Prepare embeddings
  if (!is.null(query.control)) {
    q.labs <- c(rep("Control", dim(query.control)[2]), rep("Query", dim(query)[2]))
    q.umaps <- rbind(query.control@reductions$umap@cell.embeddings, query@reductions$umap@cell.embeddings)
    query <- merge(query.control, query)
    query$queryGroup <- q.labs
    query[["umap"]] <- CreateDimReducObject(embeddings = q.umaps, key = "UMAP_", assay = query.assay)
  } else {
     query$queryGroup <- "Query"
  }
  query@meta.data[,labels.col] <- "Query"

  ref.3d@meta.data <- rbind(ref.3d@meta.data[,c(labels.col,"queryGroup","cycling.score")],
                            query@meta.data[,c(labels.col,"queryGroup","cycling.score")])

  ref.3d@reductions$umap@cell.embeddings <- rbind(ref.3d@reductions$umap@cell.embeddings,
                                                            query@reductions$umap@cell.embeddings)

  if (reduction=="ica"){
    projected <- apply.ica.obj(query=query, query.assay=query.assay, ica.obj=ref.3d@misc$ica)
  } else if (reduction=="pca") {
    projected <- apply.pca.obj.2(query=query, query.assay=query.assay, pca.obj=ref.3d@misc$pca_object)
  } else {
    stop(paste0("Unrecognized reduction slot: ", reduction))
  }
  ref.3d@reductions[[reduction]]@cell.embeddings <- rbind(ref.3d@reductions[[reduction]]@cell.embeddings, projected)

  ref.3d@meta.data <- cbind(ref.3d@meta.data, ref.3d@reductions$umap@cell.embeddings)

  if (extra.dim=="cycling" | extra.dim=="cycling.score") {
    ref.3d@meta.data <- cbind(ref.3d@meta.data,Discriminant=ref.3d$cycling.score)
  } else {
     ref.3d@meta.data <- cbind(ref.3d@meta.data,Discriminant=ref.3d@reductions[[reduction]]@cell.embeddings[,extra.dim])
  }

  ref.3d@meta.data$Group <- factor(ifelse(ref.3d@meta.data$queryGroup=="Reference",
                                          as.character(ref.3d@meta.data[labels.col]),
                                          ref.3d@meta.data$queryGroup))

  if (is.null(query.control)) {
    cols <- c(Reference="gray50",Query="red")
  } else {
    cols <- c(Reference="gray50",Control="green",Query="red")
  }

  plotting.data <- ref.3d@meta.data[,c("UMAP_1", "UMAP_2","Discriminant", labels.col,"queryGroup","Group")]
  plotting.data$size <- ifelse(plotting.data$queryGroup == "Reference",0.3,6)

  g <- plot_ly(data = plotting.data,
          x = ~UMAP_1, y = ~UMAP_2, z = ~Discriminant,
          color = ~queryGroup,
          type = "scatter3d",
          mode = "markers",
          text=~queryGroup,
          hoverinfo="text",
          alpha=0.6,
          alpha_stroke=0.6,
          size=~size,
          colors=cols
  )

  g <- g %>% layout(
    title = paste0("Projection of query on reference map + dimension ",extra.dim),
    scene = list(
      xaxis = list(title = "UMAP_1"),
      yaxis = list(title = "UMAP_2"),
      zaxis = list(title = extra.dim)
    ))
  print(g)
  return(g)
}

# Function convert count matrix to Seurat object


