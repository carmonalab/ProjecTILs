#' Load Reference Atlas
#'
#' Load or download the reference map for dataset projection. By the default it downloads the reference atlas of tumour-infiltrating lymphocytes (TILs).
#'
#' @param ref Reference atlas as a Seurat object (by default downloads a mouse reference TIL atlas).
#'     To use a custom reference atlas, provide a .rds object or a URL to a .rds object, storing a Seurat object
#'     prepared using \link{make.reference}
#' @examples
#' load.reference.map()
#' @export load.reference.map
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
      tryCatch(download.file(refUrl, refFileName), error = function(e){ stop("Sorry, download failed.")})
      tryCatch(ref <- readRDS(refFileName), error = function(e){ stop(paste("Reference object",refFileName,"is invalid"))})
    }
    tryCatch( print(paste0("Loaded Reference map ",ref@misc$projecTILs)),error = function(e){stop("Invalid Reference object")}   )

  } else {
    if (grepl("^[ftp|http]", ref, perl = T)) {
      refUrl <- ref
      refFileName <- paste0(getwd(),"/custom_reference.rds")
      print(sprintf("Trying to download custom reference from %s...", refUrl))
      
      tryCatch(download.file(refUrl, refFileName), error = function(e){ stop("Sorry, download failed.")}) 
    } else if (file.exists(ref)) {
      refFileName <- ref
    } else {
      stop("Provide ref is not a valid reference or a valid URL.")
    }
    print("Loading Custom Reference Atlas...")
    tryCatch(ref <- readRDS(refFileName), error = function(e){ stop(paste("Reference object",ref,"is invalid"))})
    tryCatch(print(paste0("Loaded Custom Reference map ",ref@misc$projecTILs)),error = function(e){stop("Invalid Reference object")})
  }
  return(ref)
}

#' Read to memory a query expression matrix
#'
#' Load a query expression matrix to be projected onto the reference atlas. Several formats (10x, hdf5, raw and log counts) 
#' are supported - see \code{type} parameter for details
#'
#' @param filename Path to expression matrix file or folder
#' @param type Expression matrix format (10x, hdf5, raw, raw.log2)
#' @param project.name Title for the project
#' @param min.cells Only keep genes represented in at least min.cells number of cells
#' @param min.features Only keep cells expressing at least min.features genes
#' @param gene.column.10x For 10x format - which column of genes.tsv or features.tsv to use for gene names
#' @param raw.rownames For raw matrix format - A vector of row names, or a single number giving the column of the table which contains the row names
#' @param raw.sep For raw matrix format - Separator for raw expression matrix
#' @param raw.header For raw matrix format - Use headers in expression matrix
#' @param use.readmtx Use ReadMtx function to read in 10x files with custom names
#' @return A Seurat object populated with raw counts and normalized counts for single-cell expression
#' @examples
#' fname <- "./sample_data"
#' querydata <- read.sc.query(fname, type="10x")
#' @export read.sc.query

read.sc.query <- function(filename,
                          type=c("10x","hdf5","raw","raw.log2"),
                          project.name="Query",
                          min.cells = 3,
                          min.features = 50,
                          gene.column.10x=2,
                          raw.rownames=1,
                          raw.sep=c("auto"," ","\t",","),
                          raw.header=TRUE,
                          use.readmtx=TRUE) {
  
  if (is.null(filename)) {stop("Please provide a query dataset in one of the supported formats")}
  type = tolower(type[1])
  
  if (type == "10x") {
    fl <- list.files(filename)
    matrix.file <- grep("matrix.mtx", fl, value=TRUE)[1]
    feature.file <- grep("features.tsv|genes.tsv", fl, value=TRUE)[1]
    barcode.file <- grep("barcodes.tsv", fl, value=TRUE)[1]
    
    if (is.na(matrix.file)) stop("Cannot find matrix file")
    if (is.na(feature.file)) stop("Cannot find genes file")
    if (is.na(barcode.file)) stop("Cannot find barcode file")
    
    matrix.file <- sprintf("%s/%s", filename, matrix.file)
    feature.file <- sprintf("%s/%s", filename, feature.file)
    barcode.file <- sprintf("%s/%s", filename, barcode.file)

    if (use.readmtx) {
      query.exp <- ReadMtx.fix(mtx=matrix.file,
                               cells=barcode.file,
                               features=feature.file,
                               feature.column=gene.column.10x)
    } else {
      query.exp <- Read10X(filename, gene.column = gene.column.10x)
    }
  } else if (type == "hdf5") {
    query.exp <- Read10X_h5(filename)
  } else if (type == "raw" | type == "raw.log2") {
    
    raw.sep <- raw.sep[1]
    if (raw.sep == "auto") {
      raw.sep <- guess_raw_separator(f=filename)
      if (is.null(raw.sep)) {
        stop("Could not guess separator for raw matrix format. Try specifying manually with raw.sep parameter")
      }
    }
    
    p <- regexpr("\\.([[:alnum:]]+)$", filename)
    extension <- ifelse(p > -1L, substring(filename, p + 1L), "")
    if (extension == "gz") {
      query.exp <- read.table(gzfile(filename),row.names=raw.rownames, sep=raw.sep, header=raw.header)
    } else {
      query.exp <- read.table(filename,row.names=raw.rownames, sep=raw.sep, header=raw.header)
    }
    query.exp[is.na(query.exp)] <- 0
    if (type == "raw.log2") {
      query.exp <- 2^(query.exp)-1
    }
    
    #Also try to determine whether genes are on rows or columns
    data(Hs2Mm.convert.table)
    gnames <- c(Hs2Mm.convert.table$Gene.MM, Hs2Mm.convert.table$Gene.stable.ID.HS, Hs2Mm.convert.table$Gene.HS)
    gr <- length(intersect(rownames(query.exp), gnames))
    gc <- length(intersect(colnames(query.exp), gnames))
    gmax <- max(gr, gc)
    if (gmax==0) {
      stop("Could not find gene names in matrix. Check matrix format")
    }
    if (gc>gr) {  #flip rows and columns
      query.exp <- t(query.exp)
    }
    
  } else {
     stop("Please provide a query dataset in one of the supported formats")
  } 
  query.seurat <- CreateSeuratObject(counts=query.exp, project=project.name, min.cells=min.cells, min.features=min.features)
  query.seurat<- NormalizeData(query.seurat)
  return(query.seurat)
}

#' Project a query scRNA-seq dataset onto a reference atlas
#'
#' This function allows projecting ("query") single-cell RNA-seq datasets onto a reference map
#' (i.e. a curated and annotated scRNA-seq dataset). 
#' To project multiple datasets, submit a list of Seurat objects with the query parameter.
#' The projection consists of 3 steps: 
#' \itemize{
#'  \item{pre-processing: optional steps which might include pre-filtering of cells by markers using `scGate`,
#' data normalization, and ortholog conversion.}
#'  \item{batch-effect correction: uses built-in STACAS algorithm to detect and correct for batch effects
#' (this step assumes that at least a fraction of the cells in the query are in the same state than cells in
#' the reference)}
#'  \item{embedding of corrected query data in the reduced-dimensionality spaces (PCA and UMAP) of the reference map.}
#' }
#' 
#' See \link{load.reference.map} to load or download a reference atlas. See
#' also \link{ProjecTILs.classifier} to use ProjecTILs as a cell type classifier.
#'
#' @param query Query data, either as single Seurat object or as a list of Seurat object
#' @param ref Reference Atlas - if NULL, downloads the default TIL reference atlas
#' @param query.assay Which assay slot to use for the query (defaults to DefaultAssay(query))
#' @param direct.projection If true, apply PCA transformation directly without alignment
#' @param STACAS.anchor.coverage Focus on few robust anchors (low STACAS.anchor.coverage) or on a large amount
#'     of anchors (high STACAS.anchor.coverage). Must be number between 0 and 1.
#' @param STACAS.correction.scale Slope of sigmoid function used to determine strength of batch effect correction.
#' @param STACAS.k.anchor Integer. For alignment, how many neighbors (k) to use when picking anchors.
#' @param STACAS.k.weight Number of neighbors to consider when weighting anchors.
#'     Default is "max", which disables local anchor weighting.
#' @param skip.normalize By default, log-normalize the count data.
#'     If you have already normalized your data, you can skip normalization.
#' @param filter.cells Pre-filter cells using `scGate`. Only set to FALSE if the dataset has 
#'     been previously subset to cell types represented in the reference.
#' @param scGate_model scGate model used to filter target cell type from query data
#'     (if NULL use the model stored in \code{ref@@misc$scGate})
#' @param ortholog_table Dataframe for conversion between ortholog genes
#'     (by default package object \code{Hs2Mm.convert.table})
#' @param fast.umap.predict Fast approximation for UMAP projection. Uses coordinates of nearest neighbors in 
#'     PCA space to assign UMAP coordinates (credits to Changsheng Li for the implementation)
#' @param ncores Number of cores for parallel execution (requires \link{BiocParallel})
#' @param progressbar Whether to show a progress bar for projection process or not (requires \link{BiocParallel})
#' @return An augmented Seurat object with projected UMAP coordinates on the reference map
#' @examples
#' data(query_example_seurat)
#' ref <- load.reference.map()
#' make.projection(query_example_seurat, ref=ref)
#' @import Seurat
#' @importFrom STACAS FindAnchors.STACAS IntegrateData.STACAS
#' @importFrom BiocParallel MulticoreParam bplapply
#' @importFrom stats aggregate quantile sd
#' @export make.projection
make.projection <- function(query, ref=NULL,
                            filter.cells=TRUE,
                            query.assay=NULL,
                            direct.projection=FALSE,
                            STACAS.anchor.coverage=0.7,
                            STACAS.correction.scale=100,
                            STACAS.k.anchor=5,
                            STACAS.k.weight="max",
                            skip.normalize=FALSE,
                            fast.umap.predict=FALSE,
                            ortholog_table=NULL,
                            scGate_model=NULL,
                            ncores=1,
                            progressbar = TRUE) {
   
  
  if(is.null(ref)){
    print("Loading Default Reference Atlas...")
    refFileName <- paste0(getwd(),"/ref_TILAtlas_mouse_v1.rds")
    refUrl <- "https://ndownloader.figshare.com/files/23136746"
    if (file.exists(refFileName)){
      print(refFileName)
      tryCatch(ref <- readRDS(refFileName),
               error = function(e){ stop(paste("Reference object",refFileName,"is invalid"))})

    } else {
      print(paste0(refFileName," not found; I will try to download it and proceed, wish me luck..."))
      tryCatch(download.file(refUrl, refFileName),
               error = function(e){ stop("Sorry, it didn't work.")})
      tryCatch(ref <- readRDS(refFileName),
               error = function(e){ stop(paste("Reference object",refFileName,"is invalid"))})
    }
    tryCatch( print(paste0("Loaded Reference map ",ref@misc$projecTILs)),
              error = function(e){stop("Invalid Reference object")}   )

  }
  projected.list <- list()
  if (is.null(ortholog_table)) {
     data(Hs2Mm.convert.table)
     ortholog_table <- Hs2Mm.convert.table 
  } 
  
  if(!is.list(query)) {
     query.list <- list(query=query)
  } else {
    query.list <- query
    if (is.null(names(query.list))) {
       names(query.list) <- paste0("query",c(1:length(query.list)))
    }
  }
  rm(query)
  
  #Parallelize (ncores>1)
  if (ncores > length(query.list)) {
    ncores <- length(query.list)
  }
  param <- BiocParallel::MulticoreParam(workers=ncores, progressbar = progressbar)
  
  #Projection over list of datasets
  projected.list <- BiocParallel::bplapply(
    X = 1:length(query.list), 
    BPPARAM =  param,
    FUN = function(i) {
         projection.helper(query=query.list[[i]], ref=ref,
                           filter.cells=filter.cells,
                           query.assay=query.assay,
                           direct.projection=direct.projection,
                           fast.umap.predict=fast.umap.predict,
                           STACAS.k.anchor=STACAS.k.anchor,
                           STACAS.k.weight=STACAS.k.weight,
                           STACAS.anchor.coverage=STACAS.anchor.coverage,
                           STACAS.correction.scale=STACAS.correction.scale,
                           remove.thr=0,
                           alpha=0.5,
                           ncores=ncores,
                           ortholog_table=ortholog_table,
                           skip.normalize=skip.normalize,
                           id=names(query.list)[i],
                           scGate_model=scGate_model)
      }
  )
      
  names(projected.list) <- names(query.list)

  #De-list if single object was submitted
  if(length(projected.list)==1)
     projected.list <- projected.list[[1]]

  return(projected.list)
}

#' Predict cell states of a projected dataset
#'
#' This function uses a nearest-neighbor algorithm to predict a feature (e.g. the cell state) of the query cells. Distances between
#' cells in the reference map and cells in the query are calculated in a reduced space (PCA or UMAP) and the feature is assigned to
#' query cells based on a consensus of its nearest neighbors in the reference object.
#'
#' @param ref Reference Atlas
#' @param query Seurat object with query data
#' @param reduction The dimensionality reduction used to calculate pairwise distances. One of "pca" or "umap"
#' @param ndim How many dimensions in the reduced space to be used for distance calculations
#' @param k Number of neighbors to assign the cell type
#' @param labels.col The metadata field of the reference to annotate the clusters (default: functional.cluster)
#' @return The query object submitted as parameter, with two additional metadata slots for predicted state and its confidence score
#' @examples
#' data(query_example_seurat)
#' ref <- load.reference.map()
#' q <- make.projection(query_example_seurat, ref=ref)
#' q <- cellstate.predict(ref, query=q)
#' table(q$functional.cluster)
#' @import Seurat
#' @export cellstate.predict
cellstate.predict = function(ref, query,
                             reduction="pca",
                             ndim=NULL,
                             k=20,
                             labels.col="functional.cluster") {
  
  if (is.null(ndim)) {
    if (!is.null(ref@misc$umap_object$data)) {
      ndim <- ncol(ref@misc$umap_object$data)
    } else {
      stop("Please specify ndim parameter.")
    }
  }  
  
  pca.dim=dim(ref@misc$umap_object$data)[2]
  
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

  for (r in 1:dim(nn.ranked@nn.idx)[1]) {
    top.k <- nn.ranked@nn.idx[r,]
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
#' @param ref Reference Atlas
#' @param query Seurat object with query data
#' @param labels.col The metadata field to annotate the clusters (default: functional.cluster)
#' @param cols Custom color palette for clusters
#' @param linesize Contour line thickness for projected query
#' @param pointsize Point size for cells in projected query
#' @param ref.alpha Transparency parameter for reference cells
#' @param ref.size Adjust point size for reference cells
#' @param ... Additional parameters for \code{DimPlot}, e.g. raster=T to
#'    limit image size
#' @return UMAP plot of reference map with projected query set in the same space
#' @examples
#' data(query_example_seurat)
#' ref <- load.reference.map()
#' q <- Run.ProjecTILs(query_example_seurat, ref=ref, fast.umap.predict=TRUE)
#' plot.projection(ref=ref, query=q)
#' @import Seurat
#' @import ggplot2
#' @importFrom scales alpha
#' @importFrom grDevices rainbow
#' @export plot.projection

plot.projection = function(ref, query=NULL, labels.col="functional.cluster",
                          cols=NULL, linesize=1, pointsize=1,
                          ref.alpha=0.3, ref.size=NULL, ...) {
  
  labels <- ref[[labels.col]][,1]
  
  states_all <- levels(factor(labels))
  nstates <- length(states_all)
  
  if (!is.null(cols)) {  #custom palette
    if (nstates>length(cols)) {
      warning("Not enough colors provided. Making an automatic palette")
      palette <- rainbow(n=nstates)
    } else {
      palette <- cols
    } 
  } else {   #default palette
    if (!is.null(ref@misc$atlas.palette)) {  #read directly from atlas, if stored
      palette <- ref@misc$atlas.palette
    } else {
      palette <- c("#edbe2a","#A58AFF","#53B400","#F8766D","#00B6EB","#d1cfcc","#FF0000","#87f6a5","#e812dd")
      if (nstates > length(palette)) {
        palette <- rainbow(n=nstates)
      }
    }  
  }
  #apply transparency to ref cells
  cols_use <- scales::alpha(palette, ref.alpha)
  
  if (is.null(query)) {
    p <- DimPlot(ref, reduction="umap", label = FALSE, group.by = labels.col, 
                 repel = TRUE, pt.size=ref.size, cols=cols_use, ...) +
      ggtitle ("Reference map") + theme(aspect.ratio=1)
  } else {
    p <- DimPlot(ref, reduction="umap", label = FALSE, group.by = labels.col,
                 repel = TRUE, pt.size=ref.size, cols=cols_use, ...) +
      geom_point(data.frame(query@reductions$umap@cell.embeddings), 
                 mapping=aes(x=UMAP_1,y=UMAP_2),alpha=0.6, size=pointsize,shape=17, color="gray10") +
      geom_density_2d(data=data.frame(query@reductions$umap@cell.embeddings), 
                      mapping=aes(x=UMAP_1,y=UMAP_2),color="black",n=200,h=2,size=linesize) +
      ggtitle ("Projection of query on reference map") + theme(aspect.ratio=1)
  }
  return(p)
}

# Function plot.statepred.composition
#' Summarize the predicted cell states of an object
#'
#' Makes a barplot of the frequency of cell states in a query object.
#'
#' @param ref Seurat object with the reference object
#' @param query Seurat object with query data
#' @param labels.col The metadata field used to annotate the clusters (default: functional.cluster)
#' @param metric One of `Count` or `Percent`. `Count` plots the absolute number of cells, `Percent` the fraction on the total number of cells.
#' @param cols Custom color palette for clusters
#' @return Barplot of predicted state composition
#' @examples
#' plot.statepred.composition(query_example.seurat)
#' @importFrom reshape2 melt
#' @import ggplot2
#' @export plot.statepred.composition
plot.statepred.composition = function(ref, query, labels.col="functional.cluster",cols=NULL, metric=c("Count","Percent")) {
  
  metric <- tolower(metric[1])
  
  labels <- ref[[labels.col]][,1]

  states_all <- levels(factor(labels))
  nstates <- length(states_all)
  
  if (!is.null(cols)) {  #custom palette
    if (nstates<=length(cols)) {
      palette <- cols[1:nstates]
    } else {  
      warning("Not enough colors provided. Making an automatic palette")
      palette <- rainbow(n=nstates)
    }  
  } else {   #default palette
    
    if (!is.null(ref@misc$atlas.palette)) {  #read directly from atlas, if stored
      palette <- ref@misc$atlas.palette
    } else {
      palette <- c("#edbe2a","#A58AFF","#53B400","#F8766D","#00B6EB","#d1cfcc","#FF0000","#87f6a5","#e812dd")
      if (nstates<=length(palette)) {
        palette <- palette[1:nstates]
      } else {   #make a new palette
        palette <- rainbow(n=nstates)
      }
    }
  }
  names(palette) <- states_all
  cols_use <- palette[states_all]

  tb <- table(factor(query[[labels.col]][,1], levels=states_all))
  
  if (metric=="percent") {  #normalize
    tb <- tb*100/sum(tb)
    tb.m <- reshape2::melt(tb)
    colnames(tb.m) <- c("Cell_state","Perc_cells")
    p <- ggplot(tb.m, aes(x=Cell_state, y=Perc_cells, fill=Cell_state)) + geom_bar(stat="identity") +
      theme_bw() + scale_fill_manual(values=cols_use) +
      theme(axis.text.x=element_blank(), legend.position="left")
  } else if (metric=="count") {
    tb.m <- reshape2::melt(tb)
    colnames(tb.m) <- c("Cell_state","Ncells")
    p <- ggplot(tb.m, aes(x=Cell_state, y=Ncells, fill=Cell_state)) + geom_bar(stat="identity") +
      theme_bw() + scale_fill_manual(values=cols_use) +
      theme(axis.text.x=element_blank(), legend.position="left")
  } else {
    stop("Unknown metric specified (Must be either 'count' or 'percent')")
  }
  
  return(p)
}

# Function plot.states.radar
#' Show expression level of key genes
#'
#' Makes a radar plot of the expression level of a set of genes. It can be useful to compare
#' the gene expression profile of different cell states in the reference atlas vs. a projected set.
#'
#' @param ref Reference Atlas
#' @param query Query data, either as a Seurat object or as a list of Seurat objects
#' @param labels.col The metadata field used to annotate the clusters
#' @param genes4radar Which genes to use for plotting
#' @param meta4radar Which metadata columns (numeric) to use for plotting. If not NULL, \code{genes4radar} are ignored
#' @param min.cells Only display cell states with a minimum number of cells
#' @param cols Custom color palette for samples in radar plot
#' @param ref.assay The assay to pull the reference expression data
#' @param query.assay The assay to pull the query expression data
#' @param return Return the combined plots instead of printing them to the default device (deprecated)
#' @param return.as.list Return plots in a list, instead of combining them in a single plot
#' @return Radar plot of gene expression of key genes by cell subtype
#' @usage plot.states.radar(ref)
#' @examples
#' ref <- load.reference.map()
#' plot.states.radar(ref)
#' @import ggplot2
#' @importFrom scales hue_pal
#' @importFrom patchwork wrap_plots plot_annotation
#' @export plot.states.radar
plot.states.radar = function(ref, query=NULL,
                             labels.col="functional.cluster",
                             ref.assay='RNA', query.assay='RNA',
                             genes4radar=c("Foxp3","Cd4","Cd8a","Tcf7","Ccr7","Gzmb",
                                              "Gzmk","Pdcd1","Havcr2","Tox","Mki67"),
                             meta4radar=NULL,
                             min.cells=50, cols=NULL,
                             return=FALSE, return.as.list=FALSE) {
  
  #Make sure query is a list
  if(!is.null(query) & !is.list(query)) {
    query <- list(Query=query)
  }
  
  #Check assays exist
  if (!ref.assay %in% Assays(ref)) {
    stop(sprintf("Assay %s not found in reference object. Please check ref.assay parameter", ref.assay))
  }
  
  #Whether to use gene expression or metadata
  if (!is.null(meta4radar)) {
    refmat <- t(ref[[]])
    feat.use <- intersect(meta4radar, row.names(refmat))
    if (length(feat.use)==0) {
      stop("None of the provided meta columns were found - check option 'meta4radar'")
    }
    feat.missing <- setdiff(meta4radar, feat.use)
    if (length(feat.missing)>0) {
      to.print <- paste(feat.missing, sep=",", collapse = ",")
      warning(sprintf("Some metadata columns were not found:\n%s", to.print))
    }
  } else {
    refmat <- ref@assays[[ref.assay]]@data

    #Check gene names/feature names
    feat.use <- intersect(genes4radar, row.names(refmat))
    #If overlap is zero, first check whether wrong species was used (upper case to human)
    if (length(feat.use)==0) {
      genes4radar <- toupper(genes4radar)
      feat.use <- intersect(genes4radar, row.names(refmat))
      if (length(feat.use)==0) {
        stop("None of the provided genes were found - check option 'genes4radar'")
      }
    }
    feat.missing <- setdiff(genes4radar, feat.use)
    if (length(feat.missing)>0) {
      to.print <- paste(feat.missing, sep=",", collapse = ",")
      warning(sprintf("Some gene symbols were not found:\n%s", to.print))
    }
  }
  
  order <- match(feat.use, row.names(refmat))
  
  rr <- as.matrix(refmat[order,])
  rr <- matrix(as.numeric(rr), ncol=ncol(rr))

  #Set colors
  ncolors <- 1+length(query)
  if (ncolors==1) {
    radar.colors <- "black"
  } else {
    if (is.null(cols)) {
      radar.colors <- c("black", hue_pal()(ncolors-1))
    } else {
      cols <- c("black", cols)
      if (ncolors <= length(cols)) {
        radar.colors <- cols[1:ncolors]
      } else {  
        warning("Not enough colors provided. Making an automatic palette")
        radar.colors <- c("black", hue_pal()(ncolors-1))
      }
    }
  }  
  names(radar.colors) <- c("Reference", names(query))
  

  labels <- ref[[labels.col]][,1]
  states_all <- levels(factor(labels))
  nstates <- length(states_all)
  
  if (!is.null(query)) {
    if (is.null(names(query))) {
      for (i in 1:length(query)) {
        names(query)[[i]] <- paste0("Query",i)
      }
    }
    labels.q <- list()
    qq <- list()
    
    for (i in 1:length(query)) {
      if (!labels.col %in% colnames(query[[i]]@meta.data)) {
         message1 <- sprintf("Could not find %s column in query object metadata.",labels.col)
         message2 <- "Did you run cellstate.predict() on this object to predict cell states?"
         stop(paste(message1, message2, sep="\n"))
      }
      labels.q[[i]] <- query[[i]][[labels.col]][,1]
      
      if (!is.null(meta4radar)) {
        qmat <- t(query[[i]][[]])
      } else {
        if (!query.assay %in% Assays(query[[i]])) {
          stop(sprintf("Assay %s not found in query object. Please check ref.assay parameter", query.assay))
        }
        qmat <- query[[i]]@assays[[query.assay]]@data
      }
      order <- match(feat.use, row.names(qmat))
      
      qq[[i]] <- as.matrix(qmat[order,])
      qq[[i]] <- matrix(as.numeric(qq[[i]]), ncol=ncol(qq[[i]]))
    }
  }
  
  #Get raw expression means, to normalize by gene
  m <- matrix(, nrow = length(states_all), ncol = length(feat.use))
  rownames(m) <- states_all
  colnames(m) <- feat.use
  for (i in 1:length(states_all)) {
    s <- states_all[i]
    m[i,] <- apply(rr[, labels == s], MARGIN=1, function(x){mean(x, na.rm=T)})
  }
  normfacs <- apply(m, MARGIN=2, function(x) {max(c(1,x), na.rm=T)})
  
  pll <- list()
  for (j in 1:length(states_all)) {
    s <- states_all[j]
    
    this.mean <- apply(rr[, labels == s], MARGIN=1, function(x){mean(x, na.rm=T)})
    this.mean <- this.mean/normfacs
    
    this.df <- data.frame(t(rbind(names(this.mean), this.mean, "Reference")))
    colnames(this.df) <- c("Gene","Expression","Dataset")
    this.df$Expression <- as.numeric(as.character(this.df$Expression))
    
    i <- 1
    while (i <= length(query)) {
      ll <- labels.q[[i]]
      m <- as.matrix(qq[[i]][, !is.na(ll) & ll == s])
      if (ncol(m) >= min.cells) {
        q.mean <- apply(m, MARGIN=1, function(x){mean(x, na.rm=T)})
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
    
    levs <- unique(this.df$Dataset)
    this.df$Dataset <- factor(this.df$Dataset, levels=levs)
    this.df$Gene <- factor(this.df$Gene, levels=feat.use)
    
    pll[[j]] <- ggplot(data=this.df,  aes(x=Gene, y=Expression, group= Dataset, colour=Dataset, fill=Dataset)) +
      geom_point(size=2) +
      geom_polygon(size = 0.75, alpha= 0.1) +
      ylim(ymin, ymax) + ggtitle(s)  +
      scale_x_discrete() +
      scale_fill_manual(values= radar.colors) +
      scale_colour_manual(values= radar.colors) +
      theme_light() +
      theme(axis.text.x=element_blank()) +
      annotate(geom="text", x=seq(1,length(feat.use)), y=ymax-0.05*ymax, label=feat.use, size=3) +
      coord_polar()
    
  }
  #Return plots
  if (return.as.list) {
    return(pll)
  } else {
    g <- wrap_plots(pll) + plot_annotation(paste0("Radar plots for ", labels.col))
    return(g)
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
#' @param print.n The number of top dimensions to return to STDOUT
#' @param verbose Print results to STDOUT
#' @return A dataframe, where rows are ICA/PCA dimensions. ICA/PCAs are ranked by statistical significance when comparing their distribution between query and control (or query vs. reference map)
#' @examples
#' find.discriminant.dimensions(ref, query=query.set)
#' find.discriminant.dimensions(ref, query=query.set, query.control=control.set)
#' @importFrom stats t.test ks.test
#' @export find.discriminant.dimensions
find.discriminant.dimensions <- function(ref, query, query.control=NULL, query.assay="RNA",
                                         state="largest", labels.col="functional.cluster",
                                         reduction="ICA", test=c("ks","t.test"), ndim=50, print.n=3, verbose=T) {

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
    if (!requireNamespace("fastICA", quietly = TRUE)) {
      stop("Please install package 'fastICA' to run this function.", call. = FALSE)
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
  
  df <- data.frame(matrix(ncol = 3, nrow = ndim))
  colnames(df) <- c("stat","stat_abs","p_val")
  rownames(df) <- colnames(perturb_dimRed)[1:ndim]
  
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
    
    ttest <- t.test(d1, d2)
    if (test=="ks") {
       this.test <- ks.test(d1, d2, alternative="two.sided")
       this.test.signed <- this.test$statistic * sign(ttest$statistic)  #KS test statistic has no sign
    } else {
       this.test <- ttest
       this.test.signed <- this.test$statistic
    }

    df[pc, "stat"] <- this.test.signed
    df[pc, "stat_abs"] <- abs(this.test.signed)
    df[pc, "p_val"] <- this.test$p.value * ndim   #multiple testing
  }
  df <- df[with(df, order(stat_abs, decreasing=T)), ]
  
  buffer <- ""
  for (i in 1:print.n) {
    topPC <- rownames(df)[i]
    pc.index <- match(topPC, colnames(perturb_dimRed))
    feats <- ref@reductions[[reduction]]@feature.loadings[,pc.index]
    topgenes <- names(head(sort(abs(feats), decreasing = T), 10))
    topgenes.sign <- feats[topgenes]
    if (df[i, "stat"]>0) {
       topgenes.p <- topgenes.sign[topgenes.sign>=0]
       topgenes.n <- topgenes.sign[topgenes.sign<0]
    } else {
      topgenes.p <- topgenes.sign[topgenes.sign<0]
      topgenes.n <- topgenes.sign[topgenes.sign>=0]
    }

    pval2print <- ifelse(df[i, "p_val"]<0.001, sprintf("%.1e",df[i, "p_val"]), sprintf("%.4f",df[i, "p_val"]))
    buffer <- paste0(buffer, sprintf("-----------------------\nTop %s component %i: %s Stat %.3f p.val %s\n",
                    reduction, i, topPC, df[i, "stat"], pval2print))
    
    buffer <- paste0(buffer, "Driver genes for this component:\n")
    buffer <- paste0(buffer, paste(c("Higher in query +++ ", names(topgenes.p)), collapse = " "), "\n")
    buffer <- paste0(buffer, paste(c("Lower in query  --- ", names(topgenes.n)), collapse = " "), "\n")
  }

  if (verbose) {
     message(buffer)
  }
  
  return(df)
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
#' @param labels.col The metadata field used to annotate the clusters
#' @param query.state Only plot the query cells from this specific state
#' @param extra.dim The additional dimension to be added on the z-axis of the plot. Can be either:
#' \itemize{
#'   \item{An ICA or PCA dimension (e.g. ICA_10). See `find.discriminant.dimensions`}
#'   \item{Any numeric metadata field associated to the cells (e.g. 'cycling.score')}
#' }
#' @return A three dimensional plot with UMAP_1 and UMAP_2 on the x and y axis respectively, and the specified `extra.dim` on the z-axis.
#' @examples
#' plot.discriminant.3d(ref, query=query, extra.dim="ICA_19")
#' plot.discriminant.3d(ref, query=treated.set, query.control=control.set, extra.dim="ICA_2")
#' @export plot.discriminant.3d

plot.discriminant.3d <- function(ref, query, query.control=NULL, query.assay="RNA",
                                 labels.col="functional.cluster", extra.dim="ICA_1", query.state=NULL) {
  
  require(plotly)
  reduction=NULL
  message(paste0("Generating UMAP with 3rd dimension on ",extra.dim))
  if (grepl("^ica_\\d+", tolower(extra.dim), perl=T)) {
    reduction="ica"
  } else if (grepl("^pca_\\d+", tolower(extra.dim), perl=T)) {
    reduction="pca"
  }
  
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
  
  
  metacols <- intersect(names(ref.3d@meta.data),names(query@meta.data))
  ref.3d@meta.data <- rbind(ref.3d@meta.data[,metacols], query@meta.data[,metacols])
  
  ref.3d@reductions$umap@cell.embeddings <- rbind(ref.3d@reductions$umap@cell.embeddings, 
                                                  query@reductions$umap@cell.embeddings)
  
  
  ref.3d@meta.data <- cbind(ref.3d@meta.data, ref.3d@reductions$umap@cell.embeddings)
  
  
  #Calculate ICA or PCA embeddings for query
  if (!is.null(reduction)) {
    if (reduction=="ica"){
      projected <- apply.ica.obj(query=query, query.assay=query.assay, ica.obj=ref.3d@misc$ica)
      ref.3d@reductions[[reduction]]@cell.embeddings <- rbind(ref.3d@reductions[[reduction]]@cell.embeddings, projected)
    } else if (reduction=="pca") {
      projected <- apply.pca.obj.2(query=query, query.assay=query.assay, pca.obj=ref.3d@misc$pca_object)
      ref.3d@reductions[[reduction]]@cell.embeddings <- rbind(ref.3d@reductions[[reduction]]@cell.embeddings, projected)
    }
  }
  
  # Add metadata column 
  if (is.null(reduction)) {
    if (extra.dim %in% colnames(ref.3d@meta.data)){
      ref.3d@meta.data <- cbind(ref.3d@meta.data, Discriminant = ref.3d@meta.data[,extra.dim])
    } else {
      stop(sprintf("extra.dim %s not present in meta.data", extra.dim))
    }
  } else {
    ref.3d@meta.data <- cbind(ref.3d@meta.data,Discriminant=ref.3d@reductions[[reduction]]@cell.embeddings[,extra.dim])
  }
  
  if (is.null(query.control)) {
    cols <- c(Reference="gray50",Query="red")
  } else {
    cols <- c(Reference="gray50",Control="green",Query="red")
  }
  
  plotting.data <- ref.3d@meta.data[,c("UMAP_1", "UMAP_2","Discriminant", labels.col,"queryGroup")]
  
  plotting.data$size <- ifelse(plotting.data$queryGroup == "Reference",0.3,6)
  
  g <- plotly::plot_ly(data = plotting.data,
               x = ~UMAP_1, y = ~UMAP_2, z = ~Discriminant,
               color = ~queryGroup,
               type = "scatter3d",
               mode = "markers",
               text=~queryGroup,
               hoverinfo="text",
               alpha=0.6,
               alpha_stroke=0.6,
               size=~size,
               colors=cols ) |> plotly::layout(
                                  title = sprintf("Projection of query on reference map + dimension %s", extra.dim),
                                  scene = list(
                                     xaxis = list(title = "UMAP_1"),
                                     yaxis = list(title = "UMAP_2"),
                                     zaxis = list(title = extra.dim)
                                  ))
  print(g)
  return(g)
}

# Function find.discriminant.genes
#' Find discriminant genes
#'
#' Based on `FindMarkers`. It performs differential expression analysis between a projected query and a control (either the reference map or a control sample), for
#' a given cell type. Useful to detect whether specific cell states over/under-express genes between conditions or with respect to the reference.
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
#' @param test Type of test for DE analysis. See help for `FindMarkers` for implemented tests.
#' @param min.cells Minimum number of cells in the cell type to proceed with analysis.
#' @param genes.use What subset of genes to consider for DE analysis:
#' \itemize{
#'   \item{"variable" - Only consider variable genes of the reference}
#'   \item{"all" - Use intersection of all genes in query and control}
#'   \item{A custom list of genes}
#' }
#' @param ... Adding parameters for `FindMarkers`
#' @return A dataframe with a ranked list of genes as rows, and statistics as columns (e.g. log fold-change, p-values). See help for `FindMarkers` for more details.
#' @examples
#' # Discriminant genes between query and reference in cell type "Tex"
#' markers <- find.discriminant.genes(ref, query=query.set, state="Tex")
#' 
#' # Discriminant genes between query and control sample in most represented cell type
#' markers <- find.discriminant.genes(ref, query=query.set, query.control=control.set)
#' 
#' # Pass results to EnhancedVolcano for visual results
#' library(EnhancedVolcano)
#' EnhancedVolcano(markers, lab = rownames(markers), x = 'avg_logFC', y = 'p_val')
#' 
#' @import Seurat
#' @export find.discriminant.genes
#' 
find.discriminant.genes <- function(ref, query, query.control=NULL, query.assay="RNA",
                                    state="largest", labels.col="functional.cluster",
                                    test="wilcox", min.cells=10, genes.use=c("variable","all"), ...)
{  
  
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
    message(paste0("Performing differential expression analysis using the most abundant cell state in the query - ", state))
  } else if (state=="all") {
    message("Performing differential expression analysis using all cells")
  } else {
    if (!state %in% query[[labels.col]][,1]) {
      stop(sprintf("State %s not found in query metadata colum %s", state, labels.col))
    }
    if (!is.null(query.control) & !state %in% query.control[[labels.col]][,1]) {
      stop(sprintf("State %s not found in query.control metadata colum %s", state, labels.col))
    }
    message(paste0("Performing differential expression analysis with user-specified state - ", state))
  }
  #Subset query data on specific state
  if (state=="all") {
    s1.cells <- colnames(query)
    if (!is.null(query.control)) {
      s2.cells  <- colnames(query.control)
    } else {
      s2.cells <- colnames(ref)
    }
  } else {
    s1.cells <- colnames(query)[which(query[[labels.col]]==state)]
    
    if (!is.null(query.control)) {
      s2.cells <- colnames(query.control)[which(query.control[[labels.col]]==state)]
    } else {
      s2.cells <- colnames(ref)[which(ref[[labels.col]]==state)]
    }
  }
  
  #check we have enough cells
  if (length(s1.cells)<min.cells) {
    stop(sprintf("Too few cells for state %s in query. Exit.", state))
  }
  if (!is.null(query.control) & length(s2.cells)<min.cells ) {
    stop(sprintf("Too few cells for state %s in query control. Exit.", state))
  }
  
  #Subset on subtype
  DefaultAssay(query) <- query.assay
  s1 <- subset(query, cells=s1.cells)
  s1$Group <- "Query"
  
  if (!is.null(query.control)) {
    DefaultAssay(query.control) <- query.assay
    s2 <- subset(query.control, cells=s2.cells)
  } else {
    DefaultAssay(ref) <- query.assay
    s2 <- subset(ref, cells=s2.cells)
  }
  s2$Group <- "Control"
  
  s.m <- merge(s1, s2)
  Idents(s.m) <- "Group"
  
  #use all genes or only variable genes from the reference
  if (genes.use[1] == "all") {
    which.genes <- NULL
  } else if (genes.use[1] == "variable") {
    which.genes <- intersect(ref@assays$integrated@var.features, rownames(s.m))
  } else {
    which.genes <- intersect(genes.use, rownames(s.m))
  }
  
  markers <- FindMarkers(s.m, slot="data", ident.1="Query", ident.2="Control", only.pos = F, test.use=test, assay=query.assay,
                         features = which.genes, ...)
  
  
  return(markers)
}


#' Make a ProjecTILs reference
#'
#' Converts a Seurat object to a ProjecTILs reference atlas. You can preserve your low-dimensionality embeddings
#' (e.g. UMAP) in the reference atlas by setting `recalculate.umap=FALSE`, or recalculate the UMAP using one of
#' the two methods (\link[umap]{umap::umap} or  \link[uwot]{uwot::umap}). Recalculation allows exploting the
#' 'predict' functionalities of these methods for embedding of new points; skipping recalculation will 
#' make the projection use an approximation for UMAP embedding of the query.
#'
#' @param ref Seurat object with reference atlas
#' @param assay The data slot where to pull the expression data
#' @param atlas.name An optional name for your reference
#' @param annotation.column The metadata column with the cluster annotations for this atlas
#' @param recalculate.umap If TRUE, run the `umap` or `uwot` algorithm to generate embeddings.
#'   Otherwise use the embeddings stored in the `dimred` slot.
#' @param umap.method Which method to use for calculating the umap reduction
#' @param metric Distance metric to use to find nearest neighbors for UMAP
#' @param min_dist Effective minimum distance between UMAP embedded points
#' @param n_neighbors Size of local neighborhood for UMAP
#' @param ndim Number of PCA dimensions
#' @param dimred Use the pre-calculated embeddings stored at `Embeddings(ref, dimred)`
#' @param nfeatures Number of variable features (only calculated if not already present)
#' @param color.palette A (named) vector of colors for the reference plotting functions.
#'     One color for each cell type in 'functional.cluster'
#' @param scGate.model.human A human \link[scGate]{scGate} model to purify the cell types represented in the
#'     map. For example, if the map contains CD4 T cell subtype, specify an scGate model for CD4 T cells.
#' @param scGate.model.human A mouse \link[scGate]{scGate} model to purify the cell types represented in the
#'     map.
#' @param store.markers Whether to store the top differentially expressed genes in `ref@@misc$gene.panel`     
#' @param n.markers Store the top `n.markers` for each subtype given by differential
#'     expression analysis
#' @param seed Random seed
#' @return A reference atlas compatible with ProjecTILs
#' @examples 
#' custom_reference <- ProjecTILs::make.reference(myref, recalculate.umap=T)  
#' @importFrom stats prcomp
#' @importFrom uwot umap
#' @importFrom dplyr group_by top_n
#' @export make.reference
#' 
make.reference <- function(ref,
                           assay=NULL,
                           atlas.name="custom_reference",
                           annotation.column="functional.cluster",
                           recalculate.umap=FALSE,
                           umap.method=c("umap","uwot"),
                           metric="cosine",
                           min_dist=0.3,
                           n_neighbors = 30,
                           ndim=20,
                           dimred="umap",
                           nfeatures=1000,
                           color.palette=NULL,
                           scGate.model.human=NULL,
                           scGate.model.mouse=NULL,
                           store.markers=FALSE,
                           n.markers=10,
                           seed=123) {
  if (is.null(assay)) {
    assay=DefaultAssay(ref)
  }
  if (is.null(ref@assays[[assay]])) {
    stop(sprintf("Assay %s not found in reference object. Select a different assay", assay))
  }
  DefaultAssay(ref) <- assay
  
  if ("var.features" %in% slotNames(ref@assays[[assay]]) & !is.null(ref@assays[[assay]]@var.features)) {
    varfeat <- ref@assays[[assay]]@var.features
  } else {
    ref <- FindVariableFeatures(ref, assay = assay, nfeatures = nfeatures, verbose=FALSE)
    varfeat <- ref@assays[[assay]]@var.features
  } 
  
  #Recompute PCA embeddings using prcomp
  ref <- prcomp.seurat(ref, ndim=ndim, assay=assay)

  #Recalculate UMAP, or use an already-present dimensionality reduction
  if (!recalculate.umap) {
    if (dimred %in% names(ref@reductions)) {
      ref.pca <- ref@misc$pca_object
      cell.order = rownames(ref.pca$x)
      low.emb <- ref@reductions[[dimred]]@cell.embeddings[cell.order,]
      colnames(low.emb) <- c("UMAP_1","UMAP_2")
      #Save these embeddings
      ref@misc$umap_object <- list()
      ref@misc$umap_object$data <- ref.pca$x
      ref@misc$umap_object$layout <- low.emb
      
    } else {
      stop(sprintf("Dimred %s not found in reference object. Select a different dimensionality reduction, or set recalculate.umap=TRUE to compute UMAP coordinates", dimred))
    }
  } else {
    
    umap.method = umap.method[1]
    ref.pca <- ref@misc$pca_object
    
    if (umap.method == "umap") {
      #generate UMAP embeddings
      ref.umap <- run.umap.2(ref.pca, ndim=ndim, seed=seed,
                             n.neighbors=n_neighbors, min.dist=min_dist,
                             metric=metric)
      
      #Save UMAP object
      ref@misc$umap_object <- ref.umap
      ref@reductions$umap@cell.embeddings <- ref.umap$layout
    } else if (umap.method == "uwot") {
      warning("There are known issues with saving a loading uwot models. If you plan to save your reference as an .rds file, please use umap.method='umap'")
      #generate UMAP embeddings
      ref.umap <- run.umap.uwot(ref.pca, ndim=ndim, seed=seed,
                             n.neighbors=n_neighbors, min.dist=min_dist,
                             metric=metric)
      
      ref.umap$data <- ref.pca$x
      #Save UMAP object
      ref@misc$umap_object <- ref.umap
      ref@reductions$umap@cell.embeddings <- ref.umap$embedding
    } else {
      stop("Unsupported UMAP method.")
    }
  }
  
  if (!annotation.column == "functional.cluster") {
    ref$functional.cluster <- ref@meta.data[,annotation.column]
  }
  ref$functional.cluster <- factor(ref$functional.cluster)
  levs <- levels(ref$functional.cluster)
  
  #Reference name
  ref@misc$projecTILs=atlas.name
  
  #Add color palette
  if (!is.null(color.palette)) {
    if (length(color.palette) != length(levs)) {
      stop("Length of color palette must match number of cell types")
    }
    if (!is.null(names(color.palette))) {
      d <- setdiff(names(color.palette), levels(ref$functional.cluster))
      if (length(d)>0) {
        stop("Names of color palette do not match annotation levels")
      }
    } else {
      names(color.palette) <- levs
    }
  } else {
    color.palette <- rainbow(n=length(levs))
    names(color.palette) <- levs
  }
  ref@misc$atlas.palette <- color.palette
  
  #Add scGate models
  ref@misc$scGate <- list()
  if (!is.null(scGate.model.human)) {
    ref@misc$scGate$human <- scGate.model.human
  }
  if (!is.null(scGate.model.mouse)) {
    ref@misc$scGate$mouse <- scGate.model.mouse
  }
  
  #Add DEGs
  Idents(ref) <- "functional.cluster"
  
  if (store.markers) {
    DefaultAssay(ref) <- "RNA"
    
    cluster.markers <- FindAllMarkers(ref, only.pos = T, min.pct = 0.1, min.diff.pct=0.1, 
                                      logfc.threshold = 0.25, max.cells.per.ident = 500,
                                      test.use="wilcox",base=exp(1), verbose=F)
    
    all <- cluster.markers |> dplyr::group_by(cluster) |>
      dplyr::top_n(n = n.markers, wt = abs(avg_logFC))
    
    markers <- list()
    for (i in levels(ref@active.ident)) {
      markers[[i]] <- subset(all, cluster==i)$gene
    }  
    ref@misc$gene.panel <- markers
  }
  
  #Store in integrated assay, to be understood by ProjecTILs
  names(ref@assays)[names(ref@assays)==assay] = "integrated"
  DefaultAssay(ref) <- "integrated"
  
  return(ref)
}

#' Merge Seurat objects, including reductions (e.g. PCA, UMAP, ICA)
#'
#' Given two Seurat objects, merge counts and data as well as dim reductions (PCA, UMAP, ICA, etc.)
#'
#' @param x First object to merge
#' @param y Second object to merge
#' @param merge.dr Which reductions to merge. By default merges all reductions shared by the two objects.
#' @param ... More parameters to \link{merge} function
#' @return A merged Seurat object
#' @examples 
#' seurat.merged <- merge.Seurat.embeddings(obj.1, obj.2)  
#' #To merge multiple object stored in a list
#' seurat.merged <- Reduce(f=merge.Seurat.embeddings, x=obj.list)
#' @import Seurat
#' @export merge.Seurat.embeddings

merge.Seurat.embeddings <- function(x=NULL, y=NULL, merge.dr=NULL, ...)
{ 
  if (is.null(merge.dr)) {  #merge all shared reductions, by default
    merge.dr <- intersect(names(x@reductions), names(y@reductions))
  }
  #Use Seurat merge function, inheriting parameters
  m <- merge(x, y, merge.dr=merge.dr, ...)

  return(m)
}

#' Recalculate low dimensional embeddings after projection
#'
#' Given a reference object and a (list of) projected objects, recalculate low-dim embeddings accounting for the projected cells
#'
#' @param ref Reference map
#' @param projected A projected object (or list of projected objects) generated using \link{make.projection}
#' @param ref.assay Assay for reference object
#' @param proj.assay Assay for projected object(s)
#' @param ndim Number of dimensions for recalculating dimensionality reductions
#' @param n.neighbors Number of neighbors for UMAP algorithm
#' @param min.dist Tightness parameter for UMAP embedding
#' @param metric Distance metric to use to find nearest neighbors for UMAP
#' @param recalc.pca Whether to recalculate the PCA embeddings with the combined reference and projected data
#' @param umap.method Which method should be used to calculate UMAP embeddings
#' @param resol Resolution for unsupervised clustering
#' @param k.param Number of nearest neighbors for clustering
#' @param seed Random seed for reproducibility
#' @return A combined reference object of reference and projected object(s), with new low dimensional embeddings
#' @examples 
#' combined <- recalculate.embeddings(ref, projected, ndim=10)
#' @export recalculate.embeddings

recalculate.embeddings <- function(ref, projected, ref.assay="integrated", proj.assay="integrated",
                                   ndim=NULL, n.neighbors=20, min.dist=0.3, recalc.pca=FALSE,
                                   resol=0.4, k.param=15, metric="cosine",
                                   umap.method=c('umap','uwot'), seed=123){ 
  
  if (is.null(ref) | is.null(projected)) {
    stop("Please provide a reference and a projected object (or list of projected objects)")
  }
  
  if (is.list(projected) | is.vector(projected)) {
    projected <- Reduce(f=merge.Seurat.embeddings, x=projected)
  } 
  projected$ref_or_query <- "query"
  ref$ref_or_query <- "ref"
  umap.method <- umap.method[1]
  
  projected <- RenameCells(object = projected, add.cell.id = "Q")

  merged <- merge.Seurat.embeddings(ref, projected)
  
  DefaultAssay(merged) <- proj.assay
  DefaultAssay(ref) <- ref.assay
  
  if (is.null(ndim)) {
    ndim <- ncol(ref@misc$umap_object$data)
    if (is.null(ndim)) {
      stop("Please provide number of dimensions for dimensionality reduction (ndim)")
    }
  }
  
  VariableFeatures(merged) <- VariableFeatures(ref)
  merged@misc <- ref@misc
  merged@misc$pca_object$x <- merged@reductions$pca@cell.embeddings
  
  if (recalc.pca) {
    message("Recalculating PCA embeddings...")
    merged <- prcomp.seurat(merged, assay=proj.assay, ndim=ndim)
  }
  
  ref.pca <- merged@misc$pca_object
  
  if (umap.method=="uwot") {  
    message("Recalculating UMAP embeddings using uwot...")
    warning("There are known issues with saving a loading uwot models. If you plan to save your reference as an .rds file, please use umap.method='umap'")
    
    ref.umap <- run.umap.uwot(ref.pca, ndim=ndim,
                              n.neighbors = n.neighbors,
                              min.dist=min.dist,
                              seed=seed, metric=metric)
    
    ref.umap$data <- ref.pca$x
    #Save UMAP object
    merged@misc$umap_object <- ref.umap
    merged@reductions$umap@cell.embeddings <- ref.umap$embedding
    
  } else if (umap.method=="umap") {
    message("Recalculating UMAP embeddings using 'umap' package...")
    
    ref.umap <- run.umap.2(ref.pca, ndim=ndim,
                           n.neighbors = n.neighbors,
                           min.dist=min.dist,
                           seed=seed, metric=metric)
    #Save UMAP object
    merged@misc$umap_object <- ref.umap
    merged@reductions$umap@cell.embeddings <- ref.umap$layout
    
  } else {
    stop("UMAP method not supported.")
  }
  
  #Did any new clusters arise after adding projected data?
  DefaultAssay(merged) <- ref.assay
  merged <- FindNeighbors(merged, reduction = "pca", dims = 1:ndim,
                          k.param = k.param, verbose=FALSE)
  merged <- FindClusters(merged, resolution = resol, verbose=FALSE)
  
  tab <- table(merged$seurat_clusters, merged$ref_or_query)
  
  glob.freq <- table(merged$ref_or_query)["query"]/ncol(merged)
  freq <- apply(tab, 1, function(x){x/sum(x)})["query",] - glob.freq
  freq[freq<0] <- 0
  merged$newclusters <- freq[merged$seurat_clusters]
  
  Idents(merged) <- "functional.cluster"
  return(merged)
}

#' Calculate Silhouette coefficient
#'
#' Given a projected object and its reference, calculate silhouette coefficient for query cells with respect
#' to reference cells with the same cell labels.
#'
#' @param ref Reference object
#' @param query Query object. If not specified, the silhouette coefficient of only the reference will be calculated
#' @param reduction Which dimensionality reduction to use for euclidian distance calculation
#' @param ndim Number of dimensions in the dimred to use for distance calculation. If NULL, use all dimensions.
#' @param label_col Metadata column with cell type annotations. Must be present both in reference and query
#' @param normalize.scores Whether to normalize silhouette scores by the average cell type silhouettes of the reference
#' @param min.cells Only report silhouette scores for cell type with at least this number of cells
#' @return A dataframe with average silhouette coefficient for each cell type
#' @examples
#' data(query_example_seurat)
#' ref <- load.reference.map()
#' q <- Run.ProjecTILs(query_example_seurat, ref=ref, fast.umap.predict=TRUE)
#' combined <- compute_silhouette(ref, query=q)
#' @importFrom pracma distmat
#' @export compute_silhouette

compute_silhouette <- function(ref, query=NULL,
                               reduction="pca",
                               ndim=NULL,
                               label_col="functional.cluster",
                               normalize.scores=FALSE,
                               min.cells=20) {
  
  y <- Reductions(ref, slot=reduction)
  if (is.null(ndim)) {
    ndim <- ncol(y@cell.embeddings)
  }
  y <- y@cell.embeddings[,1:ndim]
  
  levs <- levels(as.factor(ref@meta.data[,label_col]))
  labs.y <- ref@meta.data[,label_col]
  
  if (is.null(query)) {
    x <- y
    labs.x <- labs.y
  } else {
    
    subtypes <- table(query@meta.data[,label_col])
    subtypes <- names(subtypes[subtypes > min.cells])
    cid <- Cells(query)[query@meta.data[,label_col] %in% subtypes]
    
    query <- subset(query, cells = cid)
    
    x <- Reductions(query, slot=reduction)
    x <- x@cell.embeddings[,1:ndim]
    
    labs.x <- factor(query@meta.data[,label_col], levels=levs)
  }
  
  dists <- pracma::distmat(x,y)
  
  sil <- silhouette_2sets(dists, labs.x, labs.y)
  
  #summarize by cluster
  means <- aggregate(sil$sil_width, list(sil$cluster), FUN=mean) 
  colnames(means) <- c("Cluster","Silhouette")
  
  if (normalize.scores) { #normalize by silhouette of the reference
    dist.ref <- pracma::distmat(y,y)
    sil.ref <- silhouette_2sets(dist.ref, labs.y, labs.y)
    means.ref <- aggregate(sil.ref$sil_width, list(sil.ref$cluster), FUN=mean) 
    colnames(means.ref) <- c("Cluster","Silhouette")
    
    for (i in 1:nrow(means)) {
      subset <- means[i,"Cluster"]
      j <- which(means.ref$Cluster == subset)
      means[i,"Silhouette.norm"] <- means[i,"Silhouette"]/means.ref[j,"Silhouette"]
      if (means[i,"Silhouette.norm"] > 1) {means[i,"Silhouette.norm"]=1}
      if (means[i,"Silhouette.norm"] < 0) {means[i,"Silhouette.norm"]=0}
    }
  }
  return(means)
}

#' Project a query scRNA-seq dataset onto a reference atlas
#'
#' This function allows projecting ("query") single-cell RNA-seq datasets onto a reference map
#' (i.e. a curated and annotated scRNA-seq dataset). 
#' To project multiple datasets, submit a list of Seurat objects with the query parameter.
#' The projection consists of 3 steps: 
#' \itemize{
#'  \item{pre-processing: optional steps which might include pre-filtering of cells by markers using `scGate`,
#' data normalization, and ortholog conversion.}
#'  \item{batch-effect correction: uses built-in STACAS algorithm to detect and correct for batch effects
#' (this step assumes that at least a fraction of the cells in the query are in the same state than cells in
#' the reference)}
#'  \item{embedding of corrected query data in the reduced-dimensionality spaces (PCA and UMAP) of the reference map.}
#' }
#' This function acts as a wrapper for \link{make.projection} and \link{cellstate.predict}
#' 
#' See \link{load.reference.map} to load or download a reference atlas. See
#' also \link{ProjecTILs.classifier} to use ProjecTILs as a cell type classifier.
#'
#' @param query Query data, either as single Seurat object or as a list of Seurat object
#' @param ref Reference Atlas - if NULL, downloads the default TIL reference atlas
#' @param filter.cells Pre-filter cells using `scGate`. Only set to FALSE if the dataset has 
#'     been previously subset to cell types represented in the reference.
#' @param split.by Grouping variable to split the query object (e.g. if the object contains multiple samples)    
#' @param reduction The dimensionality reduction used to assign cell type labels
#' @param ndim The number of dimensions used for cell type classification
#' @param k Number of neighbors for cell type classification
#' @param labels.col The metadata field of the reference to annotate the clusters
#' @param ... Additional parameters to \link[ProjecTILs]{make.projection}
#' @return An augmented Seurat object with projected UMAP coordinates on the reference map and cell classifications
#' @examples
#' data(query_example_seurat)
#' ref <- load.reference.map()
#' q <- Run.ProjecTILs(query_example_seurat, ref=ref, fast.umap.predict=TRUE)
#' plot.projection(ref=ref, query=q)
#' @export Run.ProjecTILs
Run.ProjecTILs <- function(query, ref=NULL,
                           filter.cells = TRUE,
                           split.by = NULL,
                           reduction="pca",
                           ndim=NULL, k=20,
                           labels.col="functional.cluster", ...) {
    
    if (!is.null(split.by)) {
        if (!split.by %in% colnames(query[[]])) {
            stop(sprintf("No grouping variable %s available in metadata", split.by))
        }
        query <- SplitObject(query, split.by = split.by)
    }
    
    #Run projection
    query <- make.projection(query=query, ref=ref, filter.cells=filter.cells, ...)
    
    if(!is.list(query)) {
        query <- list(query=query)
    }
    #Cell type classification
    query <- lapply(query, function(x) {
        cellstate.predict(ref=ref, query=x,
                          reduction=reduction,
                          ndim=ndim, k=k,
                          labels.col = labels.col)
    })
    #Merge embeddings
    if (length(query)==1 || !is.null(split.by)) {
       query <- Reduce(merge.Seurat.embeddings, query)
    }
    query
}

#' Annotate query dataset using a reference object
#'
#' Apply label transfer to annotate a query dataset with the cell types of a reference object.
#' Compared to \link{Run.ProjecTILs}, only cell labels are returned. The low-dim embeddings of
#' the query object (PCA, UMAP) are not modified.
#' 
#' See \link{load.reference.map} to load or download a reference atlas. 
#' See \link{Run.ProjecTILs} to embed the query in the same space of the reference
#'
#' @param query Query data stored in a Seurat object
#' @param ref Reference Atlas - if NULL, downloads the default TIL reference atlas
#' @param filter.cells Pre-filter cells using `scGate`. Only set to FALSE if the dataset has 
#'     been previously subset to cell types represented in the reference.
#' @param split.by Grouping variable to split the query object (e.g. if the object contains multiple samples)  
#' @param reduction The dimensionality reduction used to assign cell type labels
#' @param ndim The number of dimensions used for cell type classification
#' @param k Number of neighbors for cell type classification
#' @param labels.col The metadata field with label annotations of the reference, which will
#' be transfered to the query dataset
#' @param overwrite Replace any existing labels in \code{labels.col} with new labels.
#'     This may be useful for predicting cell types using multiple reference maps; run
#'     this function with \code{overwrite=FALSE} to combine existing labels
#'     with new labels from a second reference map.
#' @param ... Additional parameters to \link[ProjecTILs]{make.projection}
#' @return The query object with an additional metadata column containing predicted cell labels.
#' If cells are filtered prior to projection, they will be labeled as 'NA'
#' @examples
#' data(query_example_seurat)
#' ref <- load.reference.map()
#' q <- ProjecTILs.classifier(query_example_seurat, ref=ref)
#' table(q$functional.cluster, useNA="ifany")
#' @export ProjecTILs.classifier
ProjecTILs.classifier <- function(query, ref=NULL,
                           filter.cells = TRUE,
                           split.by = NULL,
                           reduction="pca",
                           ndim=NULL, k=20,
                           labels.col="functional.cluster",
                           overwrite=TRUE,
                           ...) {
  
  fast.umap.predict <- TRUE
  #only needed if we want to predict labels based on UMAP neighbors
  if (reduction=="umap") { fast.umap.predict <- FALSE }
  
  if(is.list(query)) {
     stop("Query must be a single Seurat object")
  }
  
  current.labs <- NULL
  if (labels.col %in% colnames(query[[]])) {
    current.labs <- query[[labels.col]]
  }
  
  if (!is.null(split.by)) {
      if (!split.by %in% colnames(query[[]])) {
          stop(sprintf("No grouping variable %s available in metadata", split.by))
      }
      q <- SplitObject(query, split.by = split.by)
  } else {
      q <- query
  }
  #Run projection
  q <- make.projection(query=q, ref=ref, filter.cells=filter.cells,
                       fast.umap.predict = fast.umap.predict, ...)
  
  if(!is.list(q)) {
      q <- list(query=q)
  }
  
  #Cell type classification
  q <- lapply(q, function(x) {
      cellstate.predict(ref=ref, query=x,
                        reduction=reduction,
                        ndim=ndim, k=k,
                        labels.col = labels.col)
  })
  
  #Merge embeddings
  q <- Reduce(merge.Seurat.embeddings, q)
  
  #Transfer labels to original query object
  labs <- q[[labels.col]]
  
  if (overwrite) {
    new.labs <- labs[[labels.col]]
    names(new.labs) <- rownames(labs)
  } else {
    new.labs <- combine_labels(current.labs, labs)
  }
  
  query@meta.data[,labels.col] <- NA
  query@meta.data[names(new.labs),labels.col] <- new.labs
  
  query
}


#' Plot a averaged expression heatmap from a Seurat object
#'
#' This function allows to plot a averaged expression heatmap starting from a Seurat object,
#' split by possibly any metadata present in the starting Seurat object
#'
#' The function first loads the pheatmap package, sets a seed for reproducibility, and selects the first element of the "method" parameter as the clustering method to be used
#' It then selects the wanted metadata from the data set, removing any missing values if the "remove.NA.meta" parameter is set to TRUE. <- 
#' Next, the function filters the data set to only include samples that have at least "min.cells" in the "metaSubset" variable. It then calculates the mean expression of the selected genes by the "metaSubset" variable.
#' The function then reorders the data frame if the "order.by" parameter is not null, and sets up the colors for the heatmap.
#' Finally, the function creates the heatmap using the pheatmap package, with the "method", "scale", "flip", "cluster_genes", "cluster_samples", and "show_samplenames" parameters controlling various aspects of the heatmap.
#'
#' @param data A Seurat object to be used for the heatmap
#' @param assay A string indicating the assay type, default is "RNA"
#' @param genes A vector of genes to be used in the heatmap
#' @param ref A ProjecTILs reference Seurat object to define the order of functional.cluster
#' @param scale A string indicating the scale of the heatmap, default is "row"
#' @param method A string or vector of strings indicating the clustering method to be used, default is "ward.D2"
#' @param brewer.palette A string indicating the color palette to be used, default is "RdBu"
#' @param palette_reverse A boolean indicating if color palette should be reversed, default is FALSE
#' @param cluster.col A string indicating the column name of the functional cluster to be used
#' @param metadata A vector of metadata to be used in the heatmap
#' @param order.by A string indicating the column name to reorder the heatmap
#' @param flip A boolean indicating if the heatmap should be flipped, default is FALSE
#' @param show_samplenames A boolean indicating whether the heatmap should display the sample names or not, default is FALSE
#' @param cluster_genes A boolean indicating if genes should be clustered, default is FALSE
#' @param cluster_samples A boolean indicating if samples should be clustered, default is FALSE
#' @param min.cells A value defining the minimum number of cells a sample should have to be kept, default is 10
#' @param remove.NA.meta A boolean indicating if missing samples with missing metadata should be plotted, default is TRUE
#' @param palette A named list containing colors vectors compatible with pheatmap. The list is named by the metadata names, default is taking these palettes to plot metadata: "Paired","Set2","Accent","Dark2","Set1","Set3".
#' @return A pheatmap plot, displaying averaged expression values across genes for each selected genes and samples.
#' @import pheatmap
#' @importFrom tidyr drop_na
#' @import RColorBrewer
#' @examples
#' library(Seurat)
#' ref <- load.reference.map(ref = "https://figshare.com/ndownloader/files/38921366")
#' celltype.heatmap(ref, assay = "RNA", genes = c("LEF1","SELL","GZMK","FGFBP2","HAVCR2","PDCD1","XCL1","KLRB1"), ref = ref, cluster.col = "functional.cluster", metadata = c("orig.ident", "Tissue"), order.by = "Tissue")
#' @export celltype.heatmap
celltype.heatmap <- function(data, assay="RNA", genes, ref = NULL, scale="row", 
                         method=c("ward.D2","ward.D", "average"), brewer.palette="RdBu", palette_reverse=F, palette = NULL,
                         cluster.col = "functional.cluster", metadata = NULL, order.by = NULL, flip=FALSE,
                         cluster_genes = FALSE, cluster_samples=FALSE, min.cells = 10,
                         show_samplenames = FALSE, remove.NA.meta = TRUE, breaks = seq(-2, 2, by = 0.1), ...
                         ) {
  
  set.seed(123)
  
  # Select clustering method to be used
  method = method[1]

  
  # Select desired metadata
  if (is.null(metadata)) {
    stop("Must at least provide one metadata")
  } else {
    meta.sub <-
      data@meta.data[, which(colnames(data@meta.data) %in% cluster.col), drop = F]
    meta.sub <-
      cbind(meta.sub, data@meta.data[, metadata, drop = F])
  }
  
  # Transform "NA" into true NAs
  meta.sub[meta.sub=="NA"] = NA
  # Remove NAs from metadata
  if(remove.NA.meta == TRUE){
    meta.sub <- meta.sub |> drop_na()
  }
  
  # Filters the data set to only include samples that have at least "min.cells" in the "metaSubset" variable.
  data$metaSubset <- factor(apply(meta.sub,1,paste,collapse="!"))
  
  t <- table(data$metaSubset)
  accept <- names(t)[t>min.cells]
  
  data <- subset(data, subset=metaSubset %in% accept)
  
  # Calculate mean expression by cluster
  m <- c()
  genes.removed <- c()
  for( g in unique(genes)){
    if(g %in% rownames(data@assays[[assay]])){
      m[[g]] <- tapply(data@assays[[assay]][g,],data$metaSubset, mean)
    }
    else{
      genes.removed <- c(genes.removed, g)
    }
  }
  if(length(genes.removed)>0){
    cat("These genes were not found in the assay and were excluded from plotting:", genes.removed)
  }
  
  m <- as.data.frame(m)
  
  m <- m[accept,]
  
  # Compute metadata for the annotation colors
  m.subset <- factor(unlist(lapply(strsplit(rownames(m),"!",perl = T),function(x) x[[1]])))
  m.meta <- list()
  for (i in 1:length(metadata)){
    m.meta[[i]] <- factor(unlist(lapply(strsplit(rownames(m),"!",perl = T),function(x) x[[i+1]])))
  }
  names(m.meta) <-  metadata
  m.meta <- as.data.frame(m.meta)
  
  # Reorder dataframe if "ref" is defined
  m <- cbind(m, m.subset, m.meta) 
  if (!is.null(ref)) {
    m$m.subset <- factor(m$m.subset, levels = levels(ref$functional.cluster))
    m <- m |> arrange(m.subset)
    
    # Reappend good annotation order
    m.subset <-
      factor(unlist(lapply(strsplit(rownames(m), "!", perl = T), function(x)
        x[[1]])))
  }
  
  # Reorder dataframe if "order.by" is defined
  if (!is.null(order.by)) {
    m <- m |> arrange(m.subset, get(order.by))
  
    # Reappend good annotation order for metadata
    m.meta <- list()
    for (i in 1:length(metadata)) {
      m.meta[[i]] <-
        factor(unlist(lapply(strsplit(rownames(m), "!", perl = T), function(x)
          x[[i + 1]])))
    }
    names(m.meta) <-  metadata
  }
  m <- m[1:(length(m) - length(metadata) - 1)]
  
  
  # Setup color palette list
  
  require(RColorBrewer)
  
  if(palette_reverse){
    color = colorRampPalette(rev(brewer.pal(n = 7, name = brewer.palette)))(length(breaks))
  } else{
    color = colorRampPalette(brewer.pal(n = 7, name = brewer.palette))(length(breaks))  
  }
  
  palettes.default <-  c("Paired","Set2","Accent","Dark2","Set1","Set3")
  if (is.null(palette)) {
    palette <- list()
    palette[["Subtype"]] <- colorRampPalette(brewer.pal(n=8, name="Set1"))(length(unique(unlist(m.subset))))
    names(palette[["Subtype"]]) <- c(unique(m.subset))
    for (i in 1:length(metadata)){
      meta <- metadata[i]
      palette[[meta]] <- colorRampPalette(brewer.pal(n=6, name=palettes.default[i]))(length(unique(m.meta[[meta]])))
      names(palette[[meta]]) <- levels(m.meta[[meta]])
    }
  }
  
  # Define the colors of functional.cluster if ref is given
  if (!is.null(ref)) {
    palette[["Subtype"]] <- ref@misc$atlas.palette
  }
  
  # Compute annotation dataframe
  annotation_col = data.frame(
    Subtype = m.subset
  )
  annotation_col <- cbind(annotation_col, m.meta)
  rownames(annotation_col) = rownames(m)
  
  # Plot heatmap
  if (flip) { 
    h <- pheatmap::pheatmap(m, cluster_rows = cluster_samples,
                            cluster_cols = cluster_genes,scale = scale,
                            breaks = breaks, color=color, 
                            annotation_row = annotation_col, 
                            show_rownames = show_samplenames,
                            border_color = NA,
                            annotation_colors = palette, 
                            fontsize_row=6,fontsize = 7, 
                            clustering_method=method, ...)
  } else {
    h <- pheatmap::pheatmap(t(m),cluster_rows = cluster_genes,
                            cluster_cols = cluster_samples,scale = scale,
                            breaks = breaks, color=color, 
                            annotation_col = annotation_col, 
                            show_colnames = show_samplenames,
                            border_color = NA,
                            annotation_colors = palette, 
                            fontsize_row=6,fontsize = 7, 
                            clustering_method=method, ... )
  }
  return(h)
}

#' Gene expression markers shared by multiple groups of cells
#'
#' This function expands \link[Seurat]{FindAllMarkers} to find markers that are differentially expressed across multiple
#' datasets or samples. Given a Seurat object with identity classes (for example annotated clusters) and a grouping
#' variable (for example a Sample ID), it calculate differentially expressed genes (DEGs) individually for each sample. 
#' Then it determines the fraction of samples for which the gene was found to be differentially expressed.
#' 
#' This function can be useful to find marker genes that are specific for individual cell types, and that are found
#' to be so consistently across multiple samples.
#'
#' @param object A Seurat object
#' @param split.by A metadata column name - the data will be split by this column to calculate \link[Seurat]{FindAllMarkers}
#'     separately for each data split
#' @param only.pos Only return positive markers (TRUE by default)   
#' @param min.cells.group Minimum number of cells in the group - if lower the group is skipped
#' @param min.freq Only return markers which are differentially expressed in at least this fraction of datasets.
#' @param ... Additional paramters to \link[Seurat]{FindAllMarkers}
#' @return A list of marker genes for each identity class (typically clusters), with an associated numerical value.
#'     This value represents the fraction of datasets for which the marker was found to be differentially expressed.
#' @importFrom dplyr filter select
#' @examples
#' ref <- load.reference.map(ref = "https://figshare.com/ndownloader/files/38921366")
#' Idents(ref) <- "functional.cluster"
#' FindAllMarkers.bygroup(ref, split.by = "Sample", min.cells.group=30, min.freq=0.8)
#' @export FindAllMarkers.bygroup

FindAllMarkers.bygroup <- function(object,
                                   split.by = NULL,
                                   only.pos = TRUE,
                                   min.cells.group = 10,
                                   min.freq = 0.5,
                                   ...) {
  if (is.null(split.by)) {
    stop("Please provide a grouping variable with 'split.by' parameter")
  }
  ids <- names(table(Idents(object)))
  meta <- object[[]]
  
  obj.list <- SplitObject(object, split.by = split.by)
  deg <- lapply(obj.list, function(x){
    FindAllMarkers(x, only.pos = only.pos,
                   min.cells.group = min.cells.group, ...)
  })
  
  genes <- lapply(ids, function(i) {
    
    #count groups with at least min cells
    cells <- Idents(object) == i
    t <- table(meta[cells, split.by])
    max.c <- sum(t>=min.cells.group)
    
    #count occurrences
    this <- lapply(deg, function(x) {
      dplyr::filter(x, cluster==i) |> select(gene)
    })
    freqs <- table(unlist(this)) / max.c
    
    #minimum frequency
    freqs <- freqs[freqs >= min.freq]
    sort(freqs, decreasing = T)
  })
  names(genes) <- ids
  
  genes
}