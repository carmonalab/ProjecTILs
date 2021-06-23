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
#' @param raw.rawnames For raw matrix format - A vector of row names, or a single number giving the column of the table which contains the row names
#' @param raw.sep For raw matrix format - Separator for raw expression matrix
#' @param raw.header For raw matrix format - Use headers in expression matrix
#' @return A Seurat object populated with raw counts and normalized counts for single-cell expression
#' @examples
#' fname <- "./sample_data"
#' querydata <- read.sc.query(fname, type="10x")
#' @export

read.sc.query <- function(filename, type=c("10x","hdf5","raw","raw.log2"), project.name="Query",
                          min.cells = 3, min.features = 50, gene.column.10x=2,
                          raw.rownames=1, raw.sep=c("auto"," ","\t",","), raw.header=T) {
  
  if (is.null(filename)) {stop("Please provide a query dataset in one of the supported formats")}
  type = tolower(type[1])
  
  if (type == "10x") {
    query.exp <- Read10X(filename, gene.column = gene.column.10x)
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
#' This function allows projecting single-cell RNA-seq datasets onto a reference map of cellular states. To project multiple dataset, submit a list of
#' Seurat objects with the query parameter.
#' 
#' See `load.reference.map()` to load or download a reference atlas.
#'
#' @param query Query data, either as single Seurat object or as a list of Seurat object
#' @param ref Reference Atlas - if NULL, downloads the default TIL reference atlas
#' @param filter.cells Pre-filter T cells using `scGate`. Only set to FALSE if the dataset has been previously subset to T cells only.
#' @param query.assay Which assay slot to use for the query (defaults to DefaultAssay(query))
#' @param direct.projection If true, apply PCA transformation directly without alignment
#' @param seurat.k.filter Integer. For alignment, how many neighbors (k) to use when picking anchors. Default is 200; try lower value in case of failure
#' @param skip.normalize By default, log-normalize the count data. If you have already normalized your data, you can skip normalization.
#' @param scGate_model scGate model to filter T cells from datasets (if NULL use the model stored in \code{ref@@misc$scGate})
#' @param human.ortho Project human data on murine reference atlas, using mouse orthologs (deprecated from v.0.9.9)
#' @param ncores Number of cores for parallel execution (requires \code{future.apply})
#' @param future.maxSize For multi-core functionality, maximum allowed total size (in Mb) of global variables. To increment if required from \code{future.apply}
#' @return An augmented Seurat object with projected UMAP coordinates on the reference map and cell classifications
#' @examples
#' data(query_example_seurat)
#' make.projection(query_example_seurat)
#' @export
make.projection <- function(query, ref=NULL, filter.cells=T, query.assay=NULL, direct.projection=FALSE,
                             seurat.k.filter=200, skip.normalize=FALSE, scGate_model=NULL, human.ortho=FALSE, 
                             ncores=1, future.maxSize=3000) {
   
  
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
  projected.list <- list()
  data(Hs2Mm.convert.table)
  
  if(!is.list(query)) {
     query.list <- list(query=query)
  } else {
    query.list <- query
    if (is.null(names(query.list))) {
       names(query.list) <- paste0("query",c(1:length(query.list)))
    }
  }
  rm(query)
  
  #Parallel or sequential
  if (ncores>1) {
    require(future.apply)
    options(future.globals.maxSize= future.maxSize*1024^2)
    future_ncores <<- ncores
    
    future::plan(future::multisession(workers=future_ncores))
    
    projected.list <- future_lapply(
      X = 1:length(query.list),
      FUN = function(i) {
         res <- projection.helper(query=query.list[[i]], ref=ref, filter.cells=filter.cells, query.assay=query.assay,
                                        direct.projection=direct.projection, seurat.k.filter=seurat.k.filter, ncores=ncores, 
                                        skip.normalize=skip.normalize, id=names(query.list)[i], scGate_model=scGate_model)
         return(res)
      }, future.seed = 1
    )
    plan(strategy = "sequential")
      
  } else {
    projected.list <- lapply(
      X = 1:length(query.list),
      FUN = function(i) {
        res <- projection.helper(query=query.list[[i]], ref=ref, filter.cells=filter.cells, query.assay=query.assay,
                                 direct.projection=direct.projection, seurat.k.filter=seurat.k.filter, ncores=ncores,
                                 skip.normalize=skip.normalize, id=names(query.list)[i], scGate_model=scGate_model)
        return(res)
      }
    )
  }    
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
#' @return UMAP plot of reference map with projected query set in the same space
#' @examples
#' plot.projection(ref, query_example.seurat)
#' @export plot.projection

plot.projection= function(ref, query=NULL, labels.col="functional.cluster", cols=NULL, linesize=1, pointsize=1) {
  require(Seurat)
  require(ggplot2)
  
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
  
  if (is.null(query)) {
    p <- DimPlot(ref, reduction="umap", label = F, group.by = labels.col, repel = T, cols=cols_use) +
      ggtitle ("Reference map") + theme(aspect.ratio=1)
  } else {
    p <- DimPlot(ref, reduction="umap", label = F, group.by = labels.col, repel = T, cols=cols_use) +
      geom_point(data.frame(query@reductions$umap@cell.embeddings), mapping=aes(x=UMAP_1,y=UMAP_2),alpha=0.6, size=pointsize,shape=17, color="gray10") +
      geom_density_2d(data=data.frame(query@reductions$umap@cell.embeddings), mapping=aes(x=UMAP_1,y=UMAP_2),color="black",n=200,h=2,size=linesize) +
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
#' @export plot.statepred.composition
plot.statepred.composition = function(ref, query, labels.col="functional.cluster",cols=NULL, metric=c("Count","Percent")) {
  require(reshape2)
  require(ggplot2)
  
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
      scale_fill_manual(values=cols_use) +
      theme(axis.text.x=element_blank(), legend.position="left")
  } else if (metric=="count") {
    tb.m <- reshape2::melt(tb)
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
#' @param ref Reference Atlas
#' @param query Query data, either as a Seurat object or as a list of Seurat objects
#' @param labels.col The metadata field used to annotate the clusters (default: functional.cluster)
#' @param genes4radar Which genes to use for plotting (default: c("Foxp3","Cd4","Cd8a","Tcf7","Ccr7","Gzmb","Gzmk","Pdcd1","Havcr2","Tox,"Mki67")
#' @param min.cells Only display cell states with a minimum number of cells
#' @param cols Custom color palette for samples in radar plot
#' @param return Return the combined grobs instead of printing it to the default device
#' @param return.as.list Return plots in a list, instead of combining them in a single plot
#' @return Radar plot of gene expression of key genes by cell subtype
#' @examples
#' plot.states.radar(ref)
#' @export plot.states.radar
plot.states.radar = function(ref, query=NULL, labels.col="functional.cluster",
                                  genes4radar=NULL, min.cells=10, cols=NULL, return=F, return.as.list=F) {
  require(ggplot2)
  require(scales)
  require(gridExtra)
  
  #Make sure query is a list
  if(!is.list(query)) {
    query <- list(Query=query)
  }
  
  #Set genes
  if (is.null(genes4radar)) {
    genes4radar <- c("Foxp3","Cd4","Cd8a","Tcf7","Ccr7","Gzmb","Gzmk","Pdcd1","Havcr2","Tox","Mki67")
  }
  
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
  
  genes4radar <- intersect(genes4radar, row.names(ref@assays$RNA@data))
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
      qq[[i]] <- as.matrix(query[[i]]@assays$RNA@data[order,])
    }
  }
  
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
    
    levs <- unique(this.df$Dataset)
    this.df$Dataset <- factor(this.df$Dataset, levels=levs)
    
    pll[[j]] <- ggplot(data=this.df,  aes(x=Gene, y=Expression, group= Dataset, colour=Dataset, fill=Dataset)) +
      geom_point(size=2) +
      geom_polygon(size = 0.75, alpha= 0.1) +
      ylim(ymin, ymax) + ggtitle(s)  +
      scale_x_discrete() +
      scale_fill_manual(values= radar.colors) +
      scale_colour_manual(values= radar.colors) +
      theme_light() +
      theme(axis.text.x=element_blank()) +
      annotate(geom="text", x=seq(1,length(genes4radar)), y=ymax-0.05*ymax, label=genes4radar, size=3) +
      coord_polar()
    
  }
  #Return plots
  if (return.as.list) {
    return(pll)
  } else {
    g <- do.call("arrangeGrob", c(pll, ncol=3, top=paste0("Radar plots for ", labels.col)))
    if(return){
      return(g)
    } else{
      return(plot(g))
    }
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
               colors=cols ) %>% plotly::layout(
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
#' @param all.genes Whether to consider all genes for DE analysis (default is variable genes of the reference)
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
#' @export find.discriminant.genes
#' 
find.discriminant.genes <- function(ref, query, query.control=NULL, query.assay="RNA",
                                    state="largest", labels.col="functional.cluster",
                                    test="wilcox", min.cells=10, all.genes=F, ...)  ##use ellipsis to pass parameters to FindMarkers
  
{  
  require(Seurat)
  
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
  s1 <- subset(query, cells=s1.cells)
  s1$Group <- "Query"
  
  if (!is.null(query.control)) {
    s2 <- subset(query.control, cells=s2.cells)
  } else {
    s2 <- subset(ref, cells=s2.cells)
  }
  s2$Group <- "Control"
  
  s.m <- merge(s1, s2)
  Idents(s.m) <- "Group"
  
  #use all genes or only variable genes from the reference
  if (all.genes) {
    which.genes <- NULL
  } else {
    which.genes <- intersect(ref@assays$integrated@var.features, rownames(s.m))
  } 
  
  markers <- FindMarkers(s.m, slot="data", ident.1="Query", ident.2="Control", only.pos = F, test.use=test, assay=query.assay,
                         features = which.genes, ...)
  
  
  return(markers)
}
