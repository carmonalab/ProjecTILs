# A set of utils functions adapted and simplified from Seurat v4.0.1
# Hao et al. Cell 2021 - https://github.com/satijalab/seurat

ReadMtx.fix <- function(
  mtx,
  cells,
  features,
  cell.column = 1,
  feature.column = 2,
  skip.cell = 0,
  skip.feature = 0,
  unique.features = TRUE,
  strip.suffix = FALSE
) {
  all.files <- list(
    "expression matrix" = mtx,
    "barcode list" = cells,
    "feature list" = features
  )
  for (i in seq_along(along.with = all.files)) {
    all.files[[i]] <- normalizePath(all.files[[i]], mustWork = FALSE)
  }
  
  cell.barcodes <- read.table(
    file = all.files[['barcode list']],
    header = FALSE,
    sep = '\t',
    row.names = NULL,
    skip = skip.cell
  )
  feature.names <- read.table(
    file = all.files[['feature list']],
    header = FALSE,
    sep = '\t',
    row.names = NULL,
    skip = skip.feature
  )
  # read barcodes
  bcols <- ncol(x = cell.barcodes)
  if (bcols < cell.column) {
    stop(
      "cell.column was set to ",
      cell.column,
      " but ",
      cells,
      " only has ",
      bcols,
      " columns.",
      " Try setting the cell.column argument to a value <= to ",
      bcols,
      "."
    )
  }
  cell.names <- cell.barcodes[, cell.column]
  if (all(grepl(pattern = "\\-1$", x = cell.names)) & strip.suffix) {
    cell.names <- as.vector(x = as.character(x = sapply(
      X = cell.names,
      FUN = ExtractField,
      field = 1,
      delim = "-"
    )))
  }
  # read features
  fcols <- ncol(x = feature.names)
  if (fcols < feature.column) {
    stop(
      "feature.column was set to ",
      feature.column,
      " but ",
      features,
      " only has ",
      fcols, " column(s).",
      " Try setting the feature.column argument to a value <= to ",
      fcols,
      "."
    )
  }
  if (any(is.na(x = feature.names[, feature.column]))) {
    na.features <- which(x = is.na(x = feature.names[, feature.column]))
    replacement.column <- ifelse(test = feature.column == 2, yes = 1, no = 2)
    if (replacement.column > fcols) {
      stop(
        "Some features names are NA in column ",
        feature.column,
        ". Try specifiying a different column.",
        call. = FALSE
      )
    } else {
      warning(
        "Some features names are NA in column ",
        feature.column,
        ". Replacing NA names with ID from column ",
        replacement.column,
        ".",
        call. = FALSE
      )
    }
    feature.names[na.features, feature.column] <- feature.names[na.features, replacement.column]
  }
  feature.names <- feature.names[, feature.column]
  if (unique.features) {
    feature.names <- make.unique(names = feature.names)
  }
  data <- readMM(file = all.files[['expression matrix']])
  if (length(x = cell.names) != ncol(x = data)) {
    stop(
      "Matrix has ",
      ncol(data),
      " columns but found ", length(cell.names),
      " barcodes. ",
      ifelse(
        test = length(x = cell.names) > ncol(x = data),
        yes = "Try increasing `skip.cell`. ",
        no = ""
      ),
      call. = FALSE
    )
  }
  if (length(x = feature.names) != nrow(x = data)) {
    stop(
      "Matrix has ",
      ncol(data),
      " rows but found ", length(feature.names),
      " features. ",
      ifelse(
        test = length(x = feature.names) > nrow(x = data),
        yes = "Try increasing `skip.feature`. ",
        no = ""
      ),
      call. = FALSE
    )
  }
  
  colnames(x = data) <- cell.names
  rownames(x = data) <- feature.names
  data <- as(data, Class = "dgCMatrix")
  return(data)
}

#Find integration anchors using reciprocal PCA
FindIntegrationAnchors_local <- function(
    object.list = NULL,
    assay = NULL,
    anchor.coverage = 1,  #level of anchor filtering by distance [0,1]
    correction.scale = 100, #slope of the correction
    alpha=0.5,
    anchor.features = 2000,
    sct.clip.range = NULL,
    l2.norm = TRUE,
    dims = 1:30,
    k.anchor = 5,
    k.filter = NA,
    k.score = 30,
    remove.thr = 0,
    max.features = 200,
    nn.method = "annoy",
    n.trees = 50,
    eps = 0,
    verbose = TRUE
) {
  
  normalization.method <- "LogNormalize"
  reference <- NULL
  reduction <- "pca"
  
  object.ncells <- sapply(X = object.list, FUN = function(x) dim(x = x)[2])
  if (any(object.ncells <= max(dims))) {
    bad.obs <- which(x = object.ncells <= max(dims))
    stop("Max dimension too large: objects ", paste(bad.obs, collapse = ", "),
         " contain fewer than ", max(dims), " cells. \n Please specify a",
         " maximum dimensions that is less than the number of cells in any ",
         "object (", min(object.ncells), ").")
  }
  if (!is.null(x = assay)) {
    if (length(x = assay) != length(x = object.list)) {
      stop("If specifying the assay, please specify one assay per object in the object.list")
    }
    object.list <- sapply(
      X = 1:length(x = object.list),
      FUN = function(x) {
        DefaultAssay(object = object.list[[x]]) <- assay[x]
        return(object.list[[x]])
      }
    )
  } else {
    assay <- sapply(X = object.list, FUN = DefaultAssay)
  }
  object.list <- CheckDuplicateCellNames_local(object.list = object.list)
  
  slot <- "data"
  
  nn.reduction <- reduction
  internal.neighbors <- list()
  
  if (verbose) {
    message("Computing within dataset neighborhoods")
  }
  k.neighbor <- max(k.anchor, k.score)
  internal.neighbors <- lapply(
    X = 1:length(x = object.list),
    FUN = function(x) {
      Seurat:::NNHelper(
        data = Embeddings(object = object.list[[x]][[nn.reduction]])[, dims],
        k = k.neighbor + 1,
        method = nn.method,
        n.trees = n.trees,
        eps = eps
      )
    }
  )
  # determine the proper offsets for indexing anchors
  objects.ncell <- sapply(X = object.list, FUN = ncol)
  offsets <- as.vector(x = cumsum(x = c(0, objects.ncell)))[1:length(x = object.list)]
  
  if (verbose) {
    message("Finding all pairwise anchors")
  }
  
  i <- 1
  j <- 2
  object.1 <- DietSeurat(
    object = object.list[[i]],
    assays = assay[i],
    features = anchor.features,
    counts = FALSE,
    scale.data = TRUE,
    dimreducs = reduction
  )
  object.2 <- DietSeurat(
    object = object.list[[j]],
    assays = assay[j],
    features = anchor.features,
    counts = FALSE,
    scale.data = TRUE,
    dimreducs = reduction
  )
  # suppress key duplication warning
  suppressWarnings(object.1[["ToIntegrate"]] <- object.1[[assay[i]]])
  DefaultAssay(object = object.1) <- "ToIntegrate"
  if (reduction %in% Reductions(object = object.1)) {
    slot(object = object.1[[reduction]], name = "assay.used") <- "ToIntegrate"
  }
  object.1 <- DietSeurat(object = object.1, 
                         assays = "ToIntegrate",
                         counts = FALSE,
                         scale.data = TRUE, 
                         dimreducs = reduction)
  suppressWarnings(object.2[["ToIntegrate"]] <- object.2[[assay[j]]])
  
  DefaultAssay(object = object.2) <- "ToIntegrate"
  if (reduction %in% Reductions(object = object.2)) {
    slot(object = object.2[[reduction]], name = "assay.used") <- "ToIntegrate"
  }
  object.2 <- DietSeurat(object = object.2, 
                         assays = "ToIntegrate",
                         counts = FALSE,
                         scale.data = TRUE, 
                         dimreducs = reduction)
  
  #Reciprocal PCA
  common.features <- intersect(
    x = rownames(x = Loadings(object = object.1[["pca"]])),
    y = rownames(x = Loadings(object = object.2[["pca"]]))
  )
  common.features <- intersect(
    x = common.features,
    y = anchor.features
  )
  object.pair <- merge(x = object.1, y = object.2, merge.data = TRUE)
  projected.embeddings.1<- t(x = GetAssayData(object = object.1, slot = "scale.data")[common.features, ]) %*%
    Loadings(object = object.2[["pca"]])[common.features, ]
  object.pair[['projectedpca.1']] <- CreateDimReducObject(
    embeddings = rbind(projected.embeddings.1, Embeddings(object = object.2[["pca"]])),
    assay = DefaultAssay(object = object.1),
    key = "projectedpca1_"
  )
  projected.embeddings.2 <- t(x = GetAssayData(object = object.2, slot = "scale.data")[common.features, ]) %*%
    Loadings(object = object.1[["pca"]])[common.features, ]
  object.pair[['projectedpca.2']] <- CreateDimReducObject(
    embeddings = rbind(projected.embeddings.2, Embeddings(object = object.1[["pca"]])),
    assay = DefaultAssay(object = object.2),
    key = "projectedpca2_"
  )
  object.pair[["pca"]] <- CreateDimReducObject(
    embeddings = rbind(
      Embeddings(object = object.1[["pca"]]),
      Embeddings(object = object.2[["pca"]])),
    assay = DefaultAssay(object = object.1),
    key = "pca_"
  )
  reduction <- "projectedpca.1"
  reduction.2 <- "projectedpca.2"
  if (l2.norm){
    slot(object = object.pair[["projectedpca.1"]], name = "cell.embeddings") <- Sweep_local(
      x = Embeddings(object = object.pair[["projectedpca.1"]]),
      MARGIN = 2,
      STATS = apply(X = Embeddings(object = object.pair[["projectedpca.1"]]), MARGIN = 2, FUN = sd),
      FUN = "/"
    )
    slot(object = object.pair[["projectedpca.2"]], name = "cell.embeddings") <- Sweep_local(
      x = Embeddings(object = object.pair[["projectedpca.2"]]),
      MARGIN = 2,
      STATS = apply(X = Embeddings(object = object.pair[["projectedpca.2"]]), MARGIN = 2, FUN = sd),
      FUN = "/"
    )
    object.pair <- L2Dim(object = object.pair, reduction = "projectedpca.1")
    object.pair <- L2Dim(object = object.pair, reduction = "projectedpca.2")
    reduction <- paste0(reduction, ".l2")
    reduction.2 <- paste0(reduction.2, ".l2")
  }
  
  internal.neighbors <- internal.neighbors[c(i, j)]
  
  anchors <- FindAnchors_local(
    object.pair = object.pair,
    assay = c("ToIntegrate", "ToIntegrate"),
    slot = slot,
    cells1 = colnames(x = object.1),
    cells2 = colnames(x = object.2),
    internal.neighbors = internal.neighbors,
    reduction = reduction,
    reduction.2 = reduction.2,
    nn.reduction = nn.reduction,
    dims = dims,
    k.anchor = k.anchor,
    k.filter = k.filter,
    k.score = k.score,
    max.features = max.features,
    nn.method = nn.method,
    n.trees = n.trees,
    eps = eps,
    verbose = verbose
  )
  anchors[, 1] <- anchors[, 1] + offsets[i]
  anchors[, 2] <- anchors[, 2] + offsets[j]
  
  #Average distances
  anchors <- as.data.frame(anchors)
  anchors$dist.mean <- apply(anchors[,c("dist1.2","dist2.1")], MARGIN=1, mean)
  message(sprintf("  SD on anchor distances: %.3f",sd(anchors$dist.mean)))
  
  if (anchor.coverage < 1) {

    #Combine anchor distance with anchor score
    sigmoid_center <- unname(quantile(anchors$dist.mean, probs = anchor.coverage, na.rm = T))
    
    distance_factors <-  sigmoid(x = anchors$dist.mean, center = sigmoid_center, scale = correction.scale)
    
    #anchors$score <- alpha*distance_factors + (1-alpha)*anchors$score
    
    #Multiply distance factors by score
    anchors$score <- anchors$score * distance_factors
    
    ##Remove distant anchors
    anchors <- anchors[distance_factors > remove.thr,] 
    
  }
  nanchors <- nrow(anchors)
  #message(sprintf("    Retaining %i anchors after filtering by rPCA distance", nanchors))
  
  ##Include reciprocal anchors
  anchors <- rbind(anchors[, c("cell1","cell2","score","dist.mean")],
                   anchors[, c("cell2","cell1","score","dist.mean")])  
  anchors <- AddDatasetID_local(anchor.df = anchors, offsets = offsets, obj.lengths = objects.ncell)
  
  command <- LogSeuratCommand(object = object.list[[1]], return.command = TRUE)
  anchor.set <- new(Class = "IntegrationAnchorSet",
                    object.list = object.list,
                    reference.objects = seq_along(object.list),
                    anchors = anchors,
                    offsets = offsets,
                    anchor.features = anchor.features,
                    command = command
  )
  
  return(anchor.set)
}

sigmoid <- function(x, scale, center){
  sigm <- 1/(1 + exp(scale*(x-center)))
  return(sigm)
}

#Add dataset ID
AddDatasetID_local <- function(
    anchor.df,
    offsets,
    obj.lengths
) {
  ndataset <- length(x = offsets)
  row.offset <- rep.int(x = offsets, times = obj.lengths)
  dataset <- rep.int(x = 1:ndataset, times = obj.lengths)
  
  anchor.df <- data.frame(
    'cell1' = anchor.df[, 'cell1'] - row.offset[anchor.df[, 'cell1']],
    'cell2' = anchor.df[, 'cell2'] - row.offset[anchor.df[, 'cell2']],
    'score' = anchor.df[, 'score'],
    'dataset1' = dataset[anchor.df[, 'cell1']],
    'dataset2' = dataset[anchor.df[, 'cell2']],
    'dist.mean' = anchor.df[, 'dist.mean']
  )
  return(anchor.df)
}

#Find anchors between a pair of objects
FindAnchors_local <- function(
  object.pair,
  assay,
  slot,
  cells1,
  cells2,
  internal.neighbors,
  reduction,
  reduction.2 = character(),
  nn.reduction = reduction,
  dims = 1:10,
  k.anchor = 5,
  k.filter = NA,
  k.score = 30,
  max.features = 200,
  nn.method = "annoy",
  n.trees = 50,
  nn.idx1 = NULL,
  nn.idx2 = NULL,
  eps = 0,
  verbose = TRUE
) {
  # compute local neighborhoods, use max of k.anchor and k.score if also scoring to avoid
  # recomputing neighborhoods
  k.neighbor <- k.anchor
  if (!is.na(x = k.score)) {
    k.neighbor <- max(k.anchor, k.score)
  }
  object.pair <- FindNN_local(
    object = object.pair,
    cells1 = cells1,
    cells2 = cells2,
    internal.neighbors = internal.neighbors,
    dims = dims,
    reduction = reduction,
    reduction.2 = reduction.2,
    nn.reduction = nn.reduction,
    k = k.neighbor,
    nn.method = nn.method,
    n.trees = n.trees,
    nn.idx1 = nn.idx1,
    nn.idx2 = nn.idx2,
    eps = eps,
    verbose = verbose
  )
  object.pair <- FindAnchorPairs_local(
    object = object.pair,
    integration.name = "integrated",
    k.anchor = k.anchor,
    verbose = verbose
  )
  if (!is.na(x = k.score)) {
    object.pair = ScoreAnchors_local(
      object = object.pair,
      assay = DefaultAssay(object = object.pair),
      integration.name = "integrated",
      verbose = verbose,
      k.score = k.score
    )
  }
  
  ###Return distances
  anc.tab <- object.pair@tools$integrated@anchors
  d1.2 <- numeric(length = dim(anc.tab)[1])
  d2.1 <- numeric(length = dim(anc.tab)[1])
  for (r in 1:dim(anc.tab)[1]) {
    c1 <- anc.tab[r,"cell1"]
    c2 <- anc.tab[r,"cell2"]
    d1.2[r] <- object.pair@tools$integrated@neighbors$nnab@nn.dist[c1, which(object.pair@tools$integrated@neighbors$nnab@nn.idx[c1,] == c2 )]
    d2.1[r] <- object.pair@tools$integrated@neighbors$nnba@nn.dist[c2, which(object.pair@tools$integrated@neighbors$nnba@nn.idx[c2,] == c1 )]
  }
  
  object.pair@tools$integrated@anchors <- cbind(object.pair@tools$integrated@anchors, dist1.2=d1.2)
  object.pair@tools$integrated@anchors <- cbind(object.pair@tools$integrated@anchors, dist2.1=d2.1)
  
  anchors <- GetIntegrationData(
    object = object.pair,
    integration.name = 'integrated',
    slot = 'anchors'
  )
  return(anchors)
}

#Find anchor pairs
FindAnchorPairs_local <- function(
  object,
  integration.name = 'integrated',
  k.anchor = 5,
  verbose = TRUE
) {
  neighbors <- GetIntegrationData(object = object, integration.name = integration.name, slot = 'neighbors')
  max.nn <- c(ncol(x = neighbors$nnab), ncol(x = neighbors$nnba))
  if (any(k.anchor > max.nn)) {
    message(paste0('warning: requested k.anchor = ', k.anchor, ', only ', min(max.nn), ' in dataset'))
    k.anchor <- min(max.nn)
  }
  if (verbose) {
    message("Finding anchors")
  }
  # convert cell name to neighbor index
  nn.cells1 <- neighbors$cells1
  nn.cells2 <- neighbors$cells2
  cell1.index <-  suppressWarnings(which(colnames(x = object) == nn.cells1, arr.ind = TRUE))
  ncell <- 1:nrow(x = neighbors$nnab)
  ncell <- ncell[ncell %in% cell1.index]
  anchors <- list()
  # pre allocate vector
  anchors$cell1 <- rep(x = 0, length(x = ncell) * 5)
  anchors$cell2 <- anchors$cell1
  anchors$score <- anchors$cell1 + 1
  idx <- 0
  indices.ab <- Indices(object = neighbors$nnab)
  indices.ba <- Indices(object = neighbors$nnba)
  for (cell in ncell) {
    neighbors.ab <- indices.ab[cell, 1:k.anchor]
    mutual.neighbors <- which(
      x = indices.ba[neighbors.ab, 1:k.anchor, drop = FALSE] == cell,
      arr.ind = TRUE
    )[, 1]
    for (i in neighbors.ab[mutual.neighbors]){
      idx <- idx + 1
      anchors$cell1[idx] <- cell
      anchors$cell2[idx] <- i
      anchors$score[idx] <- 1
    }
  }
  anchors$cell1 <- anchors$cell1[1:idx]
  anchors$cell2 <- anchors$cell2[1:idx]
  anchors$score <- anchors$score[1:idx]
  anchors <- t(x = do.call(what = rbind, args = anchors))
  anchors <- as.matrix(x = anchors)
  object <- SetIntegrationData(
    object = object,
    integration.name = integration.name,
    slot = 'anchors',
    new.data = anchors
  )
  if (verbose) {
    message(paste0("\tFound ", nrow(x = anchors), " anchors"))
  }
  return(object)
}

#Calculate top feautures across a set of dimensions
TopDimFeatures_local <- function(
  object,
  reduction,
  dims = 1:10,
  features.per.dim = 100,
  max.features = 200,
  projected = FALSE
) {
  dim.reduction <- object[[reduction]]
  max.features <- max(length(x = dims) * 2, max.features)
  num.features <- sapply(X = 1:features.per.dim, FUN = function(y) {
    length(x = unique(x = as.vector(x = sapply(X = dims, FUN = function(x) {
      unlist(x = TopFeatures(object = dim.reduction, dim = x, nfeatures = y, balanced = TRUE, projected = projected))
    }))))
  })
  max.per.pc <- which.max(x = num.features[num.features < max.features])
  features <- unique(x = as.vector(x = sapply(X = dims, FUN = function(x) {
    unlist(x = TopFeatures(object = dim.reduction, dim = x, nfeatures = max.per.pc, balanced = TRUE, projected = projected))
  })))
  features <- unique(x = features)
  return(features)
}

#Score anchors
ScoreAnchors_local <- function(
  object,
  assay = NULL,
  integration.name = 'integrated',
  verbose = TRUE,
  k.score = 30
) {
  if (is.null(assay)) {
    assay <- DefaultAssay(object)
  }
  anchor.df <- as.data.frame(x = GetIntegrationData(object = object, integration.name = integration.name, slot = 'anchors'))
  neighbors <- GetIntegrationData(object = object, integration.name = integration.name, slot = "neighbors")
  offset <- length(x = neighbors$cells1)
  indices.aa <- Indices(object = neighbors$nnaa)
  indices.bb <- Indices(object = neighbors$nnbb)
  indices.ab <- Indices(object = neighbors$nnab)
  indices.ba <- Indices(object = neighbors$nnba)
  nbrsetA <- function(x) c(indices.aa[x, 1:k.score], indices.ab[x, 1:k.score] + offset)
  nbrsetB <- function(x) c(indices.ba[x, 1:k.score], indices.bb[x, 1:k.score] + offset)
  # score = number of shared neighbors
  anchor.new <- data.frame(
    'cell1' = anchor.df[, 1],
    'cell2' = anchor.df[, 2],
    'score' = mapply(
      FUN = function(x, y) {
        length(x = intersect(x = nbrsetA(x = x), nbrsetB(x = y)))},
      anchor.df[, 1],
      anchor.df[, 2]
    )
  )
  # normalize the score
  max.score <- quantile(anchor.new$score, 0.9)
  min.score <- quantile(anchor.new$score, 0.01)
  anchor.new$score <- anchor.new$score - min.score
  anchor.new$score <- anchor.new$score / (max.score - min.score)
  anchor.new$score[anchor.new$score > 1] <-  1
  anchor.new$score[anchor.new$score < 0] <- 0
  anchor.new <- as.matrix(x = anchor.new)
  object <- SetIntegrationData(
    object = object,
    integration.name = integration.name,
    slot = 'anchors',
    new.data = anchor.new
  )
  return(object)
}

#Ensure no duplicate cell names
CheckDuplicateCellNames_local <- function(object.list, verbose = TRUE, stop = FALSE) {
  cell.names <- unlist(x = lapply(X = object.list, FUN = colnames))
  if (any(duplicated(x = cell.names))) {
    if (stop) {
      stop("Duplicate cell names present across objects provided.")
    }
    if (verbose) {
      warning("Some cell names are duplicated across objects provided. Renaming to enforce unique cell names.")
    }
    object.list <- lapply(
      X = 1:length(x = object.list),
      FUN = function(x) {
        return(RenameCells(
          object = object.list[[x]],
          new.names = paste0(Cells(x = object.list[[x]]), "_", x)
        ))
      }
    )
  }
  return(object.list)
}

# Find nearest neighbors
FindNN_local <- function(
  object,
  cells1 = NULL,
  cells2 = NULL,
  internal.neighbors,
  grouping.var = NULL,
  dims = 1:10,
  reduction = "cca.l2",
  reduction.2 = character(),
  nn.dims = dims,
  nn.reduction = reduction,
  k = 300,
  nn.method = "annoy",
  n.trees = 50,
  nn.idx1 = NULL,
  nn.idx2 = NULL,
  eps = 0,
  integration.name = 'integrated',
  verbose = TRUE
) {
  if (xor(x = is.null(x = cells1), y = is.null(x = cells2))) {
    stop("cells1 and cells2 must both be specified")
  }
  if (!is.null(x = cells1) && !is.null(x = cells2) && !is.null(x = grouping.var)) {
    stop("Specify EITHER grouping.var or cells1/2.")
  }
  if (is.null(x = cells1) && is.null(x = cells2) && is.null(x = grouping.var)) {
    stop("Please set either cells1/2 or grouping.var")
  }
  if (!is.null(x = grouping.var)) {
    if (nrow(x = unique(x = object[[grouping.var]])) != 2) {
      stop("Number of groups in grouping.var not equal to 2.")
    }
    groups <- names(x = sort(x = table(object[[grouping.var]]), decreasing = TRUE))
    cells1 <- colnames(x = object)[object[[grouping.var]] == groups[[1]]]
    cells2 <- colnames(x = object)[object[[grouping.var]] == groups[[2]]]
  }
  if (verbose) {
    message("Finding neighborhoods")
  }
  dim.data.self <- Embeddings(object = object[[nn.reduction]])[, nn.dims]
  if (!is.null(x = internal.neighbors[[1]])) {
    nnaa <- internal.neighbors[[1]]
  } else {
    dims.cells1.self <- dim.data.self[cells1, ]
    nnaa <- Seurat:::NNHelper(
      data = dims.cells1.self,
      k = k + 1,
      method = nn.method,
      n.trees = n.trees,
      eps = eps,
      index = nn.idx1
    )
  }
  if (!is.null(x = internal.neighbors[[2]])) {
    nnbb <- internal.neighbors[[2]]
  } else {
    dims.cells2.self <- dim.data.self[cells2, ]
    nnbb <- Seurat:::NNHelper(
      data = dims.cells2.self,
      k = k + 1,
      method = nn.method,
      n.trees = n.trees,
      eps = eps,
      index = nn.idx1
    )
  }
  if (length(x = reduction.2) > 0) {
    nnab <- Seurat:::NNHelper(
      data = Embeddings(object = object[[reduction.2]])[cells2, ],
      query = Embeddings(object = object[[reduction.2]])[cells1, ],
      k = k,
      method = nn.method,
      n.trees = n.trees,
      eps = eps,
      index = nn.idx2
    )
    nnba <- Seurat:::NNHelper(
      data = Embeddings(object = object[[reduction]])[cells1, ],
      query = Embeddings(object = object[[reduction]])[cells2, ],
      k = k,
      method = nn.method,
      n.trees = n.trees,
      eps = eps,
      index = nn.idx1
    )
  } else {
    dim.data.opposite <- Embeddings(object = object[[reduction]])[ ,dims]
    dims.cells1.opposite <- dim.data.opposite[cells1, ]
    dims.cells2.opposite <- dim.data.opposite[cells2, ]
    nnab <- Seurat:::NNHelper(
      data = dims.cells2.opposite,
      query = dims.cells1.opposite,
      k = k,
      method = nn.method,
      n.trees = n.trees,
      eps = eps,
      index = nn.idx2
    )
    nnba <- Seurat:::NNHelper(
      data = dims.cells1.opposite,
      query = dims.cells2.opposite,
      k = k,
      method = nn.method,
      n.trees = n.trees,
      eps = eps,
      index = nn.idx1
    )
  }
  object <- SetIntegrationData(
    object = object,
    integration.name = integration.name,
    slot = 'neighbors',
    new.data = list('nnaa' = nnaa, 'nnab' = nnab, 'nnba' = nnba, 'nnbb' = nnbb, 'cells1' = cells1, 'cells2' = cells2)
  )
  return(object)
}

Sweep_local <- function(x, MARGIN, STATS, FUN = '-', check.margin = TRUE, ...) {
  if (any(grepl(pattern = 'X', x = names(x = formals(fun = sweep))))) {
    return(sweep(
      X = x,
      MARGIN = MARGIN,
      STATS = STATS,
      FUN = FUN,
      check.margin = check.margin,
      ...
    ))
  } else {
    return(sweep(
      x = x,
      MARGIN = MARGIN,
      STATS = STATS,
      FUN = FUN,
      check.margin = check.margin,
      ...
    ))
  }
}

