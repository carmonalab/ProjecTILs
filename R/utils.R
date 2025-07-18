filterCells <- function(query.object, species="mouse", gating.model=NULL){

  ncells <- ncol(query.object)
  if (ncells <= 1) {
    return(NULL)
  }
  if (is.null(gating.model)) {
    return(query.object)
  }
  pca.dim <- 30
  ncells <- ncol(query.object)
  if (ncells <= pca.dim) {
    pca.dim <- ncells - 1
  }

  data(cell.cycle.obj)
  query.object <- suppressWarnings(scGate::scGate(data=query.object,
                                          model = gating.model,
                                          pca.dim = pca.dim,
                                          verbose=FALSE,
                                          assay=DefaultAssay(query.object),
                         additional.signatures = cell.cycle.obj[[species]]))

  ncells.keep <- sum(query.object$is.pure == 'Pure')

  message <- sprintf("%i out of %i ( %i%% ) non-pure cells removed. Use filter.cells=FALSE to avoid pre-filtering",
                     ncells - ncells.keep, ncells, round(100*(ncells-ncells.keep)/ncells))
  print(message)

  if (ncells.keep <= 1) {
    return(NULL)
  }

  query.object <- subset(query.object, subset=is.pure=='Pure')

  #Parse metadata columns
  query.object$cycling.score <- query.object$cycling_UCell
  query.object$cycling.score.G1_S <- query.object$cycling_G1.S_UCell
  query.object$cycling.score.G2_M <- query.object$cycling_G2.M_UCell

  to_remove <- grep("is.pure", colnames(query.object@meta.data))
  to_remove <- c(to_remove, grep("_UCell$", colnames(query.object@meta.data), perl=T))

  query.object@meta.data <- query.object@meta.data[,-to_remove]
  return(query.object)
}

#Internal function to randomly split an object into subsets
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

guess_raw_separator <- function(f, sep=c(" ","\t",",")) {

  lines <- readLines(f, n=10)
  if (length(lines) == 0) {
     return(NULL)
  }
  spl <- lapply(sep, grep, x=lines)
  counts <- unlist(lapply(spl, length))
  if (max(counts)==0) {
    return(NULL)
  }
  sep.index <- which(counts==max(counts))[1]
  return(sep[sep.index])

}

#Automatically determine species and gene ID column
get.species <- function(genes, table=Hs2Mm.convert.table) {

  g.mm <- length(intersect(genes, table$Gene.MM))
  g.hs1 <- length(intersect(genes, table$Gene.stable.ID.HS))
  g.hs2 <- length(intersect(genes, table$Gene.HS))
  gg <- c(g.mm, g.hs1, g.hs2)

  if (max(gg)==g.mm) {
    species='mouse'
    col.id <- "Gene.MM"
  } else {
    species='human'
    col.id <- ifelse(g.hs1 > g.hs2, "Gene.stable.ID.HS", "Gene.HS")
  }
  res <- list("species"=species, "col.id"=col.id)
  return(res)
}


#Internal function for mouse-human ortholog conversion
convert.orthologs <- function(obj, table, from="Gene.HS", to="Gene.MM",
                              query.assay="RNA", slot="counts") {

  exp.mat <- GetAssayData(obj, assay=query.assay, layer=slot)
  genes.select <- rownames(exp.mat)[rownames(exp.mat) %in% table[[from]]]

  if (length(genes.select) < 100) {
      message("Warning: fewer than 100 genes with orthologs were found. Check your matrix format and gene names")
  }

  if (length(genes.select) > 0) {
    exp.mat <- exp.mat[genes.select, ]
  } else {
    stop(paste0("Error: No genes found in column ", from))
  }

  #Convert
  ortho.genes <- table[[to]][match(row.names(exp.mat), table[[from]])]

  #Update matrix gene names
  row.names(exp.mat) <- ortho.genes

  #Re-generate object
  if (slot=="counts") {
    this <- CreateAssayObject(counts=exp.mat)
  } else {
    this <- CreateAssayObject(data=exp.mat)
  }
  suppressWarnings(obj[[query.assay]] <- this)
  return(obj)
}

#Helper for projecting individual data sets
projection.helper <- function(query, ref=NULL, filter.cells=TRUE, query.assay=NULL,
                              direct.projection=FALSE, fast.umap.predict=FALSE,
                              ortholog_table=NULL,
                              STACAS.k.weight=100, STACAS.k.anchor=5,
                              STACAS.anchor.coverage=1, STACAS.correction.scale=100,
                              skip.normalize=FALSE, id="query1",
                              alpha=0.5, remove.thr=0,
                              scGate_model=NULL, ncores=1) {

  retry.direct <- FALSE
  do.orthology <- FALSE

  #Reference
  DefaultAssay(ref) <- "integrated"
  ref.var.features <- VariableFeatures(ref)

  #If query.assay not specified, use the default
  if (is.null(query.assay)) {
    query.assay <- DefaultAssay(query)
  } else {
     DefaultAssay(query) <- query.assay
  }
  print(paste0("Using assay ",query.assay," for ",id))
  
  if (!is.null(ref@misc$umap_object$data)) {
     pca.dim=ncol(ref@misc$umap_object$data) #use the number of PCs used to build the reference
  } else {
     pca.dim=10
  }

  species.ref <- get.species(genes=row.names(ref), table=ortholog_table)
  species.query <- get.species(genes=row.names(query), table=ortholog_table)

  if (species.ref$species != species.query$species) {
     do.orthology <- TRUE
  }
  
  #Check if slots are populated, and normalize data.
  if (skip.normalize) {
    gr <- grep("^data", Layers(query))
    if (length(gr) == 0) {
      stop("Data slot not found in your Seurat object. Please normalize the data")
    } else if (length(gr) > 1) {
      query <- JoinLayers(query)
    }
    query <- convert_to_v3(query, assay=query.assay, layer="data")
    
  } else {
    gr <- grep("^counts", Layers(query))
    if (length(gr) == 0) {
      stop("Counts slot not found in your Seurat object. If you already normalized your data, re-run with option skip.normalize=TRUE")
    } else if (length(gr) > 1) {
      query <- JoinLayers(query)
    }
    query <- convert_to_v3(query, assay=query.assay, layer="counts")
    query <- NormalizeData(query)
  }
  
  if(filter.cells){
    message("Pre-filtering cells with scGate...")
    if (is.null(scGate_model)) {  #read filter model from atlas
      if (!is.null(ref@misc$scGate[[species.query$species]])) {
        scGate_model <- ref@misc$scGate[[species.query$species]]
      } else {
        scGate_model <- NULL
        message("No scGate model specified: all cells will be projected")
      }
    }
    query <- filterCells(query, species=species.query$species, gating.model=scGate_model)
  }
  if (is.null(query)) {
    message(sprintf("Warning! Skipping %s - all cells were removed by cell filter", id))
    return(NULL)
  }

  if (do.orthology) {
    print("Transforming expression matrix into space of orthologs")
    query <- convert.orthologs(query, table=ortholog_table, query.assay=query.assay, slot="data",
                               from=species.query$col.id, to=species.ref$col.id)
  }  
    
  query <- RenameCells(query, add.cell.id = "Q")
  query.metadata <- query@meta.data   #back-up metadata (and re-add it after projection)

  genes4integration <- intersect(ref.var.features, row.names(query))

  if(length(genes4integration)/length(ref.var.features)<0.5) {
    stop("Too many genes missing. Check input object format") }
  #TODO implement ID mapping? e.g. from ENSEMBLID to symbol?

  if (length(genes4integration)/length(ref.var.features)<0.8) {
    print("Warning! more than 20% of variable genes not found in the query")
  }

  if (direct.projection) {
    projected <- query

    print("DIRECTLY projecting query onto Reference PCA space")
    query.pca.proj <-apply.pca.obj.2(query, pca.obj=ref@misc$pca_object,
                                     query.assay=query.assay)
    projected[["pca"]] <- CreateDimReducObject(embeddings = query.pca.proj,
                                               key = "PC_", assay = query.assay)

    print("DIRECTLY projecting query onto Reference UMAP space")
    query.umap.proj <- make.umap.predict(ref.umap=ref@misc$umap_obj,
                                         query.assay=query.assay,
                                         pca.query.emb = query.pca.proj,
                                         fast.umap.predict=fast.umap.predict)
    projected[["umap"]] <- CreateDimReducObject(embeddings = query.umap.proj,
                                                key = "UMAP_", assay = query.assay)

    DefaultAssay(projected) <- query.assay
  } else {
    tryCatch(    #Try to do alignment, if it fails (too few cells?) do direct projection
      expr = {

        print(paste0("Aligning ", id, " to reference map for batch-correction..."))

        #for compatibility with older versions of STACAS
        is_x <- 'min.sample.size' %in% names(formals(FindAnchors.STACAS))

        if (is_x) {
          proj.anchors <- FindAnchors.STACAS(object.list = list(ref, query),
                                             assay = c("integrated", query.assay),
                                             anchor.features = genes4integration,
                                             dims = 1:pca.dim, alpha = alpha,
                                             k.anchor = STACAS.k.anchor,
                                             anchor.coverage = STACAS.anchor.coverage,
                                             correction.scale = STACAS.correction.scale,
                                             verbose = FALSE, min.sample.size = 1)
        } else {
          proj.anchors <- FindAnchors.STACAS(object.list = list(ref, query),
                                             assay = c("integrated", query.assay),
                                             anchor.features = genes4integration,
                                             dims = 1:pca.dim, alpha = alpha,
                                             k.anchor = STACAS.k.anchor,
                                             anchor.coverage = STACAS.anchor.coverage,
                                             correction.scale = STACAS.correction.scale,
                                             verbose = FALSE)
        }
        #always integrate query into reference
        tree <- matrix(c(-1,-2), nrow=1, ncol=2)

        projected <- suppressWarnings(IntegrateData.STACAS(proj.anchors, k.weight = STACAS.k.weight,
                             dims=1:pca.dim, sample.tree = tree,
                             features.to.integrate = genes4integration,
                             verbose = FALSE))

        #Subset query data from integrated space
        cells_query <- colnames(query)
        projected <- suppressMessages(subset(projected, cells = cells_query))

        projected@meta.data <- query.metadata

        rm(proj.anchors)

        #Make PCA and UMAP projections
        cat("\nProjecting corrected query onto Reference PCA space\n")
        query.pca.proj <- apply.pca.obj.2(projected,
                                          pca.obj=ref@misc$pca_object,
                                          query.assay="integrated")
        projected[["pca"]] <- CreateDimReducObject(embeddings = query.pca.proj, key = "PC_", assay = "integrated")

        cat("\nProjecting corrected query onto Reference UMAP space\n")
        query.umap.proj <- make.umap.predict(ref.umap=ref@misc$umap_obj,
                                             pca.query.emb=query.pca.proj,
                                             query.assay="integrated",
                                             fast.umap.predict=fast.umap.predict)
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
      tryCatch(    #Try Direct projection
        expr = {
          projected <- query

          print("DIRECTLY projecting query onto Reference PCA space")
          query.pca.proj <- apply.pca.obj.2(query, pca.obj=ref@misc$pca_object,
                                            query.assay=query.assay)
          projected[["pca"]] <- CreateDimReducObject(embeddings = query.pca.proj,
                                                     key = "PC_", assay = query.assay)

          print("DIRECTLY projecting query onto Reference UMAP space")
          query.umap.proj <- make.umap.predict(ref.umap=ref@misc$umap_obj,
                                               pca.query.emb = query.pca.proj,
                                               query.assay=query.assay,
                                               fast.umap.predict=fast.umap.predict)
          projected[["umap"]] <- CreateDimReducObject(embeddings = query.umap.proj,
                                                      key = "UMAP_", assay = query.assay)

          DefaultAssay(projected) <- query.assay
        },
        error = function(e) {
          message(paste("Direct projection failed due to:", e, "\n"))
          message(sprintf("Warning: failed to project dataset %s...", id))
          projected <- NULL
        }
      )
    }
  }

  if (!is.null(projected)) {
      VariableFeatures(projected, assay=query.assay) <- ref.var.features
      cellnames <- gsub("^Q_","",colnames(projected))  #remove prefix from cell names
      projected <- RenameCells(projected, new.names=cellnames)
  }
  return(projected)
}

#Utility to convert Seurat objects from v5 to v3
convert_to_v3 <- function(object, assay="RNA", layer="counts") {
  
  if (inherits(object[[assay]], "Assay5")) {
    if (layer == "data") {
      assay_v3 <- CreateAssayObject(
        data = object[[assay]]$data
      )
    } else {
      assay_v3 <- CreateAssayObject(
        counts = object[[assay]]$counts
      )
    } 
    suppressWarnings(object[[assay]] <- assay_v3)
  }
  object
}

#calculate Silhouette coefficient between for cells in rows compared to set in columns with same labels
silhouette_2sets <- function(dist, labs.x, labs.y) {

  labs.x <- as.character(labs.x)
  labs.y <- as.character(labs.y)

  ids <- sort(unique(c(labs.x, labs.y)))
  k <- length(ids)

  if(k <= 1)  #must give at least two classes
    return(NA)

  if (nrow(dist) != length(labs.x)) {
    stop(sprintf("Distance matrix has %i rows but %i row cluster labels are given", nrow(dist), length(labs.x)))
  }
  if (ncol(dist) != length(labs.y)) {
    stop(sprintf("Distance matrix has %i columns but %i column cluster labels are given", ncol(dist), length(labs.y)))
  }

  res <- data.frame(matrix(NA, nrow(dist), 2, dimnames = list(rownames(dist), c("cluster","sil_width"))))

  for (j in 1:k) {
    lab <- ids[j]
    ix <- labs.x == lab
    iy <- labs.y == lab

    Nx <- sum(ix)
    Ny <- sum(iy)
    Ny.n <- sum(!iy)
    if (Nx > 1) {
      a.i <- rowSums(dist[ix, iy])/Ny
      b.i <- rowSums(dist[ix, !iy])/Ny.n

      s.i <- (b.i - a.i) / pmax(b.i, a.i)

      res[ix, "cluster"] <- lab
      res[ix,"sil_width"] <- s.i
    }
  }
  res
}

#Combine labels from two runs of the classifier to return a consensus label and confidence score
combine_labels_and_confidence <- function(labs1, labs2,
                                          labels.col = "functional.cluster",
                                          labels.col.conf = "functional.cluster.conf") {
  if (is.null(labs1)) {
    return(labs2)
  }
  if (is.null(labs2)) {
    return(labs1)
  }
  l1 <- labs1[[labels.col]]
  names(l1) <- rownames(labs1)
  l2 <- labs2[[labels.col]]
  names(l2) <- rownames(labs2)
  new.labs <- combine_labels(l1, l2)

  c1 <- labs1[[labels.col.conf]]
  names(c1) <- rownames(labs1)
  c2 <- labs2[[labels.col.conf]]
  names(c2) <- rownames(labs2)
  new.conf <- combine_confidence(c1, c2)

  new.conf[is.na(new.labs)] <- NA

  comb <- as.data.frame(new.labs)
  comb[,2] <- new.conf
  colnames(comb) <- c(labels.col, labels.col.conf)
  comb
}

#Combine labels from two runs of the classifier to return a consensus label
combine_labels <- function(labs1, labs2) {

  #No prior labels
  if (is.null(labs1)) {
    consensus <- labs2
    return(consensus)
  } else if (is.null(labs2)) {
    consensus <- labs1
    return(consensus)
  }

  #Combine labels
  comb <- as.data.frame(labs1)
  comb[,"l2"] <- NA
  colnames(comb) <- c("l1","l2")

  comb[names(labs2),"l2"] <- labs2

  consensus <- apply(comb, 1, function(x) {
    if (is.na(x[["l1"]]) & is.na(x[["l2"]])) {
      NA
    } else if (is.na(x[["l1"]]) & !is.na(x[["l2"]])) {
      x[["l2"]]
    } else if (is.na(x[["l2"]]) & !is.na(x[["l1"]])) {
      x[["l1"]]
    } else if (x[["l1"]] == x[["l2"]]) {
      x[["l1"]]
    } else {
      NA
    }
  })
  return(consensus)
}

#Combine labels from two runs of the classifier to return a consensus label
combine_confidence <- function(conf1, conf2) {

  #Combine labels
  comb <- as.data.frame(conf1)
  comb[,"l2"] <- NA
  colnames(comb) <- c("l1","l2")

  comb[names(conf2),"l2"] <- conf2

  consensus <- apply(comb, 1, function(x) {
    if (is.na(x[["l1"]]) & is.na(x[["l2"]])) {
      NA
    } else if (is.na(x[["l1"]]) & !is.na(x[["l2"]])) {
      x[["l2"]]
    } else if (is.na(x[["l2"]]) & !is.na(x[["l1"]])) {
      x[["l1"]]
    } else {
      (x[["l1"]] + x[["l2"]])/2
    }
  })
  return(consensus)
}

#Run ProjecTILs.classifier on a single object
classifier.singleobject <- function(query,
                                    ref=NULL,
                                    filter.cells = TRUE,
                                    reduction="pca",
                                    ndim=NULL, k=5,
                                    nn.decay=0.1,
                                    min.confidence=0.2,
                                    labels.col="functional.cluster",
                                    overwrite=TRUE,
                                    ncores=1,
                                    ...) {
  #UMAP emb. only needed if we want to predict labels based on UMAP neighbors
  if (reduction=="umap") {
    fast.umap.predict <- FALSE
  } else {
    fast.umap.predict <- TRUE
  }

  if(is.list(query)) {
    stop("Query must be a single Seurat object")
  }
  labels.col.conf <- paste0(labels.col, ".conf")

  current.labs <- NULL
  if (labels.col %in% colnames(query[[]])) {
    current.labs <- query[[c(labels.col, labels.col.conf)]]
  }

  query <- make.projection(query=query, ref=ref, filter.cells=filter.cells,
                       fast.umap.predict = fast.umap.predict, ncores=ncores, ...)

  query <- cellstate.predict(ref=ref, query=query,
                          reduction=reduction,
                          ndim=ndim, k=k,
                          nn.decay=nn.decay,
                          min.confidence=min.confidence,
                          labels.col = labels.col)

  #Extract new labels and combine (or overwrite) old labels
  labs <- query[[c(labels.col,labels.col.conf)]]

  if (overwrite) {
    new.labs <- labs
  } else {
    new.labs <- combine_labels_and_confidence(current.labs, labs,
                                              labels.col, labels.col.conf)
  }
  return(new.labs)
}

#Set parallelization options
set_parall <- function(ncores, progressbar=FALSE) {
  if (ncores == 1) {
    param <- SerialParam(progressbar = progressbar)
  } else if (.Platform$OS.type == "windows") {
    param <- SnowParam(workers=ncores, progressbar = progressbar)
  } else {
    param <- MulticoreParam(workers=ncores, progressbar = progressbar)
  }
  return(param)
}


# helper to load rds reference maps
# reference should be a path to a.rds object or a URL to a .rds object, storing a Seurat object prepared using \link{make.reference}
load.helper <- function(reference){
  tryCatch(ref <- readRDS(reference),
           error = function(e){
             stop(paste("Reference object",reference,"is invalid"))
           })
  tryCatch(print(paste0("Loaded Custom Reference map ",ref@misc$projecTILs)),
           error = function(e){stop("Invalid Reference object.\nConsider increasing downloading timeout running:\n `options(timeout = 1000)`\n")
           })
  return(ref)
}


# function to fetch the metadata of figshare entries
get_figshare_metadata <- function(article_id) {

  url <- paste0("https://api.figshare.com/v2/articles/",
                article_id)

  # Make the HTTP GET request
  response <- readLines(url,
                        warn = "F",
                        encoding = "UTF-8")

  # Parse the JSON response
  metadata <- jsonlite::fromJSON(paste(response, collapse = ""))

  # get only url and md5 data
  df <- as.data.frame(metadata$files)

  return(df)
}

# function to donwload object form a url and check integrity
download_integrity <- function(url,
                               destfile,
                               hash = NULL,
                               quiet = F){
  r <- TRUE

  tryCatch({
    download.file(url,
                  destfile = destfile,
                  mode = "wb",
                  quiet = quiet)
  }, error = function(e){
    r <<- FALSE
    file.remove(destfile)
    cat("Download failed for ", destfile,
        "\n Consider increasing downloading timeout running: `options(timeout = 1000)`\n")


  }
  )

  if(r && !is.null(hash)){
    # check file integrity
    downloaded_hash <- digest::digest(file = destfile)
    # return if file integrity check passed
    if(downloaded_hash == hash){
      r <- TRUE
    } else {
      r <- FALSE
    }
  }

  return(r)
}

# function to handle errors during downloading
try.download <- function(url,
                         destfile,
                         hash = NULL,
                         verbose = TRUE,
                         # whether stop function or trown warning upon failing
                         warn = FALSE){


  file_integrity <- download_integrity(url = url,
                                       destfile = destfile,
                                       hash = hash,
                                       quiet = !verbose)
  if(!file_integrity){
    message("File ", destfile, " did not pass integrity check. Redownloading file\nConsider increasing downloading timeout running:\n `options(timeout = 1000)`\n")
    file_integrity <- download_integrity(url = url,
                                        destfile = destfile,
                                        hash = hash,
                                        quiet = !verbose)
  }

  if(!file_integrity){
    if(warn){
      cat("File ", destfile, " did not pass integrity check!!\nConsider increasing downloading timeout running:\n `options(timeout = 1000)`\n")
    } else {
      stop("File ", destfile, " did not pass integrity check!!\nConsider increasing downloading timeout running:\n `options(timeout = 1000)`\n")
    }
  }

}





