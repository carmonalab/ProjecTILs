#Internal function to filter cells using scGate
filterCells <- function(query.object, species="mouse", gating.model=NULL){
  
  ncells <- ncol(query.object)
  if (ncells <= 1) {
    return(NULL)
  }
  if (is.null(gating.model)) {
    return(query.object)
  }
  
  data(cell.cycle.obj)
  query.object <- suppressWarnings(scGate::scGate(data=query.object,
                                          model = gating.model,
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
convert.orthologs <- function(obj, table, from="Gene.HS", to="Gene.MM", query.assay="RNA", slot="counts") {
  
  exp.mat <- slot(obj@assays[[query.assay]], name=slot)
  genes.select <- rownames(exp.mat) %in% table[[from]]
  
  if (length(genes.select) < 100) {
      message("Warning: fewer than 100 genes with orthologs were found. Check your matrix format and gene names")
  }
  
  if (length(genes.select) > 0) {
    exp.mat <- exp.mat[rownames(exp.mat) %in% table[[from]], ]
  } else {
    stop(paste0("Error: No genes found in column ", from))
  }
  
  #Convert
  ortho.genes <- table[[to]][match(row.names(exp.mat), table[[from]])]
  
  #Update matrix gene names
  row.names(exp.mat) <- ortho.genes
  slot(obj@assays[[query.assay]], name=slot) <- exp.mat
  
  if (slot=="data") {  #keep same size of data matrices
    slot(obj@assays[[query.assay]], name="counts") <- exp.mat
  }
  
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
                              scGate_model=NULL, ncores=ncores) {
  
  retry.direct <- FALSE
  do.orthology <- FALSE
  
  #Reference
  DefaultAssay(ref) <- "integrated"
  ref.var.features <- ref@assays$integrated@var.features
  
  #If query.assay not specified, use the default
  if (is.null(query.assay)) {
    query.assay <- DefaultAssay(query)
  } else {
     DefaultAssay(query) <- query.assay
  }
  
  print(paste0("Using assay ",query.assay," for ",id))
  
  if (!is.null(ref@misc$umap_object$data)) { 
     pca.dim=dim(ref@misc$umap_object$data)[2] #use the number of PCs used to build the reference
  } else {
     pca.dim=10
  }
  
  species.ref <- get.species(genes=row.names(ref), table=ortholog_table)
  species.query <- get.species(genes=row.names(query), table=ortholog_table)
  
  if (species.ref$species != species.query$species) {
     do.orthology <- TRUE
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
  
  #Check if slots are populated, and normalize data.
  if (skip.normalize) {
    slot <- "data"
    exp.mat <-  slot(query@assays[[query.assay]], name=slot)
    if (dim(exp.mat)[1]==0) {
      stop("Data slot not found in your Seurat object. Please normalize the data")
    }
    if (do.orthology) {
      print("Transforming expression matrix into space of orthologs") 
      query <- convert.orthologs(query, table=ortholog_table, query.assay=query.assay, slot=slot,
                                 from=species.query$col.id, to=species.ref$col.id)
    }        
  } else {
    slot <- "counts"
    exp.mat <-  slot(query@assays[[query.assay]], name=slot)
    if (dim(exp.mat)[1]==0) {
      stop("Counts slot not found in your Seurat object. If you already normalized your data, re-run with option skip.normalize=TRUE")
    }
    if (do.orthology) {
      print("Transforming expression matrix into space of orthologs") 
      query <- convert.orthologs(query, table=ortholog_table, query.assay=query.assay, slot=slot,
                                 from=species.query$col.id, to=species.ref$col.id)
    }
    query@assays[[query.assay]]@data <- query@assays[[query.assay]]@counts
    query <- NormalizeData(query) 
  }
  rm(exp.mat)
  
  query <- RenameCells(query, add.cell.id = "Q")
  query.metadata <- query@meta.data   #back-up metadata (and re-add it after projection)
  
  genes4integration <- intersect(ref.var.features, row.names(query))
  
  if(length(genes4integration)/length(ref.var.features)<0.5){ stop("Too many genes missing. Check input object format") }
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
        
        proj.integrated <- IntegrateData.STACAS(proj.anchors, k.weight = STACAS.k.weight,
                             dims=1:pca.dim, sample.tree = tree,
                             features.to.integrate = genes4integration,
                             verbose = FALSE)
        
        #Subset query data from integrated space
        cells_query <- colnames(query)
        projected <- subset(proj.integrated, cells = cells_query)
        
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
      projected@assays[[query.assay]]@var.features <- ref.var.features
      cellnames <- gsub("^Q_","",colnames(projected))  #remove prefix from cell names
      projected <- RenameCells(projected, new.names=cellnames)
  }
  return(projected)
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

#Combine labels from two runs of the classifier to return a consensus label
combine_labels <- function(labs1, labs2) {
  
  #No prior labels
  if (is.null(labs1)) {
    consensus <- labs2[,1]
    names(consensus) <- rownames(labs2)
    return(consensus)
  } else if (is.null(labs2)) {
    consensus <- labs1[,1]
    names(consensus) <- rownames(labs1)
    return(consensus)
  }  
  
  #Combine labels
  comb <- as.data.frame(labs1)
  comb[,"l2"] <- NA 
  colnames(comb) <- c("l1","l2")
  
  comb[rownames(labs2),"l2"] <- labs2
  
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
