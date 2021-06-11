#Internal function to filter T cells using scGate
filterCells <- function(query.object, species="mouse", gating.model=NULL, ncores=ncores){
  
  if (is.null(gating.model)) {
     gating.model = scGate::scGate_DB[[species]]$Tcell   
  }
  #Note 1: possibly read relevant model from atlas
  
  data(cell.cycle.obj)
  query.object <- scGate(query.object, gating.model = gating.model, ncores=ncores, quiet = TRUE,
                         additional.signatures = cell.cycle.obj[[species]])
  ncells <- ncol(query.object)
  
  ncells.keep <- sum(query.object$is.pure == 'Pure')
  if (ncells.keep > 0) {
     query.object <- subset(query.object, subset=is.pure=='Pure') 
  } else {
     query.object <- NULL
  }
  message <- sprintf("%i out of %i ( %i%% ) non-pure T cells removed. Use filter.cells=FALSE to avoid pre-filtering (NOT RECOMMENDED)",
                     ncells - ncells.keep, ncells, round(100*(ncells-ncells.keep)/ncells))
  print(message)
  
  #Parse metadata columns
  query.object$cycling.score <- query.object$cycling_scGate
  query.object$cycling.score.G1_S <- query.object$cycling_G1.S_scGate
  query.object$cycling.score.G2_M <- query.object$cycling_G2.M_scGate
  
  to_remove <- grep("_Zscore$", colnames(query.object@meta.data), perl=T)
  to_remove <- c(to_remove, grep("_scGate$", colnames(query.object@meta.data), perl=T))
  to_remove <- c(to_remove, which(colnames(query.object@meta.data) %in% c("is.pure","scGate.annotation")))
  
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

#Internal function for mouse-human ortholog conversion
convert.orthologs <- function(obj, table, id="Gene.HS", query.assay="RNA", slot="counts") {
  
  exp.mat <- slot(obj@assays[[query.assay]], name=slot)
  exp.mat <- exp.mat[rownames(exp.mat) %in% table[[id]], ]
  
  mouse.genes <- table$Gene.MM[match(row.names(exp.mat),table[[id]])]
  
  row.names(exp.mat) <- mouse.genes
  slot(obj@assays[[query.assay]], name=slot) <- exp.mat
  return(obj)
}

#Internal function to merge Seurat objects including reductions (PCA, UMAP, ICA)
merge.Seurat.embeddings <- function(x=NULL, y=NULL, ...)
{  
  require(Seurat)
  
  #first regular Seurat merge, inheriting parameters
  m <- merge(x, y, ...)
  #preserve reductions (PCA, UMAP, ...)
  
  reds <- intersect(names(x@reductions), names(y@reductions))
  for (r in reds) {
    message(sprintf("Merging %s embeddings...", r))
    
    m@reductions[[r]] <- x@reductions[[r]]
    if (dim(y@reductions[[r]]@cell.embeddings)[1]>0) {
      m@reductions[[r]]@cell.embeddings <- rbind(m@reductions[[r]]@cell.embeddings, y@reductions[[r]]@cell.embeddings)
    }
    if (dim(y@reductions[[r]]@feature.loadings)[1]>0) {
      m@reductions[[r]]@feature.loadings <- rbind(m@reductions[[r]]@feature.loadings, y@reductions[[r]]@feature.loadings)
    }
    if (dim(y@reductions[[r]]@feature.loadings.projected)[1]>0) {
      m@reductions[[r]]@feature.loadings.projected <- rbind(m@reductions[[r]]@feature.loadings.projected, y@reductions[[r]]@feature.loadings.projected)
    }
  }
  return(m)
  
}

#Helper for projecting individual data sets
projection.helper <- function(query, ref=NULL, filter.cells=T, query.assay=NULL, direct.projection=FALSE,
                              seurat.k.filter=200, skip.normalize=FALSE, id="query1", scGate_model=NULL, ncores=ncores) {
  
  retry.direct <- FALSE
  
  #Reference
  DefaultAssay(ref) <- "integrated"
  ref.var.features <- ref@assays$integrated@var.features
  
  #If query.assay not speficied, use the default
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
  
  #automatically determine gene ID column
  g.mm <- length(intersect(row.names(query), Hs2Mm.convert.table$Gene.MM))
  g.hs1 <- length(intersect(row.names(query), Hs2Mm.convert.table$Gene.stable.ID.HS))
  g.hs2 <- length(intersect(row.names(query), Hs2Mm.convert.table$Gene.HS))
  gg <- c(g.mm, g.hs1, g.hs2)
  
  if (max(gg)==g.mm) {
    human.ortho=FALSE
    species='mouse'
  } else {
    hs.id.col <- ifelse(g.hs1 > g.hs2, "Gene.stable.ID.HS", "Gene.HS")
    if (max(g.hs1, g.hs2)<100) {
      message("Warning: fewer than 100 human genes with orthologs found. Check your matrix format and gene names")
    }
    human.ortho=TRUE
    species='human'
  }
  
  if(filter.cells){
    message("Pre-filtering of T cells with scGate...")
    query <- filterCells(query, species=species, gating.model=scGate_model, ncores=ncores)
  }
  if (is.null(query)) {
    message(sprintf("Warning! Skipping %s - all cells were removed by T cell filter", id))
    return(NULL)
  }
  
  #Check if slots are populated, and normalize data.
  if (skip.normalize) {
    slot <- "data"
    exp.mat <-  slot(query@assays[[query.assay]], name=slot)
    if (dim(exp.mat)[1]==0) {
      stop("Data slot not found in your Seurat object. Please normalize the data")
    }
    if (human.ortho) {
      print("Transforming expression matrix into space of mouse orthologs") 
      query <- convert.orthologs(query, table=Hs2Mm.convert.table, id=hs.id.col, query.assay=query.assay, slot=slot)
    }        
  } else {
    slot <- "counts"
    exp.mat <-  slot(query@assays[[query.assay]], name=slot)
    if (dim(exp.mat)[1]==0) {
      stop("Counts slot not found in your Seurat object. If you already normalized your data, re-run with option skip.normalize=TRUE")
    }
    if (human.ortho) {
      print("Transforming expression matrix into space of mouse orthologs") 
      query <- convert.orthologs(query, table=Hs2Mm.convert.table, id=hs.id.col, query.assay=query.assay, slot=slot)
    }
    query@assays[[query.assay]]@data <- query@assays[[query.assay]]@counts
    query <- NormalizeData(query) 
  }
  rm(exp.mat)
  
  query <- RenameCells(query, add.cell.id = "Q")
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
        
        print(paste0("Aligning ", id, " to reference map for batch-correction..."))
        
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
        query.pca.proj <- apply.pca.obj.2(projected, pca.obj=ref@misc$pca_object, query.assay="integrated")
        projected[["pca"]] <- CreateDimReducObject(embeddings = query.pca.proj, key = "PC_", assay = "integrated")
        
        cat("\nProjecting corrected query onto Reference UMAP space\n")
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
      tryCatch(    #Try Direct projection
        expr = {
          projected <- query
          
          print("DIRECTLY projecting query onto Reference PCA space")
          query.pca.proj <-apply.pca.obj.2(query, pca.obj=ref@misc$pca_object, query.assay=query.assay)
          projected[["pca"]] <- CreateDimReducObject(embeddings = query.pca.proj, key = "PC_", assay = query.assay)
          
          print("DIRECTLY projecting query onto Reference UMAP space")
          query.umap.proj <- make.umap.predict.2(ref.umap=ref@misc$umap_obj, pca.query.emb = query.pca.proj)
          projected[["umap"]] <- CreateDimReducObject(embeddings = query.umap.proj, key = "UMAP_", assay = query.assay)
          
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
  }
  return(projected)
}

