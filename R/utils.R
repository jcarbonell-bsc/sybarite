
init_seurat_object_from_sbe <- function(mat, meta, sobj=NULL, verbose=F, batch=NULL, label_cell_type=F, nfeatures=nrow(mat)){

  print(dim(mat))

  tobj <- CreateSeuratObject(counts = mat, meta.data = meta, min.cells = 0, min.features = 0, project = "task")
  if(!is.null(batch)){
    tobj@meta.data$batch <- batch
    print(dim(tobj))
    ifnb.list <- SplitObject(tobj, split.by = "batch")
    ifnb.list <- lapply(X = ifnb.list, FUN = function(x) {
      x <- NormalizeData(x)
      x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = nrow(mat))
    })
    features <- SelectIntegrationFeatures(object.list = ifnb.list, nfeatures = nrow(mat))
    immune.anchors <- FindIntegrationAnchors(object.list = ifnb.list, anchor.features = features)
    tobj <- IntegrateData(anchorset = immune.anchors, new.assay.name = "RNA")
    print(dim(tobj))
  }
  tobj <- FindVariableFeatures(object = tobj, selection.method = "vst", nfeatures = nfeatures, verbose=verbose)
  tobj <- ScaleData(object = tobj)
  tobj <- RunPCA(object = tobj, features = VariableFeatures(object = tobj), verbose=verbose)
  # tobj <- RunTSNE(object = tobj, verbose=verbose)
  tobj <- RunUMAP(object = tobj, dims = 1:10, verbose=verbose)
  tobj@meta.data$seurat_clusters_rna <- format_cluster_names(tobj@meta.data$seurat_clusters, "rna_")
  tobj <- FindNeighbors(object = tobj, verbose=verbose)
  tobj <- FindClusters(object = tobj, verbose=verbose)

  if(label_cell_type==T){
    library(celldex)
    # hpca.se <- HumanPrimaryCellAtlasData
    library(SingleR)
    # sboj hpca.se <- HumanPrimaryCellAtlasData
    tobj_se <- as.SingleCellExperiment(tobj)
    # ref <- SingleR::HumanPrimaryCellAtlasData()
    ref <- celldex::HumanPrimaryCellAtlasData()
    pred <- SingleR(test=tobj_se, ref=ref, labels=ref$label.main)
    tobj[["cell_type"]] <- pred$labels
  }

  if(!is.null(sobj)){
    # tobj[["pca_rna"]] <- hobj[["pca_rna"]] <- CreateDimReducObject(embeddings = Embeddings(sobj, "pca"), key = "PC_", assay = DefaultAssay(sobj))
    # tobj[["umap_rna"]] <- hobj[["umap_rna"]] <- CreateDimReducObject(embeddings = Embeddings(sobj, "umap"), key = "UMAP_", assay = DefaultAssay(sobj))
    tobj[["pca_rna"]] <- CreateDimReducObject(embeddings = Embeddings(sobj, "pca"), key = "PC_", assay = DefaultAssay(sobj))
    tobj[["umap_rna"]] <- CreateDimReducObject(embeddings = Embeddings(sobj, "umap"), key = "UMAP_", assay = DefaultAssay(sobj))
  }

  return(tobj)

}

get_quantile_matrix <- function(mat){
  qmat <- mat * 0
  for(i in 1:nrow(mat)){
    qmat[i, ] <- ecdf(mat[i, ])(mat[i, ])
  }
  # t(apply(as.matrix(mat), 1, function(x) ecdf(x)(x)))
  return(qmat)
}

get_quantile_matrix_from_assay <- function(assay){
  get_quantile_matrix(Matrix::as.array(assay[,1:ncol(assay)]))
}


get_z_score <- function(mat, NA_value=0){
  zmat <- t(apply(mat, 1, function(x) (x-mean(x))/sd(x)))
  zmat[is.na(zmat)] <- NA_value
  return(zmat)
}


get_mean_task_matrix <- function(qmat, gene2task){

  tmat <- matrix(0, nrow=length(gene2task), ncol=ncol(qmat))
  rownames(tmat) <- names(gene2task)
  colnames(tmat) <- colnames(qmat)

  for(i in 1:nrow(tmat)){
    # print(i)
    common_genes <- intersect(gene2task[[i]], rownames(qmat))
    if(length(common_genes)>0){
      tmat[i, ] <- apply(qmat[common_genes, , drop=F], 2, mean)
    }
  }

  return(tmat)
}


msj <- function(str, verbose=T){
  if(verbose==T) cat(str, "\n")
}


format_cluster_names <- function(seurat_clusters, preffix=""){
  sc2 <- as.character(seurat_clusters)
  max_char <- max(nchar(sc2))
  ids <- paste0(preffix, gsub(" ", "0", format(seurat_clusters, width=max_char, justify="right")))
  return(ids)
}


translate_gene_ids <- function(norm_counts, xref){

  entrezs <- xref[rownames(norm_counts)]
  clean_entrezs <- sapply(entrezs, function(x) if(is.null(x)) return(NA) else return(x))

  to_remove <- is.na(clean_entrezs) | duplicated(clean_entrezs)

  clean_entrezs2 <- sapply(clean_entrezs[!to_remove],"[[", 1)
  norm_counts2 <- norm_counts[!to_remove, ]
  rownames(norm_counts2) <- clean_entrezs2

  return(norm_counts2)
}

active <- function(sobj){
  Matrix::as.array(sobj@assays[[sobj@active.assay]][,1:ncol(sobj)])
}


load_me <- function(file) {
  out <- new.env()
  load(file,envir = out)
  if(length(out)==1) return(as.list(out)[[1]])
  else return(as.list(out))
}

