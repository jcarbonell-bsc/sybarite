
run_enrichment_inference <- function(sobj, gene_sets, verbose=T, method="AUCell"){
  
  msj("Loading metabolic tasks")

  # norm counts
  norm_counts <- Matrix::as.array(sobj@assays$RNA[,1:ncol(sobj)])
  rownames(norm_counts) <- sapply(strsplit(rownames(norm_counts),"\\."),"[[",1)
  
  msj("Estimating enrichment")
  if(method=="AUCell"){
    msj("   using AUCell method...")
    require(AUCell)
    auc_res <- AUCell_run(norm_counts, geneSets = gene_sets)
    mat <- getAUC(auc_res)
  } else {
    msj("   using escape method...")
    require(escape)
    mat <- t(enrichIt(obj = norm_counts, gene.sets = gene_sets, groups = 1000, cores = 2, min.size = 1))
  }
  
  msj("Initializing seurat object from SBE")
  tobj <- init_seurat_object_from_sbe(mat, sobj@meta.data, sobj=sobj)
  
  return(tobj)
  
}