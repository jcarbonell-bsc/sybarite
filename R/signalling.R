
run_signalling_inference <- function(sobj, verbose=T, hipathia_pathways_file="annots/pathways.RData", xref_file="annots/xref.RData"){


  norm_counts <- Matrix::as.array(sobj@assays$RNA[,1:ncol(sobj)])
  rownames(norm_counts) <- sapply(strsplit(rownames(norm_counts),"\\."),"[[",1)

  res <- run_hipathia(norm_counts, pathways_file = hipathia_pathways_file, xref_file = xref_file)

  path_vals <- get_paths_data(res, matrix = TRUE)
  path_vals <- path_vals[apply(path_vals, 1, function(x) !all(x==x[1])),]

  print(dim(path_vals))
  msj("Initializing seurat object from SBE")
  sb_obj <- init_seurat_object_from_sbe(path_vals, sobj@meta.data, sobj=sobj)

  return(sb_obj)

}

