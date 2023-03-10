

run_tf_activity_inference <- function(sobj, verbose=T){

  require(dorothea)
  require(Biobase)

  data(dorothea_hs, package = "dorothea")

  msj("Preparing data")
  norm_counts <- Matrix::as.array(sobj@assays$RNA[,1:ncol(sobj)])
  rownames(norm_counts) <- sapply(strsplit(rownames(norm_counts),"\\."),"[[",1)
  # TO-FIX
  norm_counts <- norm_counts[!duplicated(rownames(norm_counts)), ]

  msj("Getting regulon from DoRothEA")
  viper_regulons = df2regulon(dorothea_hs)
  # viper_regulons = df2regulon(grn)

  msj("Running Viper")
  vpres <- viper(norm_counts, viper_regulons, verbose = FALSE)

  # tf_activities <- run_viper(es, dorothea_hs, options =  list(method = "scale", minsize = 4, eset.filter = FALSE, cores = 1, verbose = FALSE))

  msj("Initializing seurat object from SBE")
  sb_obj <- init_seurat_object_from_sbe(vpres, sobj@meta.data, sobj=sobj)

  return(sb_obj)

}

