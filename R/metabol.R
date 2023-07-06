
load_metabolic_tasks <- function(tasks_file, info_file){

  rtasks <- read.table(tasks_file, sep="\t", header=T, stringsAsFactors = F)
  rtasks$task_id <- paste0("task", rtasks$task_id)

  # info <- read.table("annotations/tasks_Richelle_cellrep_2021.tsv", sep="\t", header=T, stringsAsFactors = F, quote="\"", comment.char = "", as.is = T)
  # info <- read.table("annotations/richelle__simplify_tasks__clean.csv", sep="\t", header=T, stringsAsFactors = F, quote="\"", comment.char = "", as.is = T)
  info <- read.table(info_file, sep="\t", header=T, stringsAsFactors = F, quote="\"", comment.char = "", as.is = T)
  rownames(info) <- paste0("task_", info$ID)

  gene2task <- split(rtasks$gene_id, rtasks$task_id)
  all_tasked_genes <- unique(rtasks$gene_id)

  gene2task_adj <- do.call("rbind",lapply(gene2task, function(g) {
    x <- numeric(length(all_tasked_genes))
    names(x) <- all_tasked_genes
    x[g] <- 1
    x
  }))

  return(list(
    rtasks=rtasks,
    info=info,
    gene2task=gene2task,
    gene2task_adj=gene2task_adj,
    all_tasked_genes=all_tasked_genes
  ))

}


run_metabolic_tasks_inference <- function(sobj, tasks_file, info_file, verbose=T){

  msj("Loading metabolic tasks")
  mtasks <- load_metabolic_tasks(tasks_file, info_file)

  # norm counts
  norm_counts <- Matrix::as.array(sobj@assays$RNA[,1:ncol(sobj)])
  rownames(norm_counts) <- sapply(strsplit(rownames(norm_counts),"\\."),"[[",1)
  m_norm_counts <- norm_counts[rownames(norm_counts) %in% mtasks$all_tasked_genes,]

  msj("Transforming norm counts")
  z_m_norm_counts <- get_z_score(m_norm_counts)

  msj("Computing tasks activities")
  task_mat <- get_mean_task_matrix(z_m_norm_counts, mtasks$gene2task)

  msj("Initializing seurat object from SBE")
  tobj <- init_seurat_object_from_sbe(task_mat, sobj@meta.data, sobj=sobj)

  return(tobj)

}


run_FBA <- function(sobj, scfea_path){

  msj("Preparing the data")
  norm_counts <- Matrix::as.array(sobj@assays$RNA[,1:ncol(sobj)])
  rownames(norm_counts) <- sapply(strsplit(rownames(norm_counts),"\\."),"[[",1)
  # TO-FIX
  norm_counts <- norm_counts[!duplicated(rownames(norm_counts)), ]

  msj("Running FBA with ScFEA")
  out <- scFEA_wrapper(norm_counts, scfea_path, remove_tmp = F)

  msj("Initializing seurat object from SBE")
  flux_obj <- init_seurat_object_from_sbe(out$flux, sobj@meta.data, sobj=sobj)
  balance_obj <- init_seurat_object_from_sbe(out$balance, sobj@meta.data, sobj=sobj)

  return(list(
    flux=flux_obj,
    balance=balance_obj
  ))

}







