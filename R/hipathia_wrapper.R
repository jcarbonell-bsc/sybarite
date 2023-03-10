
run_hipathia <- function(mat, pathways_file, xref_file){

  load(xref_file)
  t_mat <- translate_ids(mat, xref)

  t_mat <- t_mat/max(t_mat)
  t_mat <- t_mat*.9 + 0.05

  load(pathways_file)
  # pathways$pathigraphs <- pathways$pathigraphs[1:5]
  results <- hipathia(t_mat, pathways, decompose = FALSE, verbose=T)

  return(results)
}


translate_ids <- function(norm_counts, xref){

  entrezs <- xref$from_hgnc_to_entrez[rownames(norm_counts)]
  clean_entrezs <- sapply(entrezs, function(x) if(is.null(x)) return(NA) else return(x[1]))
  to_remove <- is.na(clean_entrezs) | duplicated(clean_entrezs)
  norm_counts <- norm_counts[!to_remove, ]
  rownames(norm_counts) <- clean_entrezs[!to_remove]

  return(norm_counts)
}

get_cannonical_tfs <- function(pathways_file, xref_file){

  load(pathways_file)
  load(xref_file)

  grn <- dorothea::entire_database
  all_tfs <- unique(grn$tf)

  last_nodes <- unlist(lapply(pathways$pathigraphs, function(x) lapply(x$effector.subgraphs, function(y) attributes(y)$last_node)), recursive = F)
  names(last_nodes) <- sapply(strsplit(names(last_nodes), "\\."), "[[", 2)

  entrezs <- unlist(lapply(pathways$pathigraphs, function(x) lapply(x$effector.subgraphs, function(y) {
    V(y)$genesList[[which(V(y)$name==attributes(y)$last_node)]]
  })), recursive = F)
  names(entrezs) <- sapply(strsplit(names(entrezs), "\\."), "[[", 2)
#
  sig_tfs <- unlist(lapply(pathways$pathigraphs, function(x) lapply(x$effector.subgraphs, function(y) {
    entrezs <- setdiff(V(y)$genesList[[which(V(y)$name==attributes(y)$last_node)]], c("/", "NA", NA))
    hgncs <- unlist(xref$from_entrez_to_hgnc[entrezs])
    hgncs[hgncs %in% all_tfs]
  })), recursive = F)
  names(sig_tfs) <- sapply(strsplit(names(sig_tfs), "\\."), "[[", 2)

  ## remove non included signalling cascades
  # last_nodes <- last_nodes[names(last_nodes) %in% rownames(sig_mat)]
  # entrezs <- entrezs[names(entrezs) %in% rownames(sig_mat)]
  # sig_tfs <- sig_tfs[names(sig_tfs) %in% rownames(sig_mat)]

  return(sig_tfs)

}








