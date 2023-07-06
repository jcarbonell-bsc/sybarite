
correlate_sbes <- function(sbe1, sbe2=NULL, cor_method="pearson", cor_thresh=0.5, adj_pvalue_thresh=0.1){

  require(reshape2)
  require(Hmisc)

  mat1 <- active(sbe1)
  if(is.null(sbe2)){
    rcor <- Hmisc::rcorr(t(mat1), type=cor_method)
  } else {
    mat2 <- active(sbe2)
    rcor <- Hmisc::rcorr(t(mat1), t(mat2), type=cor_method)
    rcor$r <- rcor$r[rownames(mat1), rownames(mat2)]
    rcor$P <- rcor$P[rownames(mat1), rownames(mat2)]
  }
  mc <- rcor$r
  mc_pvalue <- rcor$P

  # mat1 <- active(sbe1)
  # mat2 <- active(sbe2)
  # # result matrices
  # mc <- matrix(0, nrow=nrow(mat1), ncol=nrow(mat2), dimnames = list(rownames(mat1), rownames(mat2)))
  # mc_pvalue <- matrix(1, nrow=nrow(mat1), ncol=nrow(mat2), dimnames = list(rownames(mat1), rownames(mat2)))
  # for(i in 1:nrow(mat1)){
  #   for(j in 1:nrow(mat2)){
  #     ct <- cor.test(mat1[i,], mat2[j,], method = cor_method)
  #     mc[i,j] <- ct$estimate
  #     mc_pvalue[i,j] <- ct$p.value
  #   }
  # }
  # mc[is.na(mc)] <- 0
  # mc_pvalue[is.na(mc_pvalue)] <- 1
  # TO-FIX
  # mc_adj_pvalue <- matrix(p.adjust(c(mc_pvalue), method = "fdr"), nrow=nrow(mc_pvalue), ncol=ncol(mc_pvalue), dimnames=dimnames(mc_pvalue))


  if(is.null(sbe2)){
    na_template <- upper.tri(diag(ncol(mc)))+0
    na_template[na_template==0] <- NA
    mc2 <- mc * na_template
    mc_pvalue2 <- mc_pvalue * na_template
    mc_summ <- melt(mc2, na.rm = T)
    mc_pvalue_summ <- melt(mc_pvalue2, na.rm = T)
  } else {
    mc_summ <- melt(mc, na.rm = T)
    mc_pvalue_summ <- melt(mc_pvalue, na.rm = T)
  }

  # summary
  # mc_summ <- melt(mc)
  # mc_pvalue_summ <- melt(mc_pvalue)
  # # mc_adj_pvalue_summ <- melt(mc_adj_pvalue)

  summ <- data.frame(
    mc_summ,
    mc_pvalue_summ$value,
    # mc_adj_pvalue_summ$value,
    stringsAsFactors = F
  )
  # colnames(summ) <- c("sbe1", "sbe2", "cor", "pvalue", "adj.pvalue")
  colnames(summ) <- c("sbe1", "sbe2", "cor", "pvalue")

  summ$int <- sign(summ$cor)
  summ$sbe1 <- as.character(summ$sbe1)
  summ$sbe2 <- as.character(summ$sbe2)

  # key <- as.character(t(apply(summ, 1, function(x) paste0(sort(x[c("sbe1","sbe2")]), collapse="__"))))
  # summ <- summ[!duplicated(key), ]
  # summ <- summ[summ$sbe1!=summ$sbe2,]

  summ$adj.pvalue <- p.adjust(summ$pvalue, method="fdr")

  sig <- summ[abs(summ$cor)>=cor_thresh & summ$adj.pvalue<=adj_pvalue_thresh, ]

  return(
    list(
      cor_mat=mc,
      pvalue_mat=mc_pvalue,
      # adj_pvalue_mat=mc_adj_pvalue,
      summ=summ,
      sig=sig
    )
  )
}


correlate_sbes_by_cluster <- function(sbe1, sbe2=NULL, groups, cor_method="pearson", cor_thresh=0.5, adj_pvalue_thresh=0.1, min_number_of_cells=5){

  res <- list()

  res[["global"]] <- correlate_sbes(sbe1, sbe2, cor_method = cor_method, cor_thresh=cor_thresh, adj_pvalue_thresh=adj_pvalue_thresh)

  ug <- unique(groups)
  for(g in ug){
    print(g)
    sel_cells <- colnames(sbe1)[groups==g]
    if(length(sel_cells)>=min_number_of_cells){
      sub_sbe1 <- subset(sbe1, cells = sel_cells)
      if(!is.null(sbe2)){
        sub_sbe2 <- subset(sbe2, cells = sel_cells)
      } else {
        sub_sbe2 <- NULL
      }
      res[[g]] <- correlate_sbes(sub_sbe1, sub_sbe2, cor_method = cor_method, cor_thresh=cor_thresh, adj_pvalue_thresh=adj_pvalue_thresh)
    }
  }

  return(res)

}

plot_pair <- function(f1, f2, sbe1, sbe2, reduction="umap_rna"){
  require(gridExtra)
  require(ggplotify)

  f1 <- as.character(f1)
  f2 <- as.character(f2)
  ff <- function() {
    v1 <- active(sbe1)[f1, ]
    v2 <- active(sbe2)[f2, ]
    plot(v1, v2, xlab=f1, ylab=f2, pch=20, col="#0000FF77", cex=0.7)
    lines(lowess(v1, v2),col="black")
  }
  # ff()
  p0 <- as.grob(ff)

  p1 <- FeaturePlot(sbe1, f1, reduction=reduction, min.cutoff = "q25", max.cutoff = "q75")
  p2 <- FeaturePlot(sbe2, f2, reduction=reduction, min.cutoff = "q25", max.cutoff = "q75")

  gg <- grid.arrange(p1, p2, p0, layout_matrix=matrix(c(1, 2, 3,3), nrow=2))

}


plot_random_sig_pair <- function(int, sbe1, sbe2){
  sel <- int$sig[sample(1:nrow(int$sig), 1),]
  plot_pair(sel$sbe1, sel$sbe2, sbe1, sbe2)
}

select_interactions <- function(sbe_cor, adj_pvalue=0.01, cor=0.5) {
  # sbe_cor$summ <- sbe_cor$summ[sbe_cor$summ$sbe1!=sbe_cor$summ$sbe2,]
  out <- sbe_cor$summ[ sbe_cor$summ$adj.pvalue<adj_pvalue & abs(sbe_cor$summ$cor)>cor, ]
  # cat("Selected", nrow(out),"interactions\n")
  return(out)
}

create_multilayer <- function(ints, cors=rep(0.5, length(intra_ints)), layers, lcolors, paint=T){

  finals <- lapply(1:length(ints), function(i) select_interactions(ints[[i]], adj_pvalue=0.01, cor=cors[i]))
  names(finals) <- names(ints)

  all_final <- c()
  for(i in 1:length(finals)){
    if(nrow(finals[[i]])>0){
      all_final <- rbind(
        all_final,
        cbind(finals[[i]], context=names(finals)[i])
      )
    }
  }
  print(sort(table(all_final$context), decreasing = T))

  mln <- graph.data.frame(all_final, directed = F)
  V(mln)$color <- lcolors[layers[match(V(mln)$name, layers[, 1]), 2]]
  V(mln)$size <- 5

  if(paint==T){
    plot.igraph(mln, vertex.label=NA)
    legend("bottomright", legend=names(lcolors), col=lcolors, lwd=2, bty="n")
  }

  return(mln)

}

find_communities <- function(mln, method=cluster_fast_greedy, paint=T, lcolors=NULL, ...){

  comm <- method(mln, ...)

  groups <- communities(comm)
  cat("Found", length(groups), "communities\n")

  sizes <- sizes(comm)

  if(paint==T){
    hist(sizes(comm), 50)
    plot.igraph(mln, vertex.label=NA, mark.groups = comm)
    legend("bottomright", legend=names(lcolors), col=lcolors, pch=15, bty="n")
  }

  subg <- lapply(groups, function(x) induced_subgraph(mln, vids = match(x, V(mln)$name)))

  categories <- unique(E(mln)$context)
  subg_stats <- do.call("rbind", lapply(subg, function(x) table(E(x)$context)[categories]))
  subg_stats[is.na(subg_stats)] <- 0
  colnames(subg_stats) <- categories

  return(list(
    method=paste0(expression(method)),
    comm=comm,
    groups=groups,
    sizes=sizes,
    subg=subg,
    subg_stats=subg_stats
  ))
}


plot_elements <- function(el, objs, layers, reduction="umap_rna"){
  require(gridExtra)
  require(ggplotify)

  p <- list()
  for(i in 1:length(el)){

      item <- V(el)$name[i]
      type <- layers[match(item, layers[,1]), 2]

      p[[i]] <- FeaturePlot(objs[[type]], item, reduction=reduction, min.cutoff = "q25", max.cutoff = "q75", )

  }

  grid.arrange(grobs=p)

}






