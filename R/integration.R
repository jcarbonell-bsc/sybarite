
correlate_sbes <- function(sbe1, sbe2=NULL, cor_method="pearson", cor_thresh=0.5, adj_pvalue_thresh=0.1){

  require(reshape2)

  if(is.null(sbe2)){
    sbe2 <- sbe1
  }

  mat1 <- active(sbe1)
  mat2 <- active(sbe2)
  mc <- cor(t(mat1), t(active(sbe2)))

  # result matrices
  mc <- matrix(0, nrow=nrow(mat1), ncol=nrow(mat2), dimnames = list(rownames(mat1), rownames(mat2)))
  mc_pvalue <- matrix(1, nrow=nrow(mat1), ncol=nrow(mat2), dimnames = list(rownames(mat1), rownames(mat2)))
  for(i in 1:nrow(mat1)){
    for(j in 1:nrow(mat2)){
      ct <- cor.test(mat1[i,], mat2[j,], method = cor_method)
      mc[i,j] <- ct$estimate
      mc_pvalue[i,j] <- ct$p.value
    }
  }
  mc[is.na(mc)] <- 0
  mc_pvalue[is.na(mc_pvalue)] <- 1
  # TO-FIX
  # mc_adj_pvalue <- matrix(p.adjust(c(mc_pvalue), method = "fdr"), nrow=nrow(mc_pvalue), ncol=ncol(mc_pvalue), dimnames=dimnames(mc_pvalue))

  # summary
  mc_summ <- melt(mc)
  mc_pvalue_summ <- melt(mc_pvalue)
  # mc_adj_pvalue_summ <- melt(mc_adj_pvalue)

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

  key <- as.character(t(apply(summ, 1, function(x) paste0(sort(x[c("sbe1","sbe2")]), collapse="__"))))
  summ <- summ[!duplicated(key), ]

  summ$adj.pvalue <- p.adjust(summ$pvalue, method="fdr")

  sig <- summ[abs(summ$cor)>=cor_thresh & summ$adj.pvalue<=adj_pvalue_thresh, ]
  sig <- sig[sig$sbe1!=sig$sbe2,]


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








