
rm(list=ls())

library(Matrix)
library(Seurat)
library(pheatmap)
library(hipathia)
library(reshape2)
library(viper)
library(AUCell)
library(Igraph)

source("R/metabol.R")
source("R/scFEA_wrapper.R")
source("R/signalling.R")
source("R/hipathia_wrapper.R")
source("R/grns.R")
source("R/enrichment.R")
source("R/utils.R")
source("R/plots.R")
source("R/integration.R")

################################################################################
##### LOAD DATA
################################################################################

# Chron
load("../data/crohn/donor_T202.RData")
sobj <- subset(sobj, cells = sample(1:ncol(sobj), 1000))
load("annots/xref.RData")
tmat <- translate_gene_ids(sobj@assays$RNA, xref$from_ensembl_to_hgnc)
sobj <- init_seurat_object_from_sbe(tmat, sobj@meta.data)

# Dentate gyrus
# load("../data/dentate_gyrus/10X43_1.loom.RData")
# spliced_mat <- Matrix::as.array(ldat$spliced)
# rownames(spliced_mat) <- toupper(rownames(spliced_mat))
# meta <- data.frame(
#   batch = rep("b1", ncol(spliced_mat)),
#   stringsAsFactors = F,
#   row.names = colnames(spliced_mat)
# )
# sobj <- init_seurat_object_from_sbe(log(1 + spliced_mat), meta)

DimPlot(object = sobj, reduction = "umap", group.by = "seurat_clusters")
DimPlot(object = sobj, reduction = "umap", group.by = "cell_type")


### cell cycle - metabolism

sobj <- CellCycleScoring(sobj, s.features = cc.genes$s.genes, g2m.features = cc.genes$g2m.genes, set.ident = F)

DimPlot(object = sobj, reduction = "umap", group.by = "Phase")

sobj_cc <- ScaleData(sobj, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(sobj))

DimPlot(object = sobj_cc, reduction = "umap", group.by = "seurat_clusters")
DimPlot(object = sobj_cc, reduction = "umap", group.by = "cell_type")
DimPlot(object = sobj_cc, reduction = "umap", group.by = "Phase")


################################################################################
##### METABOLISM
################################################################################

m_obj <- run_metabolic_tasks_inference(
  sobj = sobj,
  tasks_file = "annots/single_task_ko_Recon2_2_1.tsv",
  info_file = "annots/richelle__simplify_tasks__clean.csv"
)

DimPlot(object = m_obj, reduction = "umap", group.by = "seurat_clusters")
DimPlot(object = m_obj, reduction = "umap", group.by = "seurat_clusters_rna")
DimPlot(object = m_obj, reduction = "umap", group.by = "cell_type")

scfea <- run_FBA(sobj, scfea_path = "~/projects/sc_systems/sybarite_dev/third_party/scFEA/")

DimPlot(object = scfea$flux, reduction = "umap", group.by = "seurat_clusters")
DimPlot(object = scfea$flux, reduction = "umap", group.by = "seurat_clusters_rna")

DimPlot(object = scfea$balance, reduction = "umap", group.by = "seurat_clusters")
DimPlot(object = scfea$balance, reduction = "umap", group.by = "seurat_clusters_rna")
DimPlot(object = scfea$balance, reduction = "umap_rna", group.by = "seurat_clusters")

flux <- active(scfea$flux)
for(i in 1:nrow(flux)){
  boxplot(flux[i,] ~ sobj@meta.data$Phase, main=rownames(flux)[i], outline=F)
  scan()
}
RidgePlot(scfea$flux, features = c("M-2", "M-4"), ncol = 2, group.by = "Phase")

balance <- active(scfea$balance)
for(i in 1:nrow(balance)){
  boxplot(balance[i,] ~ sobj@meta.data$Phase, main=rownames(balance)[i], outline=F)
  scan()
}
RidgePlot(scfea$balance, features = c("PRPP", "Pyruvate", "Cholesterol", "Glutamine"), ncol = 2, group.by = "Phase")

scfea_cc <- run_FBA(sobj_cc, scfea_path = "~/projects/sc_systems/sybarite_dev/third_party/scFEA/")

DimPlot(object = scfea_cc$flux, reduction = "umap", group.by = "seurat_clusters")
DimPlot(object = scfea_cc$flux, reduction = "umap", group.by = "seurat_clusters_rna")

################################################################################
##### SIGNALLING
################################################################################

h_obj <- run_signalling_inference(sobj, hipathia_pathways_file = "annots/cannonical_pathways.RData", xref_file = "annots/xref.RData")

DimPlot(object = h_obj, reduction = "umap", group.by = "seurat_clusters")
DimPlot(object = h_obj, reduction = "umap", group.by = "seurat_clusters_rna")
DimPlot(object = h_obj, reduction = "umap", group.by = "cell_type")

h_tfs <- get_cannonical_tfs("annots/cannonical_pathways.RData", "annots/xref.RData")

hobj <- run_signalling_inference(sobj, hipathia_pathways_file = "annots/pathways.RData", xref_file = "annots/xref.RData")
DimPlot(object = hobj, reduction = "umap", group.by = "seurat_clusters_rna")


################################################################################
##### GRNs
################################################################################


tf_obj <- run_tf_activity_inference(sobj)

DimPlot(object = tf_obj, reduction = "umap", group.by = "seurat_clusters")
DimPlot(object = tf_obj, reduction = "umap", group.by = "seurat_clusters_rna")
DimPlot(object = tf_obj, reduction = "umap", group.by = "cell_type")
DimPlot(object = tf_obj, reduction = "umap_rna", group.by = "seurat_clusters")




################################################################################
##### CELL-CELL COMMUNICATION
################################################################################



################################################################################
##### INTEGRATION
################################################################################

tf_auto_int <- correlate_sbes(tf_obj)
h_auto_int <- correlate_sbes(h_obj)
mflux_auto_int <- correlate_sbes(scfea$flux)
mbalance_auto_int <- correlate_sbes(scfea$balance)

h_tf_int <- correlate_sbes(h_obj, tf_obj)
h_mflux_int <- correlate_sbes(h_obj, scfea$flux)
h_mbalance_int <- correlate_sbes(h_obj, scfea$balance)
tf_mflux_int <- correlate_sbes(tf_obj, scfea$flux)
tf_mbalance_int <- correlate_sbes(tf_obj, scfea$balance)
mm_int <- correlate_sbes(scfea$flux, scfea$balance)


head(h_tf_int$sig)
head(h_tf_int$sig[order(h_tf_int$sig$cor, decreasing = T), ])

plot(active(h_obj)["N-hsa04380-48", ], active(tf_obj)["NME2", ])
p1 <- Seurat::FeaturePlot(h_obj, "N-hsa04380-48", reduction="umap_rna", min.cutoff = "q25", max.cutoff = "q75")
p2 <- Seurat::FeaturePlot(tf_obj, "NME2", reduction="umap_rna", min.cutoff = "q25", max.cutoff = "q75")

plot_pair("N-hsa04380-48", "NME2", h_obj, tf_obj)

sel <- h_tf_int$sig[sample(1:nrow(h_tf_int$sig), 1),]
plot_pair(sel$sbe1, sel$sbe2, h_obj, tf_obj)



plot_random_sig_pair(h_tf_int, h_obj, tf_obj)

plot_random_sig_pair(h_mflux_int, h_obj, scfea$flux)

plot_random_sig_pair(tf_mbalance_int, tf_obj, scfea$balance)

plot_random_sig_pair(h_mbalance_int, h_obj, scfea$balance)
plot_random_sig_pair(h_mflux_int, h_obj, scfea$flux)




### SBE multilayer


layers <- rbind(
  cbind(rownames(tf_obj), "TF"),
  cbind(rownames(h_obj), "SIG"),
  cbind(rownames(scfea_cc$flux), "Mflux"),
  cbind(rownames(scfea_cc$balance), "Mbalance")
)
colnames(layers) <- c("sbe", "layer")

lcolors <- RColorBrewer::brewer.pal(n = 4, "Set1")
names(lcolors) <- c("TF", "SIG", "Mflux", "Mbalance")

select_interactions <- function(sbe_cor, adj_pvalue=0.01, cor=0.5) {
  out <- sbe_cor$sig[ sbe_cor$sig$adj.pvalue<adj_pvalue & abs(sbe_cor$sig$cor)>cor, ]
  cat("Selected", nrow(out),"interactions\n")
  return(out)
}

inter_cor <- .5
intra_cor <- .25

tf_final <- select_interactions(tf_auto_int, 0.01, .95)#intra_cor)
h_final <- select_interactions(h_auto_int, 0.01, intra_cor)
mflux_final <- select_interactions(mflux_auto_int, 0.01, intra_cor)
mbalance_final <- select_interactions(mbalance_auto_int, 0.01, intra_cor)

h_tf_final <- select_interactions(h_tf_int, 0.01, intra_cor)
h_mflux_final <- select_interactions(h_mflux_int, 0.01, intra_cor)
h_mbalance_final <- select_interactions(h_mbalance_int, 0.01, intra_cor)
tf_mflux_final <- select_interactions(tf_mflux_int, 0.01, intra_cor)
tf_mbalance_final <- select_interactions(tf_mbalance_int, intra_cor, intra_cor)
mm_final <- select_interactions(mm_int, 0.01, intra_cor)


all_final <- rbind(
  # intra
  cbind(tf_final, context="TF"),
  cbind(h_final, context="SIG"),
  cbind(mbalance_final, context="Mflux"),
  # inter
  cbind(h_tf_final, context="SIG_TF"),
  cbind(h_mflux_final, context="SIG_Mflux"),
  cbind(tf_mflux_final, context="TF_Mflux")
)
sort(table(all_final$context), decreasing = T)


mln <- graph.data.frame(all_final, directed = F)
V(mln)$color <- lcolors[layers[match(V(mln)$name, layers[, 1]), 2]]
V(mln)$size <- 5
plot.igraph(mln, vertex.label=NA)
legend("bottomright", legend=names(lcolors), col=lcolors, lwd=2, bty="n")

### FIND COMMUNITIES

# comm <- cluster_fast_greedy(mln)
comm <- cluster_infomap(mln)

plot.igraph(mln, vertex.label=NA, mark.groups = comm)
legend("bottomright", legend=names(lcolors), col=lcolors, lwd=2, bty="n")

groups <- communities(comm)
subg <- lapply(groups, function(x) induced_subgraph(mln, vids = match(x, V(mln)$name)))
for(i in 1:length(subg)){
  plot.igraph(subg[[i]], vertex.label=NA)
  legend("bottomright", legend=names(lcolors), col=lcolors, lwd=2, bty="n")
  scan()
}


################################################################################
##### MECHANISTIC INTEGRATION
################################################################################











