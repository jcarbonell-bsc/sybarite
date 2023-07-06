
rm(list=ls())

library(Matrix)
library(Seurat)
library(pheatmap)
library(hipathia)
library(reshape2)
library(viper)
library(AUCell)
library(igraph)
library(hdf5r)

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

# install.packages("hdf5r", configure.args="--with-hdf5=/usr/bin/h5cc")
#
# ## Crohn disease
# patients
psobj <- load_me("../data/crohn/donor_T202.RData")
meta <- data.frame(
  batch = rep("b1", ncol(psobj)),
  stringsAsFactors = F,
  row.names = colnames(psobj)
)
psobj <- load_me("../data/crohn/donor_T189.RData")
sobj <- subset(psobj, cells = sample(1:ncol(psobj), 1000))
load("annots/xref.RData")
tmat <- translate_gene_ids(sobj@assays$RNA, xref$from_ensembl_to_hgnc)
rownames(tmat) <- as.character(rownames(tmat))
meta <- meta[colnames(tmat),, drop=F]
sobj <- init_seurat_object_from_sbe(tmat, meta, nfeatures = nrow(tmat), label_cell_type = T)

# # ovarian cancer
# load("../data/ovarian_cancer/spectrum_ov_007__epithelial.RData")
# sobj <- subset(spectrum_ov_007, cells = sample(1:ncol(spectrum_ov_007), 200))
# load("annots/xref.RData")
# tmat <- translate_gene_ids(sobj@assays$RNA, xref$from_ensembl_to_hgnc)
# sobj <- init_seurat_object_from_sbe(tmat, sobj@meta.data)

# mat <- readRDS("/home/jose/projects/data/single_cell_datasets/liver/21875532-af7b-4978-86c3-ca24f681b178/GSE124395_Normalhumanliverdata.RData")
# raw_meta <- read.table("/home/jose/projects/data/single_cell_datasets/liver/LiverCellAtlasHeterogeneity 2023-03-09 16.28.tsv", header=T, sep="\t", comment.char = "")
# sub_mat <- mat[, grep("P310_", colnames(mat))]
# meta <- data.frame(
#   batch = rep("b1", ncol(sub_mat)),
#   stringsAsFactors = F,
#   row.names = colnames(sub_mat)
# )
# sobj <- init_seurat_object_from_sbe(sub_mat, meta, label_cell_type = T)
# sobj <- subset(sobj, cells=sample(1:ncol(sobj), 1000))

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

# COVID19
# c141 <- Seurat::Read10X_h5("../../data/single_cell_datasets/GSE145926_covid19/GSM4339769_C141_filtered_feature_bc_matrix.h5")
# c141_meta <- data.frame(sample=rep("C141", ncol(c141)), group=rep("mild", ncol(c141)), row.names=colnames(c141), stringsAsFactors = F)
# c141_sel <- sample(1:nrow(c141_meta), 1000)
# c143 <- Seurat::Read10X_h5("../../data/single_cell_datasets/GSE145926_covid19/GSM4339771_C143_filtered_feature_bc_matrix.h5")
# c143_meta <- data.frame(sample=rep("C143", ncol(c143)), group=rep("severe", ncol(c143)), row.names=colnames(c143), stringsAsFactors = F)
# c143_sel <- sample(1:nrow(c143_meta), 1000)
# cmat <- cbind(c141[, c141_sel], c143[, c143_sel])
# meta <- rbind(c141_meta[c141_sel, ], c143_meta[c143_sel, ])
# sobj <- init_seurat_object_from_sbe(log(1 + cmat), meta, batch = meta$sample, label_cell_type = T, nfeatures = 2000)
# sobj@meta.data$sample_cell_type <- paste0(sobj@meta.data$sample, "__", sobj@meta.data$cell_type)
#
# # negative values from batch effect correction
# mat <- active(sobj)
# hist(mat, 100)
# th <- quantile(mat, .05)
# mat[mat<th] <- th
# hist(mat, 100)
# mat <- mat - th
# hist(mat, 100)
# sobj@assays$RNA <- Matrix(mat)




DimPlot(object = sobj, reduction = "umap", group.by = "seurat_clusters")
DimPlot(object = sobj, reduction = "umap", group.by = "sample")
DimPlot(object = sobj, reduction = "umap", group.by = "cell_type")
DimPlot(object = sobj, reduction = "umap", group.by = "Phase")
pheatmap(table(sobj@meta.data$sample, sobj@meta.data$cell_type), cluster_cols = F, cluster_rows = F)
pheatmap(log(1+table(sobj@meta.data$sample, sobj@meta.data$cell_type)), cluster_cols = F, cluster_rows = F)


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

# m_obj <- run_metabolic_tasks_inference(
#   sobj = sobj,
#   tasks_file = "annots/single_task_ko_Recon2_2_1.tsv",
#   info_file = "annots/richelle__simplify_tasks__clean.csv"
# )
#
# DimPlot(object = m_obj, reduction = "umap", group.by = "seurat_clusters")
# DimPlot(object = m_obj, reduction = "umap", group.by = "seurat_clusters_rna")
# DimPlot(object = m_obj, reduction = "umap", group.by = "cell_type")

scfea <- run_FBA(sobj, scfea_path = "~/projects/sc_systems/sybarite_dev/third_party/scFEA/")

DimPlot(object = scfea$flux, reduction = "umap", group.by = "seurat_clusters")
DimPlot(object = scfea$flux, reduction = "umap", group.by = "seurat_clusters_rna")
DimPlot(object = scfea$flux, reduction = "umap", group.by = "sample")
DimPlot(object = scfea$flux, reduction = "umap_rna", group.by = "seurat_clusters", split.by = "sample")

DimPlot(object = scfea$balance, reduction = "umap", group.by = "seurat_clusters")
DimPlot(object = scfea$balance, reduction = "umap", group.by = "seurat_clusters_rna")
DimPlot(object = scfea$balance, reduction = "umap_rna", group.by = "seurat_clusters")

# flux <- active(scfea$flux)
# for(i in 1:nrow(flux)){
#   boxplot(flux[i,] ~ sobj@meta.data$Phase, main=rownames(flux)[i], outline=F)
#   scan()
# }
# RidgePlot(scfea$flux, features = c("M-2", "M-4"), ncol = 2, group.by = "Phase")
#
# balance <- active(scfea$balance)
# for(i in 1:nrow(balance)){
#   boxplot(balance[i,] ~ sobj@meta.data$Phase, main=rownames(balance)[i], outline=F)
#   scan()
# }
# RidgePlot(scfea$balance, features = c("PRPP", "Pyruvate", "Cholesterol", "Glutamine"), ncol = 2, group.by = "Phase")

# scfea_cc <- run_FBA(sobj_cc, scfea_path = "~/projects/sc_systems/sybarite_dev/third_party/scFEA/")
#
# DimPlot(object = scfea_cc$flux, reduction = "umap", group.by = "seurat_clusters")
# DimPlot(object = scfea_cc$flux, reduction = "umap", group.by = "seurat_clusters_rna")
# DimPlot(object = scfea$flux, reduction = "umap_rna", group.by = "seurat_clusters")


################################################################################
##### SIGNALLING
################################################################################

#TO-FIX
# remove disease pathways (see /home/jose/projects/smoothie/annots/kegg/hipathia_pathways.tsv)
# h_obj <- run_signalling_inference(sobj, hipathia_pathways_file = "annots/cannonical_pathways.RData", xref_file = "annots/xref.RData")
h_obj <- run_signalling_inference(sobj, hipathia_pathways_file = "annots/cannonical_pathways__non_disease.RData", xref_file = "annots/xref.RData")

DimPlot(object = h_obj, reduction = "umap", group.by = "seurat_clusters")
DimPlot(object = h_obj, reduction = "umap", group.by = "seurat_clusters_rna")
DimPlot(object = h_obj, reduction = "umap", group.by = "cell_type")
DimPlot(object = h_obj, reduction = "umap", group.by = "sample")
DimPlot(object = h_obj, reduction = "umap_rna", group.by = "seurat_clusters")


# h_tfs <- get_cannonical_tfs("annots/cannonical_pathways.RData", "annots/xref.RData")
h_tfs <- get_cannonical_tfs("annots/cannonical_pathways__non_disease.RData", "annots/xref.RData")

# hobj <- run_signalling_inference(sobj, hipathia_pathways_file = "annots/pathways.RData", xref_file = "annots/xref.RData")
# DimPlot(object = hobj, reduction = "umap", group.by = "seurat_clusters_rna")

load("annots/cannonical_pathways__non_disease.RData")
path_labels <- paste0(pathways$all.labelids[rownames(h_obj), "path.name"], "__", sapply(h_tfs, function(x) paste0(x, collapse="_"))[rownames(h_obj)])
names(path_labels) <- rownames(h_obj)


DimHeatmap(h_obj, dims = 2, cells = 200, balanced = T)


################################################################################
##### GRNs
################################################################################


tf_obj <- run_tf_activity_inference(sobj)

DimPlot(object = tf_obj, reduction = "umap", group.by = "seurat_clusters")
DimPlot(object = tf_obj, reduction = "umap", group.by = "seurat_clusters_rna")
DimPlot(object = tf_obj, reduction = "umap", group.by = "cell_type")
DimPlot(object = tf_obj, reduction = "umap", group.by = "sample")
DimPlot(object = tf_obj, reduction = "umap_rna", group.by = "seurat_clusters")
DimPlot(object = h_obj, reduction = "umap_rna", group.by = "seurat_clusters")


# save(sobj, h_obj, tf_obj, scfea, file="/media/jose/secondary/projects/sybarite/covid19.RData")
# load("/media/jose/secondary/projects/sybarite/covid19.RData")



#
# c141_cells <- colnames(sobj)[sobj@meta.data$sample=="C143"]
# c141_cells <- c141_cells[c141_cells %in% colnames(sobj)]
# c141_sobj <- subset(sobj, cells=c141_cells)
# c141_h_obj <- subset(h_obj, cells=c141_cells)
# c141_tf_obj <- subset(tf_obj, cells=c141_cells)
# c141_mflux_obj <- subset(scfea$flux, cells=c141_cells)
#
# DimPlot(object = sobj, reduction = "umap", group.by = "seurat_clusters")
# DimPlot(object = c141_h_obj, reduction = "umap_rna", group.by = "seurat_clusters_rna")
# DimPlot(object = c141_h_obj, reduction = "umap_rna", group.by = "seurat_clusters")
# DimPlot(object = c141_tf_obj, reduction = "umap_rna", group.by = "seurat_clusters")
# DimPlot(object = c141_mflux_obj, reduction = "umap_rna", group.by = "seurat_clusters")
#
#
# DimPlot(object = sobj, reduction = "umap", group.by = "seurat_clusters", split.by = "sample")
# DimPlot(object = h_obj, reduction = "umap_rna", group.by = "seurat_clusters", split.by = "sample")
# DimPlot(object = tf_obj, reduction = "umap_rna", group.by = "seurat_clusters", split.by = "sample")
# DimPlot(object = scfea$flux, reduction = "umap_rna", group.by = "seurat_clusters", split.by = "sample")
#
# p1 <- DimPlot(object = sobj, reduction = "umap", group.by = "cell_type", cells = c141_cells) + NoLegend() + ggtitle("RNA")
# p2 <- DimPlot(object = h_obj, reduction = "umap", group.by = "cell_type", cells = c141_cells) + NoLegend() + ggtitle("SIG")
# p3 <- DimPlot(object = tf_obj, reduction = "umap", group.by = "cell_type", cells = c141_cells) + NoLegend() + ggtitle("GRN")
# p4 <- DimPlot(object = scfea$flux, reduction = "umap", group.by = "cell_type", cells = c141_cells) + NoLegend() + ggtitle("FLUX")
# grid.arrange(p1, p2, p3, p4)
#
#
# p1 <- DimPlot(object = sobj, reduction = "umap", group.by = "seurat_clusters", cells = c141_cells) + NoLegend() + ggtitle("RNA")
# p2 <- DimPlot(object = h_obj, reduction = "umap_rna", group.by = "seurat_clusters", cells = c141_cells) + NoLegend() + ggtitle("SIG")
# p3 <- DimPlot(object = tf_obj, reduction = "umap_rna", group.by = "seurat_clusters", cells = c141_cells) + NoLegend() + ggtitle("GRN")
# p4 <- DimPlot(object = scfea$flux, reduction = "umap_rna", group.by = "seurat_clusters", cells = c141_cells) + NoLegend() + ggtitle("FLUX")
# grid.arrange(p1, p2, p3, p4)
#
# p1 <- DimPlot(object = sobj, reduction = "umap", group.by = "seurat_clusters", cells = c141_cells) + NoLegend() + ggtitle("RNA")
# p2 <- DimPlot(object = h_obj, reduction = "umap", group.by = "seurat_clusters_rna", cells = c141_cells) + NoLegend() + ggtitle("SIG")
# p3 <- DimPlot(object = tf_obj, reduction = "umap", group.by = "seurat_clusters_rna", cells = c141_cells) + NoLegend() + ggtitle("GRN")
# p4 <- DimPlot(object = scfea$flux, reduction = "umap", group.by = "seurat_clusters_rna", cells = c141_cells) + NoLegend() + ggtitle("FLUX")
# grid.arrange(p1, p2, p3, p4)


################################################################################
##### DIFFERENTIAL ANALYSIS
################################################################################

# field <- "sample_cell_type"
# g1 <- "C141__Macrophage"
# g2 <- "C143__Macrophage"
field <- "cell_type"
g1 <- "B_cell"
g2 <- "Epithelial_cells"
scfea$flux@meta.data$sample_cell_type <- paste0(scfea$flux@meta.data$sample, "__", scfea$flux@meta.data$cell_type)
scfea$balance@meta.data$sample_cell_type <- paste0(scfea$balance@meta.data$sample, "__", scfea$balance@meta.data$cell_type)
h_obj@meta.data$sample_cell_type <- paste0(h_obj@meta.data$sample, "__", h_obj@meta.data$cell_type)
tf_obj@meta.data$sample_cell_type <- paste0(tf_obj@meta.data$sample, "__", tf_obj@meta.data$cell_type)

# field <- "cell_type"
# g1 <- "T_cells"
# g2 <- "B_cell"


cells1 <- colnames(sobj)[sobj@meta.data[,field]==g1]
cells2 <- colnames(sobj)[sobj@meta.data[,field]==g2]

scfea$flux <- SetIdent(scfea$flux, value = field)
mflux_macrophage_markers <- FindMarkers(scfea$flux, ident.1 = g1, ident.2 = g2, test.use = "wilcox", logfc.threshold = 0)
mflux_macrophage_markers_sig <- mflux_macrophage_markers[mflux_macrophage_markers$p_val_adj<0.05, ]
cat("Found", nrow(mflux_macrophage_markers_sig), "significant Mflux markers\n")

scfea$balance <- SetIdent(scfea$balance, value = field)
mbalance_macrophage_markers <- FindMarkers(scfea$balance, ident.1 = g1, ident.2 = g2, test.use = "wilcox", logfc.threshold = 0)
mbalance_macrophage_markers_sig <- mbalance_macrophage_markers[mbalance_macrophage_markers$p_val_adj<0.05, ]
cat("Found", nrow(mbalance_macrophage_markers_sig), "significant Mbalance markers\n")


h_obj <- SetIdent(h_obj, value = field)
h_macrophage_markers <- FindMarkers(h_obj, ident.1 = g1, ident.2 = g2, test.use = "wilcox", logfc.threshold = 0)
h_macrophage_markers_sig <- h_macrophage_markers[h_macrophage_markers$p_val_adj<0.05, ]
h_macrophage_markers_sig$label <- path_labels[rownames(h_macrophage_markers_sig)]
cat("Found", nrow(h_macrophage_markers_sig), "significant SIG markers\n")

tf_obj <- SetIdent(tf_obj, value = field)
tf_macrophage_markers <- FindMarkers(tf_obj, ident.1 = g1, ident.2 = g2, test.use = "t", logfc.threshold = 0)
tf_macrophage_markers_sig <- tf_macrophage_markers[tf_macrophage_markers$p_val_adj<0.05, ]
cat("Found", nrow(tf_macrophage_markers_sig), "significant tf markers\n")



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

# BY GROUP
groups <- h_obj@meta.data$seurat_clusters_rna
tf_auto_int_by_cluster <- correlate_sbes_by_cluster(tf_obj, groups = groups, cor_method = "pearson", min_number_of_cells = 10)
h_auto_int_by_cluster <- correlate_sbes_by_cluster(h_obj, groups = groups, cor_method = "pearson", min_number_of_cells = 10)
mflux_auto_int_by_cluster <- correlate_sbes_by_cluster(scfea$flux, groups = groups, cor_method = "pearson", min_number_of_cells = 10)
mbalance_auto_int_by_cluster <- correlate_sbes_by_cluster(scfea$balance, groups = groups, cor_method = "pearson", min_number_of_cells = 10)

h_tf_int_by_cluster <- correlate_sbes_by_cluster(h_obj, tf_obj, groups, cor_method = "pearson", min_number_of_cells = 10)
h_mflux_int_by_cluster <- correlate_sbes_by_cluster(h_obj, scfea$flux, groups, cor_method = "pearson", min_number_of_cells = 10)
h_mbalance_int_by_cluster <- correlate_sbes_by_cluster(h_obj, scfea$balance, groups, cor_method = "pearson", min_number_of_cells = 10)
tf_mflux_int_by_cluster <- correlate_sbes_by_cluster(tf_obj, scfea$flux, groups, cor_method = "pearson", min_number_of_cells = 10)
tf_mbalance_int_by_cluster <- correlate_sbes_by_cluster(tf_obj, scfea$balance, groups, cor_method = "pearson", min_number_of_cells = 10)
mm_int_by_cluster <- correlate_sbes_by_cluster(scfea$flux, scfea$balance, groups, cor_method = "pearson", min_number_of_cells = 10)

#
# for(i in 1:nrow(h_tf_int_by_cluster$global$sig)){
#   sbe1_id <- h_tf_int_by_cluster$global$sig$sbe1[i]
#   sbe2_id <- h_tf_int_by_cluster$global$sig$sbe2[i]
#   a <- do.call("rbind",lapply(h_tf_int_by_cluster, function(x) x$summ[x$summ$sbe1==sbe1_id & x$summ$sbe2==sbe2_id,]))
#   par(oma=c(1,10,1,1))
#   barplot(a$cor, names.arg = rownames(a), horiz=T, las=2, main=paste0(sbe1_id, " vs ", sbe2_id))
#   scan()
# }
#
# # sbe1_id <- "N-hsa05162-27 70"
# # sbe2_id <- "ADNP2"
# groups <- h_obj@meta.data$seurat_clusters_rna
# cell1 <- "rna_02"
# cell2 <- "rna_12"
# plot(active(h_obj)[sbe1_id, groups==cell1], active(tf_obj)[sbe2_id, groups==cell1])
# plot(active(h_obj)[sbe1_id, groups==cell2], active(tf_obj)[sbe2_id, groups==cell2])
#
# head(h_tf_int$sig)
# head(h_tf_int$sig[order(h_tf_int$sig$cor, decreasing = T), ])
#
# plot(active(h_obj)["N-hsa04380-48", ], active(tf_obj)["NME2", ])
# p1 <- Seurat::FeaturePlot(h_obj, "N-hsa04380-48", reduction="umap_rna", min.cutoff = "q25", max.cutoff = "q75")
# p2 <- Seurat::FeaturePlot(tf_obj, "NME2", reduction="umap_rna", min.cutoff = "q25", max.cutoff = "q75")
#
# plot_pair("N-hsa04380-48", "NME2", h_obj, tf_obj)
#
# sel <- h_tf_int$sig[sample(1:nrow(h_tf_int$sig), 1),]
# plot_pair(sel$sbe1, sel$sbe2, h_obj, tf_obj)



plot_random_sig_pair(h_tf_int, h_obj, tf_obj)

plot_random_sig_pair(h_mflux_int, h_obj, scfea$flux)

plot_random_sig_pair(tf_mbalance_int, tf_obj, scfea$balance)

plot_random_sig_pair(h_mbalance_int, h_obj, scfea$balance)
plot_random_sig_pair(h_mflux_int, h_obj, scfea$flux)




### SBE multilayer


layers <- rbind(
  cbind(rownames(tf_obj), "TF"),
  cbind(rownames(h_obj), "SIG"),
  cbind(rownames(scfea$flux), "Mflux"),
  cbind(rownames(scfea$balance), "Mbalance")
)
colnames(layers) <- c("sbe", "layer")

lcolors <- RColorBrewer::brewer.pal(n = 4, "Set1")
names(lcolors) <- c("TF", "SIG", "Mflux", "Mbalance")


#
# inter_cor <- .5
# intra_cor <- .75
#
# tf_final <- select_interactions(tf_auto_int, 0.01, intra_cor)
# h_final <- select_interactions(h_auto_int, 0.01, intra_cor)
# mflux_final <- select_interactions(mflux_auto_int, 0.01, intra_cor)
# mbalance_final <- select_interactions(mbalance_auto_int, 0.01, intra_cor)
#
# h_tf_final <- select_interactions(h_tf_int, 0.01, inter_cor)
# h_mflux_final <- select_interactions(h_mflux_int, 0.01, inter_cor)
# h_mbalance_final <- select_interactions(h_mbalance_int, 0.01, inter_cor)
# tf_mflux_final <- select_interactions(tf_mflux_int, 0.01, inter_cor)
# tf_mbalance_final <- select_interactions(tf_mbalance_int, 0.01, inter_cor)
# mm_final <- select_interactions(mm_int, 0.01, inter_cor)
#
# all_final <- rbind(
#   # intra
#   cbind(tf_final, context="TF"),
#   cbind(h_final, context="SIG"),
#   cbind(mbalance_final, context="Mflux"),
#   # inter
#   cbind(h_tf_final, context="SIG_TF"),
#   cbind(h_mflux_final, context="SIG_Mflux"),
#   cbind(tf_mflux_final, context="TF_Mflux")
# )
# sort(table(all_final$context), decreasing = T)
#
#
# mln <- graph.data.frame(all_final, directed = F)
# V(mln)$color <- lcolors[layers[match(V(mln)$name, layers[, 1]), 2]]
# V(mln)$size <- 5
# plot.igraph(mln, vertex.label=NA)
# legend("bottomright", legend=names(lcolors), col=lcolors, lwd=2, bty="n")
#
# ### FIND COMMUNITIES
#
# comm <- cluster_fast_greedy(mln)
# comm <- cluster_infomap(mln, nb.trials = 100)
# groups <- communities(comm)
# cat("Found", length(groups), "communities\n")
# hist(sizes(comm), 50)
#
# plot.igraph(mln, vertex.label=NA, mark.groups = comm)
# legend("bottomright", legend=names(lcolors), col=lcolors, lwd=2, bty="n")





######

# ints <- list(
#   SIG=h_auto_int, TF=tf_auto_int, Mflux=mflux_auto_int,
#   SIG_TF=h_tf_int, SIG_Mflux=h_mflux_int, TF_Mflux=tf_mflux_int
# )
ints <- list(
  SIG=h_auto_int, TF=tf_auto_int, Mbalance=mbalance_auto_int,
  SIG_TF=h_tf_int, SIG_Mbalance=h_mbalance_int, TF_Mbalance=tf_mbalance_int
)
cors <- c(0.85, 0.95, 0.5,  0.25, 0.25, 0.25)
cors <- c(0.85, 0.85, 0.85,  0.5, 0.5, 0.5)
mln <- create_multilayer(ints, cors, layers = layers, lcolors = lcolors)

# comm <- find_communities(mln, method=cluster_fast_greedy, lcolors = lcolors)
# comm <- find_communities(mln, method=cluster_edge_betweenness, lcolors = lcolors)
comm <- find_communities(mln, method=cluster_infomap, lcolors = lcolors, nb.trials=1, modularity = 1)
# comm <- find_communities(mln, method=cluster_fluid_communities, lcolors = lcolors, no.of.communities = 10)
comm$subg_stats

for(i in 1:length(comm$subg)){
  plot.igraph(comm$subg[[i]], main=i)
  legend("bottomright", legend=names(lcolors), col=lcolors, lwd=2, bty="n")
  scan()
}

plot(comm$subg[[26]])

objs <- list(SIG=h_obj, TF=tf_obj, Mflux=scfea$flux, Mbalance=scfea$balance)

plot_elements(comm$subg[[7]], objs=objs, layers = layers)
plot_elements(comm$subg[[28]], objs=objs, layers = layers)


a <- comm$subg[[28]]
b <- igraph::delete.vertices(a, c("ETV7", "SNAI3"))
plot.igraph(b, vertex.size=10, vertex.label=paste("\n\n", V(b)$name))
plot_elements(b, objs=objs, layers = layers)

###### COMPARISON




################################################################################
##### MECHANISTIC INTEGRATION
################################################################################











