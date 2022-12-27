
rm(list=ls())

library(Matrix)
library(Seurat)
library(pheatmap)
library(hipathia)
library(reshape2)
library(viper)
library(AUCell)

source("R/metabol.R")
source("R/signalling.R")
source("R/hipathia_wrapper.R")
source("R/grns.R")
source("R/enrichment.R")
source("R/utils.R")
source("R/plots.R")

load("annots/xref.RData")

# load("../GSM4339771_C143_filtered_feature_bc_matrix__primary.RData")

# load("~/projects/analysis/cll_davide/CLL5__seurat_object.RData")
# pbmc <- subset(pbmc, cells = sample(1:ncol(pbmc), 1000))

# load("../data/ovarian_cancer/spectrum_ov_007__epithelial.RData")
# pbmc <- subset(spectrum_ov_007, cells = sample(1:ncol(spectrum_ov_007), 1000))

load("../data/crohn/donor_T202.RData")
pbmc <- subset(sobj, cells = sample(1:ncol(sobj), 1000))
tmat <- translate_ids(pbmc@assays$RNA, xref$from_ensembl_to_hgnc)
pbmc <- init_seurat_object_from_sbe(tmat, pbmc@meta.data)



#
DimPlot(object = pbmc, reduction = "umap", group.by = "seurat_clusters")
DimPlot(object = pbmc, reduction = "umap", group.by = "cell_type")


tobj <- run_metabolic_tasks_inference(
  sobj = pbmc,
  tasks_file = "../annotations/single_task_ko_Recon2_2_1.tsv",
  info_file = "../annotations/richelle__simplify_tasks__clean.csv"
)

DimPlot(object = tobj, reduction = "umap", group.by = "seurat_clusters")
DimPlot(object = tobj, reduction = "umap", group.by = "seurat_clusters_rna")
DimPlot(object = tobj, reduction = "umap", group.by = "cell_type")
DimPlot(object = tobj, reduction = "umap_rna", group.by = "seurat_clusters_rna")
DimPlot(object = tobj, reduction = "umap_rna", group.by = "seurat_clusters")

FeaturePlot(tobj, features = c("task9"), reduction="umap_rna")


mmat <- as.array(tobj@assays$RNA[])


# norm counts
norm_counts <- Matrix::as.array(pbmc@assays$RNA[,1:ncol(pbmc)])
rownames(norm_counts) <- as.character(sapply(strsplit(rownames(norm_counts),"\\."),"[[",1))
auc_res <- AUCell_run(norm_counts, geneSets = mtasks$gene2task)
auc_mat <- getAUC(auc_res)
auc_obj <- init_seurat_object_from_sbe(auc_mat, pbmc@meta.data)

DimPlot(object = auc_obj, reduction = "umap", group.by = "seurat_clusters")
DimPlot(object = auc_obj, reduction = "umap", group.by = "seurat_clusters_rna")
DimPlot(object = auc_obj, reduction = "umap", group.by = "cell_type")
DimPlot(object = auc_obj, reduction = "umap_rna", group.by = "seurat_clusters")


mtasks <- load_metabolic_tasks(tasks_file = "../annotations/single_task_ko_Recon2_2_1.tsv", info_file = "../annotations/richelle__simplify_tasks__clean.csv")
e_obj <- run_enrichment_inference(pbmc, mtasks$gene2task, method="escape")
mmat2 <- as.array(e_obj@assays$RNA[])

DimPlot(object = e_obj, reduction = "umap", group.by = "seurat_clusters")
DimPlot(object = e_obj, reduction = "umap", group.by = "seurat_clusters_rna")
DimPlot(object = e_obj, reduction = "umap", group.by = "cell_type")
DimPlot(object = e_obj, reduction = "umap_rna", group.by = "seurat_clusters")

## HIPATHIA

hobj <- run_signalling_inference(pbmc, hipathia_pathways_file = "annots/cannonical_pathways.RData", xref_file = "annots/xref.RData")
# hobj <- run_signalling_inference(pbmc, hipathia_pathways_file = "annots/pathways.RData", xref_file = "annots/xref.RData")

DimPlot(object = hobj, reduction = "umap", group.by = "seurat_clusters")
DimPlot(object = hobj, reduction = "umap", group.by = "seurat_clusters_rna")
DimPlot(object = hobj, reduction = "umap", group.by = "cell_type")
DimPlot(object = hobj, reduction = "umap_rna", group.by = "seurat_clusters")

FeaturePlot(hobj, features = c("P-hsa03320-37"), reduction = "prueba")

# smat <- as.array(hobj@assays$RNA[])
sig_mat <- as.array(hobj@assays$RNA[])

i <- sample(1:nrow(path_vals),1); hist(smat[i, ], 100)


Idents(hobj) <- "seurat_clusters_rna"
rna_cluster_markers <- FindAllMarkers(hobj, only.pos = TRUE, min.pct = 0, logfc.threshold = 0, test.use = "wilcox")


### VIPER

tf_obj <- run_tf_activity_inference(pbmc)

DimPlot(object = tf_obj, reduction = "umap", group.by = "seurat_clusters")
DimPlot(object = tf_obj, reduction = "umap", group.by = "seurat_clusters_rna")
DimPlot(object = tf_obj, reduction = "umap", group.by = "cell_type")
DimPlot(object = tf_obj, reduction = "umap_rna", group.by = "seurat_clusters")

tf_mat <- as.array(tf_obj@assays$RNA[])


#######


ms_cor <- cor(t(tf_mat), t(sig_mat), method="pearson")
ms_cor[is.na(ms_cor)] <- 0

hist(ms_cor, 100)
pheatmap(ms_cor)

ms_cor_df <- melt(ms_cor, as.is = T)

sig_ms_cor <- ms_cor_df[abs(ms_cor_df$value)>.5, ]
print(dim(sig_ms_cor))

for(i in 1:nrow(sig_ms_cor)){
  plot(tf_mat[sig_ms_cor$Var1[i], ], sig_mat[sig_ms_cor$Var2[i], ], xlab=sig_ms_cor$Var1[i], ylab=sig_ms_cor$Var2[i])
  scan()
}


plot(mmat["task137", ], smat["P-hsa04066-44", ])




load("annots/cannonical_pathways.RData")
load("annots/xref.RData")
grn <- dorothea::entire_database
all_tfs <- unique(grn$tf)
#
last_nodes <- unlist(lapply(pathways$pathigraphs, function(x) lapply(x$effector.subgraphs, function(y) attributes(y)$last_node)), recursive = F)
names(last_nodes) <- sapply(strsplit(names(last_nodes), "\\."), "[[", 2)
#
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
last_nodes <- last_nodes[names(last_nodes) %in% rownames(sig_mat)]
entrezs <- entrezs[names(entrezs) %in% rownames(sig_mat)]
sig_tfs <- sig_tfs[names(sig_tfs) %in% rownames(sig_mat)]
##

u_sig_tfs <- unique(unlist(sig_tfs))

viper_sig_tf_indexes <- lapply(rownames(tf_mat), function(x){
  which(sapply(sig_tfs, function(stfs) x %in% stfs))
})
viper_sig_tf_indexes_valid <- viper_sig_tf_indexes[sapply(viper_sig_tf_indexes, length)>0]


singles <- u_sig_tfs[!(u_sig_tfs %in% rownames(tf_mat))]
singles


sig_to_viper_indexes <- lapply(sig_tfs, function(x) match(x, rownames(tf_mat)))
sig_to_viper_indexes_best <- sapply(sig_to_viper_indexes, "[[", 1)

sig_viper_cross <- data.frame(
  id1 = names(sig_to_viper_indexes),
  index1 = match(names(sig_to_viper_indexes), rownames(sig_mat)),
  id2 = rownames(tf_mat)[sig_to_viper_indexes_best],
  index2 = sig_to_viper_indexes_best,
  stringsAsFactors = F
)
sig_viper_cross <- sig_viper_cross[!is.na(sig_viper_cross$id2),]

for(i in 1:nrow(sig_viper_cross)){
  print(sig_viper_cross[i, ])
  plot(sig_mat[sig_viper_cross$index1[i],], tf_mat[sig_viper_cross$index2[i],], main=paste0(sig_viper_cross$id1[i], "-", sig_viper_cross$id2[i]))
  scan()
}

all_paths <- unlist(lapply(pathways$pathigraphs, "[[", "effector.subgraphs"), recursive=F, use.names = T)
names(all_paths) <- sapply(strsplit(names(all_paths), "\\."), "[[", 2)


plot(all_paths[["N-hsa04010-68"]])

plot(sig_mat["N-hsa04010-68", ], tf_mat["ATF2", ])
