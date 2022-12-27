
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


################################################################################
##### LOAD DATA
################################################################################

load("../data/crohn/donor_T202.RData")
sobj <- subset(sobj, cells = sample(1:ncol(sobj), 1000))
load("annots/xref.RData")
tmat <- translate_gene_ids(sobj@assays$RNA, xref$from_ensembl_to_hgnc)
sobj <- init_seurat_object_from_sbe(tmat, sobj@meta.data)

DimPlot(object = sobj, reduction = "umap", group.by = "seurat_clusters")
DimPlot(object = sobj, reduction = "umap", group.by = "cell_type")


################################################################################
##### METABOLISM
################################################################################

m_obj <- run_metabolic_tasks_inference(
  sobj = sobj,
  tasks_file = "../annotations/single_task_ko_Recon2_2_1.tsv",
  info_file = "../annotations/richelle__simplify_tasks__clean.csv"
)

DimPlot(object = m_obj, reduction = "umap", group.by = "seurat_clusters")
DimPlot(object = m_obj, reduction = "umap", group.by = "seurat_clusters_rna")
DimPlot(object = m_obj, reduction = "umap", group.by = "cell_type")


################################################################################
##### SIGNALLING
################################################################################

h_obj <- run_signalling_inference(sobj, hipathia_pathways_file = "annots/cannonical_pathways.RData", xref_file = "annots/xref.RData")
h_tfs <- get_cannonical_tfs("annots/cannonical_pathways.RData", "annots/xref.RData")

DimPlot(object = h_obj, reduction = "umap", group.by = "seurat_clusters")
DimPlot(object = h_obj, reduction = "umap", group.by = "seurat_clusters_rna")
DimPlot(object = h_obj, reduction = "umap", group.by = "cell_type")



# hobj <- run_signalling_inference(pbmc, hipathia_pathways_file = "annots/pathways.RData", xref_file = "annots/xref.RData")


################################################################################
##### GRNs
################################################################################


tf_obj <- run_tf_activity_inference(sobj)

DimPlot(object = tf_obj, reduction = "umap", group.by = "seurat_clusters")
DimPlot(object = tf_obj, reduction = "umap", group.by = "seurat_clusters_rna")
DimPlot(object = tf_obj, reduction = "umap", group.by = "cell_type")
DimPlot(object = tf_obj, reduction = "umap_rna", group.by = "seurat_clusters")


################################################################################
##### INTEGRATION
################################################################################







