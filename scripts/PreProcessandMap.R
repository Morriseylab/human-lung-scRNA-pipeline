#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)


library(Seurat)
library(tidyverse)
library(patchwork)
library(scExtras)


org<-'human'
input10x = 'AHJP457_D_CD45neg/scRNA/STARSolo/AHJP457_D_CD45negSolo.out/Gene/filtered/'


outdir='AHJP457_D_CD45neg/Seurat_refmap'

dir.create(outdir,recursive = T,showWarnings = F)
plotdir <- paste0(outdir,'/plots')
dir.create(plotdir,showWarnings = F)
qcdir <-paste0(outdir,'/qc')
dir.create(qcdir,showWarnings = F)


query= RunQC(dir=outdir,org=org,name="AHJP457_D_CD45",files=input10x ,filter=T, doubletdetection = T,UpperMitoCutoff=10)


reference <- readRDS("~/dsdata/projects/Morrisey/Maria/lungMAP/Periph_5sample_Intergrate/Seuratv4_Anno.RDS")


query <- SCTransform(
  object = query,
  assay = "RNA",
  vars.to.regress = c('percent.mito','nCount_RNA','nFeature_RNA'),
  method = 'glmGamPoi'
)

anchors <- FindTransferAnchors(
  reference = reference,
  query = query,
  dims = 1:50,
  reference.reduction='pca'
)


query <- MapQuery(
  anchorset = anchors,
  query = query,
  reference = reference,
  reference.reduction = "pca", 
  reduction.model = 'umap'
)

saveRDS(object = query,file = paste0(outdir,'seurat.RDS'))


