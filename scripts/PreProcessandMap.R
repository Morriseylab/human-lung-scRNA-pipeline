#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)


library(Seurat)
library(tidyverse)
library(patchwork)
library(scExtras)

name<-args[1]
org<-args[2]
input10x = args[3]
outdir=args[4]
ref<-args[5]

# Make Dirs ---------------------------------------------------------------
dir.create(outdir,recursive = T,showWarnings = F)
plotdir <- paste0(outdir,'/plots')
dir.create(plotdir,showWarnings = F)
qcdir <-paste0(outdir,'/qc')
dir.create(qcdir,showWarnings = F)


query= RunQC(dir=outdir,org=org,name=name,files=input10x ,filter=T, doubletdetection = F,UpperMitoCutoff=10)
reference <- readRDS(ref)

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

saveRDS(object = query,file = paste0(outdir,'/seurat.RDS'))

