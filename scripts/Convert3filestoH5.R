#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)

library(Seurat)
library(DropletUtils)


dir <- 'AHKW078/Distal_CD45neg'

system(paste0('gzip ',dir,'/STARSolo/Solo.out/Gene/filtered/*'))

cts <- Read10X(paste0(dir,'/STARSolo/Solo.out/Gene/filtered/'))


write10xCounts(path=paste0(dir,'/STARSolo/gene_filtered.h5'),
               x=cts, 
               type='HDF5'
               )


