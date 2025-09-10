# Author: Aine Fairbrother-Browne
# Date: 09/24
# Aim: To normalise the pseudobulked count data for assessment of covariates in the covariate assessment pipeline

# conda activate /home/MRAineFairbrotherBrowne/miniconda3/envs/R

# --- 0. Setup -----------------------------------------------------------------------------------------------------------------------

library(tidyverse)
library(SingleCellExperiment)
library(magrittr)
library(SummarizedExperiment)
library(data.table)
library(DESeq2)

# --- 1. Split pseudobulk by tissue ---------------------------------------------------------------------------------------------------

# level 1
list.files("~/MinaRyten/Aine/wood_full/data/pseudobulk/all_tissue", "level_1.+csv", full.names=T) %>% 
  parallel::mclapply(X=., mc.cores=length(.), FUN=function(f){
    
    # f = "/home/MRAineFairbrotherBrowne/MinaRyten/Aine/wood_full/data/pseudobulk/all_tissue/pseudobulk_sum_annotation_level_1_opc.csv"
    
    pb = vroom::vroom(f)
    
    outdir = "/home/drihome/MRAineFairbrotherBrowne/MinaRyten/Aine/wood_full/data/normalised_pseudobulk/all_tissue/"
    outname = paste0("normalised_",gsub(".csv", ".rds", basename(f)))
    outfile = paste0(outdir, outname)
    
    pseudo_t = pb %>% 
      tibble::column_to_rownames("sample_id") %>% 
      .[, colSums(. != 0) > 0] %>% # remove genes with all-0
      data.table::data.table() %>% 
      data.table::transpose()
    
    new_rownames = pb %>%
      tibble::column_to_rownames("sample_id") %>%
      .[, colSums(. != 0) > 0] %>%
      colnames()
    
    colnames(pseudo_t) = pb$sample_id
    
    pseudo_t = pseudo_t %>% 
      as.data.frame() %>% 
      `rownames<-`(new_rownames)
    
    meta = colnames(pseudo_t) %>% data.frame(sample_id=.)
    
    # normalisation procedure as recommended in the variancePartition vignette
    # create DESeq2 object from pseudobulk count matrix and metadata
    dds = DESeq2::DESeqDataSetFromMatrix(countData = pseudo_t,
                                         colData = meta,
                                         design = ~ 1)
    
    # estimate library size correction scaling factors
    dds = DESeq2::estimateSizeFactors(dds)
    
    # identify genes that pass expression cutoff - then filter out genes that are extremely lowly expressed
    # here, the standard cut-off set is fpm>1 in 50% or more of samples
    isexpr = rowSums(DESeq2::fpm(dds)>1) >= 0.5 * ncol(dds)
    
    # # compute log2 Fragments Per Million
    quantlog.counts = log2(DESeq2::fpm(dds)[isexpr,] + 1)
    quantlog.counts = log2(DESeq2::fpm(dds) + 1)
    
    saveRDS(quantlog.counts, outfile)
    
  })

# level 2
list.files("~/MinaRyten/Aine/wood_full/data/pseudobulk/all_tissue", "level_2.+csv", full.names=T) %>% 
  parallel::mclapply(X=., mc.cores=length(.), FUN=function(f){
    
    # f = "/home/MRAineFairbrotherBrowne/MinaRyten/Aine/wood_full/data/pseudobulk/all_tissue/pseudobulk_sum_annotation_level_1_opc.csv"
    
    pb = vroom::vroom(f)
    
    outdir = "/home/drihome/MRAineFairbrotherBrowne/MinaRyten/Aine/wood_full/data/normalised_pseudobulk/all_tissue/"
    outname = paste0("normalised_",gsub(".csv", ".rds", basename(f)))
    outfile = paste0(outdir, outname)
    
    pseudo_t = pb %>% 
      tibble::column_to_rownames("sample_id") %>% 
      .[, colSums(. != 0) > 0] %>% # remove genes with all-0
      data.table::data.table() %>% 
      data.table::transpose()
    
    new_rownames = pb %>%
      tibble::column_to_rownames("sample_id") %>%
      .[, colSums(. != 0) > 0] %>%
      colnames()
    
    colnames(pseudo_t) = pb$sample_id
    
    pseudo_t = pseudo_t %>% 
      as.data.frame() %>% 
      `rownames<-`(new_rownames)
    
    meta = colnames(pseudo_t) %>% data.frame(sample_id=.)
    
    # normalisation procedure as recommended in the variancePartition vignette
    # create DESeq2 object from pseudobulk count matrix and metadata
    dds = DESeq2::DESeqDataSetFromMatrix(countData = pseudo_t,
                                         colData = meta,
                                         design = ~ 1)
    
    # estimate library size correction scaling factors
    dds = DESeq2::estimateSizeFactors(dds)
    
    # identify genes that pass expression cutoff - then filter out genes that are extremely lowly expressed
    # here, the standard cut-off set is fpm>1 in 50% or more of samples
    isexpr = rowSums(DESeq2::fpm(dds)>1) >= 0.5 * ncol(dds)
    
    # # compute log2 Fragments Per Million
    quantlog.counts = log2(DESeq2::fpm(dds)[isexpr,] + 1)
    quantlog.counts = log2(DESeq2::fpm(dds) + 1)
    
    saveRDS(quantlog.counts, outfile)
    
  })

# level 3
list.files("~/MinaRyten/Aine/wood_full/data/pseudobulk/all_tissue", "level_3.+csv", full.names=T) %>% 
  parallel::mclapply(X=., mc.cores=length(.), FUN=function(f){
    
    # f = "/home/MRAineFairbrotherBrowne/MinaRyten/Aine/wood_full/data/pseudobulk/all_tissue/pseudobulk_sum_annotation_level_1_opc.csv"
    
    pb = vroom::vroom(f)
    
    outdir = "/home/drihome/MRAineFairbrotherBrowne/MinaRyten/Aine/wood_full/data/normalised_pseudobulk/all_tissue/"
    outname = paste0("normalised_",gsub(".csv", ".rds", basename(f)))
    outfile = paste0(outdir, outname)
    
    pseudo_t = pb %>% 
      tibble::column_to_rownames("sample_id") %>% 
      .[, colSums(. != 0) > 0] %>% # remove genes with all-0
      data.table::data.table() %>% 
      data.table::transpose()
    
    new_rownames = pb %>%
      tibble::column_to_rownames("sample_id") %>%
      .[, colSums(. != 0) > 0] %>%
      colnames()
    
    colnames(pseudo_t) = pb$sample_id
    
    pseudo_t = pseudo_t %>% 
      as.data.frame() %>% 
      `rownames<-`(new_rownames)
    
    meta = colnames(pseudo_t) %>% data.frame(sample_id=.)
    
    # normalisation procedure as recommended in the variancePartition vignette
    # create DESeq2 object from pseudobulk count matrix and metadata
    dds = DESeq2::DESeqDataSetFromMatrix(countData = pseudo_t,
                                         colData = meta,
                                         design = ~ 1)
    
    # estimate library size correction scaling factors
    dds = DESeq2::estimateSizeFactors(dds)
    
    # identify genes that pass expression cutoff - then filter out genes that are extremely lowly expressed
    # here, the standard cut-off set is fpm>1 in 50% or more of samples
    isexpr = rowSums(DESeq2::fpm(dds)>1) >= 0.5 * ncol(dds)
    
    # # compute log2 Fragments Per Million
    quantlog.counts = log2(DESeq2::fpm(dds)[isexpr,] + 1)
    quantlog.counts = log2(DESeq2::fpm(dds) + 1)
    
    saveRDS(quantlog.counts, outfile)
    
  })

# level 4
list.files("~/MinaRyten/Aine/wood_full/data/pseudobulk/all_tissue", "level_4.+csv", full.names=T) %>% 
  parallel::mclapply(X=., mc.cores=length(.), FUN=function(f){
    
    # f = "/home/MRAineFairbrotherBrowne/MinaRyten/Aine/wood_full/data/pseudobulk/all_tissue/pseudobulk_sum_annotation_level_1_opc.csv"
    
    pb = vroom::vroom(f)
    
    outdir = "/home/drihome/MRAineFairbrotherBrowne/MinaRyten/Aine/wood_full/data/normalised_pseudobulk/all_tissue/"
    outname = paste0("normalised_",gsub(".csv", ".rds", basename(f)))
    outfile = paste0(outdir, outname)
    
    pseudo_t = pb %>% 
      tibble::column_to_rownames("sample_id") %>% 
      .[, colSums(. != 0) > 0] %>% # remove genes with all-0
      data.table::data.table() %>% 
      data.table::transpose()
    
    new_rownames = pb %>%
      tibble::column_to_rownames("sample_id") %>%
      .[, colSums(. != 0) > 0] %>%
      colnames()
    
    colnames(pseudo_t) = pb$sample_id
    
    pseudo_t = pseudo_t %>% 
      as.data.frame() %>% 
      `rownames<-`(new_rownames)
    
    meta = colnames(pseudo_t) %>% data.frame(sample_id=.)
    
    # normalisation procedure as recommended in the variancePartition vignette
    # create DESeq2 object from pseudobulk count matrix and metadata
    dds = DESeq2::DESeqDataSetFromMatrix(countData = pseudo_t,
                                         colData = meta,
                                         design = ~ 1)
    
    # estimate library size correction scaling factors
    dds = DESeq2::estimateSizeFactors(dds)
    
    # identify genes that pass expression cutoff - then filter out genes that are extremely lowly expressed
    # here, the standard cut-off set is fpm>1 in 50% or more of samples
    isexpr = rowSums(DESeq2::fpm(dds)>1) >= 0.5 * ncol(dds)
    
    # # compute log2 Fragments Per Million
    quantlog.counts = log2(DESeq2::fpm(dds)[isexpr,] + 1)
    quantlog.counts = log2(DESeq2::fpm(dds) + 1)
    
    saveRDS(quantlog.counts, outfile)
    
  })


