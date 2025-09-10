# /home/drihome/MRAineFairbrotherBrowne/MinaRyten/Aine/wood_scrnaseq/scripts/05-dreamlet_eachTissue_ControlPD.R
# Date: 09/24
# Author: Aine Fairbrother-Browne
#
# Aim: To run the dreamlet pipeline on the Wood scRNA-seq data
# Description: Based on the workflow detailed in the dreamlet tutorial: 
# https://gabrielhoffman.github.io/dreamlet/articles/dreamlet.html#voom-for-pseudobulk
#
# Note: dreamlet version 0.99.16, 2023-07-04, from Github (DiseaseNeurogenomics/dreamlet@68caeca)
# conda activate /home/MRAineFairbrotherBrowne/miniconda3/envs/R

# set pseudobulk variables
assay = "X" # this is values the aggregation will be performed on, could be cpm, or X (raw counts)
fun = "sum" # this is how the counts will be aggregated
obj_name = paste0(assay, "_", fun)

## --- 0. Setup ---------------------------------------------------------------------------------------------------------------------------------------------------------------------

library(tidyverse)
library(SingleCellExperiment)
library(dreamlet)
library(zellkonverter)
library(reticulate)
library(anndata)
library(magrittr)
library(SummarizedExperiment)
library(parallel)
library(AnnotationDbi)
library(org.Hs.eg.db)

# set base dir
base_dir = "~/MinaRyten/Aine/wood_full/"

# import Regina's ggplot theme and set as default
source("~/MinaRyten/Aine/wood_full/scripts/functions/get_rhr_theme.R")
theme_rhr = get.rhr.theme()

## --- 1. Define DE function - to run over each pseudobulk level --------------------------------------------------------------------------------------------------------------------------------

run.de = function(pseudobulk_rds_filepath, assay="X", fun="sum", annot_level="1", min.count=10){
  
  # # TEST
  # pseudobulk_rds_filepath="~/MinaRyten/Aine/wood_full/data/pseudobulk/all_tissue/"
  # assay="X"
  # fun="sum"
  # annot_level="1"
  # min.count=10
  
  outname = paste0(assay, "_", fun, "_", "annotation_level", "_", annot_level)
  
  message(paste("Running DE for the configuration:", outname))
  
  ## --- 3a. Read & wrangle pseudobulk ---------------------------------------------------------------------------------------------------------------------------------------------------------------------
  
  # read in pseudobulk SCE object (generated above) for annotation_level_1
  pseudobulk = readRDS(paste0(pseudobulk_rds_filepath, paste0("pseudobulk","_",fun,"_annotation_level_",annot_level,".rds")))
  
  # one 'assay' can now be accessed per cell type
  SummarizedExperiment::assayNames(pseudobulk) %>% 
    print()
  
  assay(pseudobulk, assayNames(pseudobulk)[1]) %>% dim()
  dreamlet::colData(pseudobulk) %>% dim()
  
  ## --- 3b. Remove sample outliers from SCE object ---------------------------------------------------------------------------------------------------------------------------------------------------------------------
  
  # # Wood 55 samples had no sample outliers, but did have low quality samples to remove
  # # read these in
  # outliers = vroom::vroom("~/MinaRyten/Aine/wood_scrnaseq/data/03c-samples_to_remove/samples_to_remove.csv", delim=",")$sample_id
  # 
  # print(paste("Removing", length(outliers), "low quality samples."))
  # 
  # # filter outliers out of slots in the pseudobulk SCE object
  # pseudobulk = pseudobulk[,!(colnames(pseudobulk) %in% outliers)]
  # 
  # metadata(pseudobulk)$aggr_means = metadata(pseudobulk)$aggr_means %>% 
  #   dplyr::filter((!sample_id %in% gsub("_.+", "", outliers)))
  
  ## --- 3d. Deal with sample naming ---------------------------------------------------------------------------------------------------------------------------------------------------------------------

  # remove _kit5 label
  colnames(pseudobulk) = gsub("kit5|_kit5", "", colnames(pseudobulk))
  colnames(pseudobulk) = gsub("new|_new", "", colnames(pseudobulk))

  metadata(pseudobulk)$aggr_means = metadata(pseudobulk)$aggr_means %>%
    dplyr::mutate(sample_id = gsub("kit5|_kit5", "", sample_id)) %>%
    dplyr::mutate(sample_id = gsub("new|_new", "", sample_id))

  names(int_colData(pseudobulk)$n_cells) = gsub("kit5|_kit5", "", names(int_colData(pseudobulk)$n_cells))
  names(int_colData(pseudobulk)$n_cells) = gsub("new|_new", "", names(int_colData(pseudobulk)$n_cells))
  
  ## --- 3e. Set covariates and group as colData ---------------------------------------------------------------------------------------------------------------------------------------------------------------------
  
  # read in final covariates
  covs = vroom::vroom("~/MinaRyten/Aine/wood_full/data/metadata/05d-final_covariates.csv") %>% 
    dplyr::mutate(batch = as.factor(as.character(batch))) %>% 
    dplyr::mutate(sample_id = gsub("kit5|_kit5", "", sample_id)) %>% 
    dplyr::mutate(sample_id = gsub("new|_new", "", sample_id))
  
  dreamlet::colData(pseudobulk) = dreamlet::colData(pseudobulk) %>% 
    as.data.frame() %>% 
    tibble::rownames_to_column("sample_id") %>% 
    dplyr::as_tibble() %>% 
    dplyr::mutate(sample_id = gsub("kit5|_kit5", "", sample_id)) %>% 
    dplyr::mutate(sample_id = gsub("new|_new", "", sample_id)) %>% 
    janitor::clean_names() %>%
    dplyr::select(sample_id, group) %>% 
    dplyr::left_join(x=., y=covs, by="sample_id") %>% 
    dplyr::relocate(sample_id, tissue, group) %>% 
    as.data.frame() %>% 
    # tibble::column_to_rownames("sample_id") %>% 
    S4Vectors::DataFrame()
  
  colnames(pseudobulk) = dreamlet::colData(pseudobulk)$sample_id
  
  colnames(pseudobulk) %>% length()
  covs$sample_id %>% length()
  
  intersect(
    colnames(pseudobulk),
    covs$sample_id
  ) %>% length()
  
  ## --- 3f. Add variable that defines group/region for norm & modelling purposes ---------------------------------------------------------------------------------------------------------------------------------------------------------------------
  
  # add new col that defines group + region 
  dreamlet::colData(pseudobulk) = dreamlet::colData(pseudobulk) %>% 
    as.data.frame() %>% 
    dplyr::mutate(groupRegion = paste0(group, ".", tissue)) %>% 
    # use this function to convert the df to an S4 object, as this is the only object type that can be stored in an SCE slot
    S4Vectors::DataFrame()
  
  ## --- 3g. Deal with gene names ---------------------------------------------------------------------------------------------------------------------------------------------------------------------
  
  # at present, these are a mix of ENS and SYM - this is problematic for downstream GO analysis
  # so, set the ENS gene_ids as the row names instead of the mixed names
  rownames(pseudobulk) = rowData(pseudobulk)$gene_ids
  
  ## --- 3h. Voom normalisation ---------------------------------------------------------------------------------------------------------------------------------------------------------------------
  
  # set up formula for normalisation
  # include variables to be controlled for as well as variables to be measured
  norm.form = ~ (1|sex) + (1|batch) + (1|downsampled) + (1|merged) + (1|redo) + scale(total_deduplicated_percentage) + (1|groupRegion)
  
  unique(colnames(pseudobulk) == rownames(dreamlet::colData(pseudobulk)))
  
  # voom-style normalization - this step may take ~ 5-6 mins per cell type
  # Normalize and apply voomWithDreamWeights to handle random effects
  if(!file.exists(paste0("~/MinaRyten/Aine/wood_full/data/dreamlet/", "voomNorm_",outname,".rds"))){
    
    res.proc = dreamlet::processAssays(
      sceObj = pseudobulk,
      formula = norm.form,
      min.count = min.count
    )
    
    # save just res.proc object
    saveRDS(res.proc, paste0("~/MinaRyten/Aine/wood_full/data/dreamlet/", "voomNorm_",outname,".rds"))
    
  }else{
    res.proc = readRDS(paste0("~/MinaRyten/Aine/wood_full/data/dreamlet/", "voomNorm_",outname,".rds"))
  }
  
  plotVoom(res.proc) %>% 
    ggsave(
      filename = paste0("voomPlot_", outname, ".png"),
      path = "~/MinaRyten/Aine/wood_full/data/dreamlet/",
      plot = .,
      device = "png",
      scale = 1,
      width = 10,
      height = 7,
      units = c("in"),
      dpi = 300
    )
  
  ## --- 3i. Differential expression (per-region) ---------------------------------------------------------------------------------------------------------------------------------------------------------------------
  
  # DE performed on the voom normalized data without an intercept term (~ 0).
  # Instead estimate the mean expression in stimulated, controls and then set Diff 
  # to the difference between the two. Group here is modeled as a random effect, 
  # as this is the thing we're measuring.
  
  # define contrasts
  contrasts = c(
    C_PDControl = "groupRegionPD.C - groupRegionControl.C",
    P_PDControl = "groupRegionPD.P - groupRegionControl.P",
    F_PDControl = "groupRegionPD.F - groupRegionControl.F"
  )
  
  # form = ~ 0 + groupRegion + (1|sex) + (1|batch) + scale(gex_q30_bases_in_read_2) + scale(gex_reads_mapped_antisense_to_gene) + scale(gex_reads_mapped_confidently_to_transcriptome) +
  #   scale(percent_gc)
  form = ~ 0 + groupRegion + (1|sex) + (1|batch) + (1|downsampled) + (1|merged) + (1|redo) + scale(total_deduplicated_percentage)
  
  n_assays_before_filtration = length(assayNames(res.proc))
  
  # remove cell types that don't have enough samples to make all region/diagnosis comparisons. Must have at least one sample in each region/diagnosis group. 
  # error message: Error in eBayes(fit, robust = robust, trend = isCounts) : Need Amean component in fit to estimate trend
  assays_to_remove_1 = lapply(X=assayNames(res.proc), FUN=function(x){
    
    # this finds unique diagnosis/region pairs in each assay of res.proc and counts
    dreamlet::colData(res.proc) %>% 
      as.data.frame() %>% 
      # tibble::rownames_to_column("sample_id") %>% 
      dplyr::as_tibble() %>% 
      dplyr::filter(sample_id %in% (assay(res.proc, x)$targets %>% rownames())) %>% 
      dplyr::select(tissue, group) %>% 
      dplyr::distinct() %>% 
      dplyr::group_by(tissue, group) %>% 
      dplyr::group_split() %>% 
      length()
    
  }) %>% 
    `names<-`(assayNames(res.proc)) %>% 
    # check if counts sum to double the number of contrasts - this is the no. of comparisons being made
    .[.!=(length(contrasts)*2)] %>% 
    names()
  
  # must also have more samples than terms in the model, so remove assays with samples <= model terms
  n_model_terms = form %>% paste0(., collapse="") %>% stringr::str_split(., " \\+ ") %>% .[[1]] %>% .[2:length(.)] %>% length()
  
  res_proc_summary = details(res.proc) %>% 
    dplyr::as_tibble() %>% 
    dplyr::select(assay, n_retained) %>% 
    dplyr::arrange(n_retained)
  
  assays_to_remove_2 = res_proc_summary %>% 
    dplyr::filter(n_retained <= n_model_terms) %>% 
    dplyr::pull(assay)
  
  # filter the assays before running dreamlet
  res.proc = res.proc[!(assayNames(res.proc) %in% c(assays_to_remove_1, assays_to_remove_2))]
  
  if(!file.exists(paste0("~/MinaRyten/Aine/wood_full/data/dreamlet/failed_eBayes_",outname,".csv"))){
    
    # pre-screen of assays (has to be done this way as assays are processed iteratively within dreamlet function)
    # if eBayes errors out and can't be used due to gene number too small, get the name of this assay for removal in next step
    assays_that_fail_eBayes_de = assayNames(res.proc) %>% 
      parallel::mclapply(X=., mc.cores=length(.), FUN=function(x){
        
        test_de = try(
          dreamlet::dreamlet(
            x = res.proc[(assayNames(res.proc) %in% c(x))],
            formula =  form,
            contrasts = contrasts, 
            robust = TRUE,
            quiet = FALSE,
            use.eBayes = TRUE
          ))
        
        if("try-error" %in% class(test_de)){
          return(x)
        } else{
          return(NULL)
        }
      }) %>% 
      unlist()
    
    # write out list of these 
    data.frame(failed_ebayes=assays_that_fail_eBayes_de) %>% 
      vroom::vroom_write(x=., paste0("~/MinaRyten/Aine/wood_full/data/dreamlet/failed_eBayes_",outname,".csv"))
    
  } else{
    assays_that_fail_eBayes_de = vroom::vroom(paste0("~/MinaRyten/Aine/wood_full/data/dreamlet/failed_eBayes_",outname,".csv"))$failed_ebayes
  }
  
  # write out list of all failed assays
  data.frame(all_assays_failed=unique(c(assays_to_remove_1, assays_to_remove_2, assays_that_fail_eBayes_de))) %>% 
    vroom::vroom_write(x=., paste0("~/MinaRyten/Aine/wood_full/data/dreamlet/assay_fail_report_",outname,".csv"))
  
  # remove assays failing eBayes
  res.proc = res.proc[!(assayNames(res.proc) %in% c(assays_that_fail_eBayes_de))]
  
  # report n assays removed
  n_assays_after_filtration = length(assayNames(res.proc))
  
  message(paste(
    "Number of assays before filtration: ", n_assays_before_filtration, 
    "Number of assays after filtration: ", n_assays_after_filtration
  ))
  
  message("Running dreamlet DE...")
  
  # res.proc = res.proc[!(assayNames(res.proc) %in% c("M2 Macrophages"))]
  
  de = dreamlet::dreamlet(
      x = res.proc,
      formula =  form,
      contrasts = contrasts, 
      robust = TRUE,
      quiet = FALSE,
      use.eBayes = TRUE
    )
  
  ## --- 3j. Save DE object -------------------------------------------------------------------------------------------------------------------------------------------------------------------
  
  message("Writing DE outdata")
  
  # save just DE object
  saveRDS(de, paste0("~/MinaRyten/Aine/wood_full/data/dreamlet/", "de_",outname,".rds"))
  
}

## --- 4. Implement run.de ----------------------------------------------------------------------------------------------------------------------------------------------------------------------

parallel::mclapply(X=c("1", "2", "3", "4"), mc.cores=4, FUN=function(annot_level){
  run.de(pseudobulk_rds_filepath="~/MinaRyten/Aine/wood_full/data/pseudobulk/all_tissue/",
         annot_level=annot_level)
})

message("DE done.")

## --- 5. Session info ----------------------------------------------------------------------------------------------------------------------------------------------------------------------

# Session info
library("sessioninfo")
options(width = 120)
session_info()

## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
