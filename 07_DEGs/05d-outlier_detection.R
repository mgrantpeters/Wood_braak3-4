# Author: Aine Fairbrother-Browne
# Date: 11/23
# Aim: To wrangle, integrate and prepare metadata for ingestion into the covariate correction pipeline

# --- 0. Setup -----------------------------------------------------------------------------------------------------------------------

library(here)
library(ggplot2)
library(tidyverse)
library(stringr)
library(DESeq2)
library(ggsci)
library(variancePartition)
library(patchwork)
library(factoextra)
library(WGCNA)

source("~/MinaRyten/Aine/wood_full/scripts/functions/get_rhr_theme.R")
theme_rhr = get.rhr.theme()

get_last_char = function(x){
  x_split = str_split(x, "")[[1]]
  x_length = length(x_split)
  last_char = str_sub(x, start=x_length)
  return(last_char)
}

# --- 1. Load data -------------------------------------------------------------------------------------------------------------------

# To understand whether outlier samples exist in the Wood snRNAseq data

pb = list.files("~/MinaRyten/Aine/wood_full/data/normalised_pseudobulk/all_tissue", pattern="level_2", full.names=T) %>% 
  lapply(X=., FUN=function(x){
    readRDS(x) %>% 
      as.data.frame() %>% 
      tibble::rownames_to_column("gene_id") %>% 
      tidyr::pivot_longer(2:ncol(.), names_to="sample_id", values_to="norm_fpm") %>% 
      dplyr::mutate(
        cell_type = stringr::str_match_all(string=x, pattern="_\\d_(.+)\\.rds")[[1]][,2],
        tissue = "all_tissues"
      ) %>% 
      return()
  }) %>% 
  dplyr::bind_rows() %>% 
  dplyr::mutate(sample_id = gsub("kit5", "", sample_id)) %>% 
  dplyr::mutate(sample_id = gsub("new", "", sample_id))

# get covariates
grouping = vroom::vroom("~/MinaRyten/Aine/wood_full/data/metadata/05a-metadata_for_input_to_covariate_pipeline.csv") %>% 
  dplyr::select(sample_id, group)

covs = vroom::vroom("~/MinaRyten/Aine/wood_full/data/metadata/05d-final_covariates.csv") %>% 
  dplyr::mutate(batch = as.factor(as.character(batch))) %>% 
  dplyr::mutate(sample_id = gsub("_kit5", "", sample_id)) %>% 
  dplyr::left_join(x=., y=grouping, by="sample_id") %>% 
  dplyr::select(-tissue)

# --- 2. Implement WGCNA/ sample connectivity method -------------------------------------------------------------------------------------------------------------------

# Adapted from https://github.com/dhglab/Broad-transcriptomic-dysregulation-across-the-cerebral-cortex-in-ASD/blob/master/code/01_RNAseqProcessing/01_02_A_CountsProcessing.R

sample_connectivity_res = pb %>% 
  dplyr::group_by(cell_type, tissue) %>% 
  dplyr::group_split() %>% 
  
  parallel::mclapply(X=., mc.cores=length(.), FUN=function(df){
    
    df = pb %>%
      dplyr::group_by(cell_type, tissue) %>%
      dplyr::group_split() %>%
      .[[1]]
    
    cell_type = unique(df$cell_type)
    tissue = unique(df$tissue)
    
    mat = df %>% 
      dplyr::select(sample_id, gene_id, norm_fpm) %>% 
      dplyr::distinct() %>% 
      tidyr::pivot_wider(names_from = "gene_id", values_from = "norm_fpm") %>% 
      tibble::column_to_rownames("sample_id") %>% 
      as.matrix() %>% 
      t()
    
    # Calculate biweight midcorrelation matrix between samples
    bicor.mat = mat %>% 
      WGCNA::bicor(., use = "pairwise.complete.obs")
    
    # Calculate signed adjacency matrix between samples
    normadj.mat = (0.5+0.5*bicor.mat)^2
    
    # Calculate connectivity
    connectivity.mat = fundamentalNetworkConcepts(normadj.mat)
    ku = connectivity.mat$Connectivity
    z.ku = (ku-mean(ku))/sqrt(var(ku))
    
    z.ku %>% 
      stack() %>% 
      dplyr::rename(sample_id = ind, 
                    wgcna_connectivity_z_score = values) %>% 
      dplyr::relocate(sample_id, wgcna_connectivity_z_score) %>% 
      dplyr::mutate(cell_type=cell_type, 
                    tissue=tissue) %>% 
      return(.)
  }) %>% 
  dplyr::bind_rows()

# how many datasets is the case+region an outlier in, /7 (one per cell type)
sample_connectivity_res %>% 
  dplyr::filter(wgcna_connectivity_z_score < -2) %>% 
  dplyr::group_by(sample_id) %>% 
  dplyr::summarise(n=n()) %>% 
  dplyr::arrange(-n)

# --- 3. Implement PCA/ Z-score method -------------------------------------------------------------------------------------------------------------------

# run PCA on each of the datasets separately
pca.ls = pb %>% 
  
  # get grouping info by binding metadata table
  dplyr::left_join(x=.,
                   y=covs %>%
                     dplyr::select(group, sample_id) %>%
                     dplyr::distinct(),
                   by=c("sample_id")) %>%
  dplyr::relocate(sample_id, group, tissue, cell_type) %>%
  
  # group and split to run PCA on each dataset separately
  dplyr::group_by(cell_type, tissue) %>% 
  dplyr::group_split() %>% 
  
  # format each dataset for PCA
  parallel::mclapply(X=., mc.cores=length(.), FUN=function(df){
    
    df %>% 
      dplyr::mutate(index = paste0(cell_type, ".", sample_id, ".", group)) %>% 
      dplyr::select(-any_of(c("cell_type", "sample_id", "tissue", "group", "log2fpm", "residual_fpm"))) %>% 
      dplyr::relocate(index, gene_id, norm_fpm) %>% 
      tidyr::pivot_wider(names_from = gene_id, values_from=norm_fpm) %>% 
      tibble::column_to_rownames("index") %>% 
      # run PCA
      .[ , apply(., 2, function(x) !any(is.na(x)))] %>% 
      stats::prcomp(.)
  })

all.pca.res = parallel::mclapply(X=c(1:length(pca.ls)), mc.cores=length(pca.ls), FUN=function(i){
  
  pca.ls[[i]]$x %>% 
    tibble::as_tibble(., rownames="index") %>% 
    tidyr::separate(col=index, into=c("cell_type", "sample_id", "group"), sep="\\.", remove=F) %>% 
    dplyr::rowwise() %>% 
    dplyr::mutate(tissue = get_last_char(sample_id)) %>% 
    dplyr::mutate(dataset = paste0(cell_type, "--", tissue)) %>% 
    dplyr::relocate(index, dataset) %>% 
    return(.)
  
}) %>% 
  dplyr::bind_rows()

# --- 4. Collect the results of the PCA run on each of the 14 datasets, marking Z-score outliers -------------------------------------------------------------------------------------------------------------------

# Collect the results of the PCA run on each of the 14 datasets, marking Z-score outliers (>3 in PCX or PCY for each PC pair).  
# Visualise these outliers for each dataset in PC1-10 space. 
# calculate Z-score to identify outliers
# Jon's code for Z-score calculation (per-PC)

plot.outlier.pca = function(all.pca.res, data.set){
  
  dat.and.plots = list(
    c("PC1", "PC2"),
    c("PC3", "PC4"),
    c("PC5", "PC6"),
    c("PC7", "PC8"),
    c("PC9", "PC10")) %>% 
    
    lapply(X=., FUN=function(i){
      
      PCX = i[[1]]
      PCY = i[[2]]
      
      X.zscore.colname = paste0("PC_z_score_", PCX)
      X.PC.colname = paste0("value_", PCX)
      
      Y.zscore.colname = paste0("PC_z_score_", PCY)
      Y.PC.colname = paste0("value_", PCY)
      
      dat = all.pca.res %>% 
        
        dplyr::filter(dataset==data.set) %>%
        
        dplyr::select(index:group, paste0("PC", c(1:10))) %>% 
        tidyr::pivot_longer(PC1:ncol(.), names_to="PC", values_to="value") %>% 
        dplyr::group_by(dataset, PC) %>% 
        dplyr::mutate(PC_mean=mean(value), 
                      PC_sd=sd(value)) %>% 
        dplyr::ungroup() %>% 
        dplyr::mutate(PC_z_score=abs((value-PC_mean)/PC_sd)) %>% 
        dplyr::select(-PC_mean, -PC_sd) %>% 
        tidyr::pivot_wider(names_from="PC", values_from=c("value", "PC_z_score")) %>% 
        
        # is the sample a Z-score outlier according to PCX or PCY? If so, mark TRUE
        dplyr::mutate(is_outlier = case_when(
          (!!sym(X.zscore.colname) > 3) | (!!sym(Y.zscore.colname) > 3) ~ TRUE,
          TRUE ~ FALSE
        )) %>%
        
        # label point with sample_id if TRUE
        dplyr::mutate(label = case_when(
          is_outlier==TRUE ~ sample_id, 
          is_outlier==FALSE ~ ""
        ))
      
      plot = dat %>% 
        {
          # plot
          ggplot(data=., aes(x=!!sym(X.zscore.colname), y=!!sym(Y.zscore.colname), colour=is_outlier, label=label)) + 
            geom_point(alpha=0.5, size=2) + 
            facet_wrap(~dataset, scales = "free") + 
            # add % buffer zone to axes so labels/ points don't get cut-off by boundaries
            scale_y_continuous(expand = expansion(mult = 0.4)) +
            scale_x_continuous(expand = expansion(mult = 0.4)) + 
            geom_text(size=3, hjust=0, vjust=1) +
            labs(x=PCX, y=PCY)
        }
      
      return(plot)
      
    })
  
  patch.plot = dat.and.plots %>% 
    patchwork::wrap_plots(.) + 
    patchwork::plot_layout(guides="collect")
  
  return(patch.plot)
  
}

z.score.plot.out = all.pca.res$dataset %>% unique() %>%
  lapply(X=., function(x){
    plot.outlier.pca(
      all.pca.res=all.pca.res,
      data.set=x)
  })

patchwork::wrap_plots(z.score.plot.out, ncol=1)

get.zscore.pca.dat = function(all.pca.res, data.set){
  
  all_dat = list(
    c("PC1", "PC2"),
    c("PC3", "PC4"),
    c("PC5", "PC6"),
    c("PC7", "PC8"),
    c("PC9", "PC10")) %>%
    
    lapply(X=., FUN=function(i){
      
      PCX = i[[1]]
      PCY = i[[2]]
      
      X.zscore.colname = paste0("PC_z_score_", PCX)
      X.PC.colname = paste0("value_", PCX)
      
      Y.zscore.colname = paste0("PC_z_score_", PCY)
      Y.PC.colname = paste0("value_", PCY)
      
      dat = all.pca.res %>%
        
        dplyr::filter(dataset==data.set) %>%
        
        dplyr::select(index:group, paste0("PC", c(1:10))) %>%
        tidyr::pivot_longer(PC1:ncol(.), names_to="PC", values_to="value") %>%
        dplyr::group_by(dataset, PC) %>%
        dplyr::mutate(PC_mean=mean(value),
                      PC_sd=sd(value)) %>%
        dplyr::ungroup() %>%
        dplyr::mutate(PC_z_score=abs((value-PC_mean)/PC_sd)) %>%
        dplyr::select(-PC_mean, -PC_sd) %>%
        tidyr::pivot_wider(names_from="PC", values_from=c("value", "PC_z_score")) %>%
        
        # is the sample a Z-score outlier according to PCX or PCY? If so, mark TRUE
        dplyr::mutate(is_outlier = case_when(
          (!!sym(X.zscore.colname) > 3) | (!!sym(Y.zscore.colname) > 3) ~ TRUE,
          TRUE ~ FALSE
        )) %>%
        
        # label point with sample_id if TRUE
        dplyr::mutate(label = case_when(
          is_outlier==TRUE ~ sample_id,
          is_outlier==FALSE ~ ""
        ))
      
      return(dat)
      
    })
  
  all_dat %>%
    dplyr::bind_rows() %>%
    return()
}

z.score.dat.out = all.pca.res$dataset %>% unique() %>% 
  lapply(X=., function(x){
    get.zscore.pca.dat(
      all.pca.res=all.pca.res,
      data.set=x)
  }) %>% 
  dplyr::bind_rows() %>% 
  dplyr::rowwise() %>% 
  dplyr::mutate(tissue = get_last_char(sample_id)) %>% 
  dplyr::select(dataset, sample_id, cell_type, tissue, contains("PC_z_score_"))

z.score.dat.out %>% 
  head()

# --- 5. Visualise the WGCNA connectivity Z-scores and the PCA -derived Z-scores for each dataset -------------------------------------------------------------------------------------------------------------------

# visualise WGCNA connectivity distribution
plot_wgcna_connectivity = z.score.dat.out %>% 
  dplyr::mutate(dataset = paste0(cell_type, "--", tissue)) %>% 
  dplyr::mutate(tissue="all_tissues") %>% 
  dplyr::left_join(
    x=., 
    y=sample_connectivity_res,
    by=c("sample_id", "cell_type", "tissue")
  ) %>% 
  dplyr::distinct() %>% 
  
  ggplot(data=., aes(x=dataset, y=wgcna_connectivity_z_score)) + 
  geom_violin(fill='grey') + 
  geom_jitter(alpha=0.1) + 
  geom_hline(yintercept = -2)

plot_pca_zscore = z.score.dat.out %>% 
  dplyr::left_join(
    x=., 
    y=sample_connectivity_res,
    by=c("sample_id", "cell_type", "tissue")
  ) %>% 
  dplyr::distinct() %>% 
  dplyr::relocate(dataset, sample_id, cell_type, tissue, wgcna_connectivity_z_score) %>% 
  tidyr::pivot_longer(contains("PC_z_score_"), values_to="PC_z_score", names_to="PC") %>% 
  dplyr::mutate(PC = as.numeric(gsub("PC_z_score_PC", "", PC))) %>% 
  
  ggplot(data=., aes(x=dataset, y=PC_z_score, size=PC)) + 
  geom_point(alpha=0.1) + 
  scale_size_identity() + 
  geom_hline(yintercept = 4)

plot_wgcna_connectivity | plot_pca_zscore

# --- 6. Visualise the WGCNA connectivity Z-scores and the PCA -derived Z-scores per sample -------------------------------------------------------------------------------------------------------------------

pca_z_score_cutoff = 3
wgcna_connectivity_cutoff = -2

z.score.dat.out %>% 
  dplyr::left_join(
    x=., 
    y=sample_connectivity_res,
    by=c("sample_id", "cell_type", "tissue")
  ) %>% 
  dplyr::distinct() %>% 
  dplyr::relocate(dataset, sample_id, cell_type, tissue, wgcna_connectivity_z_score) %>% 
  tidyr::pivot_longer(contains("PC_z_score_"), values_to="PC_z_score", names_to="PC") %>% 
  dplyr::mutate(PC = as.numeric(gsub("PC_z_score_PC", "", PC))) %>% 
  dplyr::mutate(label = case_when(
    ((PC_z_score > pca_z_score_cutoff) & (wgcna_connectivity_z_score < wgcna_connectivity_cutoff)) ~ sample_id,
    TRUE~""
  )) %>% 
  
  ggplot(data=.) + 
  geom_point(aes(x=sample_id, y=PC_z_score), colour="lightblue") +
  geom_point(aes(x=sample_id, y=wgcna_connectivity_z_score), colour="red", alpha=0.1) +
  # geom_text(aes(x=sample_id, y=PC_z_score, label=label)) + 
  # geom_text(aes(x=sample_id, y=wgcna_connectivity_z_score, label=label)) + 
  # facet_grid(tissue~"", scales = "free_x") +
  geom_hline(yintercept = pca_z_score_cutoff, linetype=2, alpha=0.5) +
  geom_hline(yintercept = wgcna_connectivity_cutoff, linetype=2, alpha=0.5) +
  scale_y_continuous() + 
  ylab("Z-score")

# --- 7. plot number of datasets sample is an outlier in -------------------------------------------------------------------------------------------------------------------

# To try to understand whether any outliers are considered outliers across multiple datasets
  
z.score.dat.out %>% 
  dplyr::mutate(tissue="all_tissues") %>% 
  dplyr::left_join(
    x=., 
    y=sample_connectivity_res,
    by=c("sample_id", "cell_type", "tissue")
  ) %>% 
  dplyr::mutate(is_outlier = case_when(((PC_z_score_PC1 > 3 |
                                           PC_z_score_PC2 > 3 |
                                           PC_z_score_PC3 > 3 |
                                           PC_z_score_PC4 > 3 |
                                           PC_z_score_PC5 > 3 |
                                           PC_z_score_PC6 > 3 |
                                           PC_z_score_PC7 > 3 |
                                           PC_z_score_PC8 > 3 |
                                           PC_z_score_PC9 > 3 |
                                           PC_z_score_PC10 > 3) & (wgcna_connectivity_z_score < -2)) ~ 'TRUE',
                                       TRUE ~ 'FALSE')) %>% 
  dplyr::select(dataset, sample_id, tissue, cell_type, is_outlier) %>% 
  dplyr::mutate(is_outlier = as.character(is_outlier)) %>% 
  dplyr::distinct() %>% 
  dplyr::group_by(sample_id, tissue) %>% 
  # count how many are outliers for each case/ brain region
  dplyr::summarise(
    n = n(),
    across(is_outlier, list(true = ~ sum(. == "TRUE"),
                            false  = ~ sum(. == "FALSE")))) %>%
  dplyr::ungroup() %>% 
  
  ggplot(., aes(x=reorder(sample_id, is_outlier_true), y=is_outlier_true)) + 
  geom_col() + 
  ylim(0,5) + 
  geom_hline(yintercept=3, linetype=2) + 
  xlab("Sample") + 
  ylab("No. of celltypes the sample is an outlier in")

# --- 8. Write file containing sample_id of outliers -------------------------------------------------------------------------------------------------------------------

# no outliers detected for Wood full sample set (as of 0924)
samples.to.remove = z.score.dat.out %>% 
  dplyr::mutate(tissue="all_tissues") %>% 
  dplyr::left_join(
    x=., 
    y=sample_connectivity_res,
    by=c("sample_id", "cell_type", "tissue")
  ) %>% 
  dplyr::mutate(is_outlier = case_when(((PC_z_score_PC1 > 3 |
                                           PC_z_score_PC2 > 3 |
                                           PC_z_score_PC3 > 3 |
                                           PC_z_score_PC4 > 3 |
                                           PC_z_score_PC5 > 3 |
                                           PC_z_score_PC6 > 3 |
                                           PC_z_score_PC7 > 3 |
                                           PC_z_score_PC8 > 3 |
                                           PC_z_score_PC9 > 3 |
                                           PC_z_score_PC10 > 3) & (wgcna_connectivity_z_score < -2)) ~ 'TRUE',
                                       TRUE ~ 'FALSE')) %>% 
  dplyr::select(dataset, sample_id, tissue, cell_type, is_outlier) %>% 
  dplyr::mutate(is_outlier = as.character(is_outlier)) %>% 
  dplyr::distinct() %>% 
  dplyr::group_by(sample_id, tissue) %>% 
  # count how many are outliers for each case/ brain region
  dplyr::summarise(
    n = n(),
    across(is_outlier, list(true = ~ sum(. == "TRUE"),
                            false  = ~ sum(. == "FALSE")))) %>%
  dplyr::ungroup() %>% 
  dplyr::filter(is_outlier_true>3) %>% 
  dplyr::pull(sample_id)

samples.to.remove %>% 
  data.frame(sample_id = .) %>% 
  vroom::vroom_write("~/MinaRyten/Aine/wood_full/data/outliers/outliers.csv")

