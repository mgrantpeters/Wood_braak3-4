# Author: Aine Fairbrother-Browne
# Date: 09/24
# Aim: To wrangle, integrate and prepare metadata for ingestion into the covariate correction pipeline

# --- 0. Setup -----------------------------------------------------------------------------------------------------------------------

library(tidyverse)
library(SingleCellExperiment)
library(dreamlet)
library(magrittr)
library(SummarizedExperiment)
library(ggsci)

# --- 1. Read SCE file ---------------------------------------------------------------------------------------------------------------

sce = readRDS("~/MinaRyten/Aine/wood_full/data/sce/wood_fulldata_wrangled_lowQualRemoved_annot1Fixed.sce")

# --- 2. Load metadata sources -------------------------------------------------------------------------------------------------------

# load coldata from SCE object - covariates that Melissa loaded into the original AnnData obj
# notes: 
# the region column has levels C, F, P, corresponding to:
#   F = Frontal cortex
#   P = Parietal cortex
#   C = Cingulate cortex
sce_coldata = colData(sce) %>% 
  as.data.frame() %>% 
  dplyr::as_tibble()

# --- 3. Write out metadata ----------------------------------------------------------------------------------------------------------

sce_coldata %>% 
  as.data.frame() %>% 
  .[1:4,]

# remove cell-level metadata, retaining only sample-level (we can't correct for cell level due to pseudobulking)
final_meta = sce_coldata %>% 
  dplyr::select(-c(
    doublet_scores, predicted_doublets, n_genes_by_counts, log1p_n_genes_by_counts, total_counts, total_counts, 
    log1p_total_counts, total_counts_hb, log1p_total_counts_hb, pct_counts_hb, total_counts_mt, log1p_total_counts_mt, 
    pct_counts_mt, total_counts_rp, log1p_total_counts_rp, pct_counts_rp, markers_neutro_score, s_score, g2m_score, 
    phase, contains("leiden"), contains("annotation_level"), cell_id, old_cell_id
  )) %>% 
  dplyr::distinct() %>% 
  dplyr::mutate(sample_id = gsub("_.*", "", sample_id_suffix)) %>% 
  dplyr::relocate(sample_id)

final_meta %>% 
  vroom::vroom_write(x=., file="~/MinaRyten/Aine/wood_full/data/metadata/05a-metadata_for_input_to_covariate_pipeline.csv", ",")
