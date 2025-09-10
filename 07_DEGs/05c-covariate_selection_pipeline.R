# Author: Aine Fairbrother-Browne (with code from Jon B and Guillermo RP integrated)
# Date: 09/24
# Aim: To understand what the major sources of variation are in the Wood snRNA-seq data and to identify 
#      covariates that need to be controlled for in downstream analyses.
# Info: 
#      Covariate selection pipeline steps: 
#        i. A priori variable filtration - removal of ID-type, 0 variance, clinical/disease variables (05a).
#        ii. Secondary variable filtration - removal of high missingness and management of co-linearity (05a).
#        iii. Sample and gene level assessment - putative covariates passed to variancePartition and to PCA (05b).
#        iv. PCA-covariate correlation and generate visualisations of PCA and variancePartition results (05c - this script)
#        v. Organise visualisations and perform final covariate selection - cut-offs tuned and covariates selected (notebook 05d)

# --- 0. Setup -----------------------------------------------------------------------------------------------------------------------

library(tidyverse)
library(SingleCellExperiment)
library(dreamlet)
library(magrittr)
library(SummarizedExperiment)
library(ggsci)
library(caret)
library(viridis)
library(patchwork)
library(Hmisc)

# import Regina's ggplot theme and set as default
source("~/MinaRyten/Aine/wood_full/scripts/functions/get_rhr_theme.R")
theme_rhr = get.rhr.theme()

# load additional functions
source("~/MinaRyten/Aine/wood_full/scripts/functions/run_stat_test_for_all_cols_get_res.R")

# set variables
annot_level = 1
cores = 5

# --- 1. Implement step iv in pipeline across datasets -------------------------------------------------------------------------------

# message("Removing existing objects")
# system("cd ~/MinaRyten/Aine/wood_scrnaseq/data/03d-final_objects/per-tissue; 
#        rm *.rds;")

# get all celltype/tissue combinations
vroom::vroom(paste0("~/MinaRyten/Aine/wood_full/data/annotation_mapping_tables/annotation",annot_level,"_to_others.csv"))$annotation_level_mapped %>%
  unique() %>%
  janitor::make_clean_names() %>%
  .[!grepl("opc", .)] %>% 
  # run PCA-covariate correlation and generation of visualisations
  parallel::mclapply(X=., mc.cores=cores, FUN=function(x){
    
    # # TEST
    # x="immune"
    
    message(paste(
      "Processing dataset:",x
    ))
    
    outfile = paste0("~/MinaRyten/Aine/wood_full/data/covariate_pipeline/final_objects/", 
                     "covariate_pipeline_final_obj", "_", 
                     "annotation_level", "_", 
                     annot_level, "_", 
                     x, ".rds")
    
    # --- 1a. Load data -------------------------------------------------------------------------------------------------------------------
    
    # read in objects output by 03d-1-covariate_selection_pipeline.R:
    # 1. variancePartition output
    # 2. PCA output
    # 3. putative covariates
    varPart = readRDS(paste0("~/MinaRyten/Aine/wood_full/data/covariate_pipeline/variancePartition_output/variancePartition_object_",annot_level,"_",x,"_ALL.rds"))
    pca = readRDS(paste0("~/MinaRyten/Aine/wood_full/data/covariate_pipeline/PCA_output/PCA_object_",annot_level,"_",x,"_ALL.rds"))
    meta = readRDS(paste0("~/MinaRyten/Aine/wood_full/data/covariate_pipeline/variancePartition_putative_covariates/variancePartition_putative_covariates_",annot_level,"_",x,"_ALL.rds")) %>% 
      tibble::rownames_to_column("sample_id")
    
    # --- 1b. Generate variancePartition visualisation and calculate summary stats --------------------------------------------------------
    
    # plot variancePartition output
    varPart_plot = varPart %>%
      variancePartition::sortCols() %>%
      variancePartition::plotVarPart(., label.angle=60)
    
    # get variancePartition summary stats
    varPart_summary_stats = varPart %>%
      as.data.frame() %>%
      tibble::rownames_to_column("gene_id") %>%
      tibble::tibble() %>%
      tidyr::pivot_longer(2:ncol(.)) %>%
      dplyr::group_by(name) %>%
      dplyr::summarise(med=median(value),
                       max=max(value)) %>%
      dplyr::filter(name!="Residuals") %>% 
      dplyr::arrange(-med)
    
    # --- 1c. Compute PC-covariate correlations -------------------------------------------------------------------------------------------
    
    # highest PC to extract - extract PCs that together explain 85% of variance
    highest_pc = which(summary(pca)$importance[3, ] > 0.85)[1]
    
    # get spearman correlation between PCs and numeric covariates [rho]
    spearman_rho_pca_numeric = pca %>%
      .[["x"]] %>%
      as.data.frame() %>%
      tibble::rownames_to_column("sample_id") %>%
      tibble::tibble() %>%
      .[, seq(highest_pc+1)] %>%
      dplyr::left_join(x=.,
                       y=meta %>% dplyr::select(sample_id, where(is.numeric)),
                       by="sample_id") %>%
      dplyr::select(-any_of(c("sample_id"))) %>%
      # dplyr::mutate(across(c(where(is.character)), as.factor),
      #               across(c(where(is.factor)), as.numeric)) %>%
      as.matrix() %>%
      Hmisc::rcorr(., type = "spearman") %>%
      .$r %>%
      as.data.frame() %>%
      dplyr::select(matches("^PC\\d+")) %>%
      # # J: need to get rid of some of the meta variables from picard that being with PCT
      # dplyr::select(-contains("PCT")) %>%
      tibble::rownames_to_column("meta_var") %>%
      dplyr::filter(!grepl("^PC\\d+", meta_var)) %>%
      tidyr::pivot_longer(2:ncol(.)) %>%
      dplyr::mutate(name = as.factor(as.numeric(gsub("PC", "", name)))) %>%
      dplyr::group_by(meta_var) %>%
      dplyr::mutate(max.PC.cor = max(value),
                    mean.PC.cor = mean(value)) %>%
      dplyr::ungroup()
    
    # get spearman correlation between PCs and numeric covariates [P-values]
    spearman_pvalue_pca_numeric = pca %>%
      .[["x"]] %>%
      as.data.frame() %>%
      tibble::rownames_to_column("sample_id") %>%
      tibble::tibble() %>%
      .[, seq(highest_pc+1)] %>%
      dplyr::left_join(x=.,
                       y=meta %>% dplyr::select(sample_id, where(is.numeric)),
                       by="sample_id") %>%
      dplyr::select(-any_of(c("sample_id"))) %>%
      # dplyr::mutate(across(c(where(is.character)), as.factor),
      #               across(c(where(is.factor)), as.numeric)) %>%
      as.matrix() %>%
      Hmisc::rcorr(., type = "spearman") %>%
      .$P %>%
      as.data.frame() %>%
      dplyr::select(matches("^PC\\d+")) %>%
      # # J: need to get rid of some of the meta variables from picard that being with PCT
      # dplyr::select(-contains("PCT")) %>%
      tibble::rownames_to_column("meta_var") %>%
      dplyr::filter(!grepl("^PC\\d+", meta_var)) %>%
      tidyr::pivot_longer(2:ncol(.)) %>%
      dplyr::mutate(name = as.factor(as.numeric(gsub("PC", "", name))))
    
    # bind P-values and rho values
    correlations_numeric = dplyr::left_join(
      x=spearman_pvalue_pca_numeric %>% dplyr::rename(p=value),
      y=spearman_rho_pca_numeric %>% dplyr::rename(stat=value),
      by=c("meta_var", "name")
    ) %>%
      dplyr::mutate(var_type = "continuous") %>%
      dplyr::group_by(meta_var) %>%
      dplyr::mutate(min.pvalue = min(p)) %>%
      dplyr::relocate(var_type, .after = last_col())
    
    # get K-W Chi-squared correlation between PCs and categorical covariates
    correlations_categorical = pca %>%
      .[["x"]] %>%
      as.data.frame() %>%
      tibble::rownames_to_column("sample_id") %>%
      tibble::tibble() %>%
      .[, seq(highest_pc+1)] %>%
      dplyr::left_join(x=.,
                       y=meta %>% dplyr::select(sample_id, !where(is.numeric)),
                       by="sample_id") %>%
      dplyr::select(-any_of(c("case", "sample_id"))) %>%
      # A: run K-W test
      run_stat_test_for_all_cols_get_res(in_df=., stat_test="kruskal") %>%
      # A: filter so var1 contains the numerics, var2 contains the categorical
      dplyr::filter(grepl("^PC\\d+", var1),
                    !grepl("^PC\\d+", var2)) %>%
      dplyr::select(-test_run) %>%
      dplyr::rename(name=var1, meta_var=var2) %>%
      dplyr::group_by(meta_var) %>%
      dplyr::mutate(max.PC.cor = max(stat),
                    mean.PC.cor = mean(stat),
                    min.pvalue = min(p)) %>%
      dplyr::mutate(name = as.factor(as.numeric(gsub("PC", "", name)))) %>%
      dplyr::ungroup() %>%
      dplyr::relocate(meta_var, name) %>%
      dplyr::mutate(var_type="categorical")
    
    # generate dataframe containing covariate-PC correlations for all putative covariates (numeric & categorical)
    correlations_main = dplyr::bind_rows(
      correlations_numeric,
      correlations_categorical
    ) %>%
      dplyr::mutate(name = as.factor(name)) %>%
      dplyr::mutate(
        padj = p.adjust(p, method="fdr"),
        is.sig = case_when(
          padj<0.05 ~ stat,
          TRUE ~ NA_integer_
        )) %>%
      # A: add variance explained
      dplyr::left_join(
        x=.,
        y=summary(pca) %>%
          .$importance %>%
          as.data.frame() %>%
          t() %>%
          as.data.frame() %>%
          tibble::rownames_to_column("name") %>%
          janitor::clean_names() %>%
          dplyr::mutate(name = as.factor(gsub("PC", "", name))),
        by="name"
      ) %>%
      # calculate weighted correlation - weights PC-covariate correlation coefficient by the % variance explained by that PC
      dplyr::mutate(variance_weighted_correlation = proportion_of_variance*(stat**2)) %>%
      # add adjusted P-value - adjust for number of putative covariates
      dplyr::group_by(name) %>%
      dplyr::mutate(padj = p.adjust(p, "fdr")) %>%
      dplyr::ungroup()
    
    # --- 1d. generate PC-covariate heatmaps ----------------------------------------------------------------------------------------------
    
    # plot PC-cov correlations
    heatmap_PCA_numeric_covariates = correlations_main %>%
      dplyr::mutate(name = as.numeric(as.character(name))) %>%
      dplyr::filter(var_type=="continuous") %>%
      dplyr::arrange(name) %>%
      {
        ggplot(data=., aes(x=name, y=reorder(meta_var, max.PC.cor), fill=stat)) +
          geom_tile() +
          geom_text(aes(label = round(is.sig, 2)), color = "black", size = 1.75) +
          scale_fill_gradient2(low = "#075AFF",
                               mid = "white",
                               high = "#FF0000",
                               limits = c(-1, 1),
                               name = "Spearman's rho") +
          xlab("PC") +
          ylab("") +
          theme(axis.text.x = element_text(angle=0, vjust = 0.5, hjust=0.5),
                legend.position = "right") +
          scale_x_continuous(breaks=c(1:max(.$name))) +
          facet_wrap(~var_type, scales = "free")
      }
    
    heatmap_PCA_categorical_covariates = correlations_main %>%
      dplyr::mutate(name = as.numeric(as.character(name))) %>%
      dplyr::filter(var_type=="categorical") %>%
      dplyr::arrange(name) %>%
      {
        ggplot(data=., aes(x=name, y=reorder(meta_var, max.PC.cor), fill=stat)) +
          geom_tile() +
          geom_text(aes(label = ifelse(p<0.05, round(p, 2), "")), color = "black", size = 1.75) +
          scale_fill_gradient2(low = "#075AFF",
                               mid = "white",
                               high = "#FF0000",
                               name = "K-W chi-squared") +
          xlab("PC") +
          ylab("") +
          theme(axis.text.x = element_text(angle=0, vjust = 0.5, hjust=0.5),
                legend.position = "right") +
          scale_x_continuous(breaks=c(1:max(.$name))) +
          facet_wrap(~var_type, scales = "free")
      }
    
    # use patchwork to get cateorical and numeric/cont heatmaps in one figure
    heatmap_main_fig = (heatmap_PCA_numeric_covariates / heatmap_PCA_categorical_covariates)
    
    # --- 1e. Generate variance explained plot --------------------------------------------------------------------------------------------
    
    # plot prop of variance explained
    barplot_variance_explained_pca = summary(pca) %>%
      .$importance %>%
      as.data.frame() %>%
      t() %>%
      as.data.frame() %>%
      tibble::rownames_to_column("PC") %>%
      janitor::clean_names() %>%
      dplyr::mutate(pc = as.numeric(gsub("PC", "", pc))) %>%
      dplyr::filter(pc %in% seq(1,highest_pc)) %>% 
      {
        ggplot(data=.) +
          geom_col(aes(x=pc, y=proportion_of_variance)) +
          geom_line(aes(x=pc, y=cumulative_proportion), linetype=2) +
          xlab("Principle component #") +
          ylab("Proportion of variance") +
          scale_x_continuous(breaks = round(seq(min(.$pc), max(.$pc), by = 1),1)) +
          theme(axis.text.x = element_text(angle=0)) +
          ylim(0,1)
        # geom_text(aes(x=pc, y=proportion_of_variance+0.025, label=paste0(round(proportion_of_variance, 3)*100, "%")), size=2.5)
      }
    
    # --- 1f. Write out final object (list of objects) ------------------------------------------------------------------------------------
    
    final_object_list = list(
      varPart_plot, # 1. variancePartition distribution (violin plot)
      heatmap_main_fig, # 2. PCA-covariate correlation (heatmap)
      barplot_variance_explained_pca, # 3. Prop variance explained by each PC (barplot)
      correlations_main # 5. PC-meta correlations with variance explained and variance_weighted_correlation (table)
    )
    
    final_object_list %>% 
      saveRDS(., outfile)
    
  })



