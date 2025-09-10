# Author: Aine Fairbrother-Browne (with code from Jon B and Guillermo RP integrated)
# Date: 09/24
# Aim: To understand what the major sources of variation are in the Wood snRNA-seq data and to identify 
#      covariates that need to be controlled for in downstream analyses
# Info: 
#      Covariate selection pipeline steps: 
#        i. A priori variable filtration - removal of ID-type, 0 variance, clinical/disease variables (05a).
#        ii. Secondary variable filtration - removal of high missingness and management of co-linearity (05a).
#        iii. Sample and gene level assessment - putative covariates passed to variancePartition and to PCA (05b - this script).
#        iv. PCA-covariate correlation and generate visualisations of PCA and variancePartition results (05c)
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

# import Regina's ggplot theme and set as default
source("~/MinaRyten/Aine/wood_full/scripts/functions/get_rhr_theme.R")
theme_rhr = get.rhr.theme()

# load additional functions
source("~/MinaRyten/Aine/wood_full/scripts/functions/run_stat_test_for_all_cols_get_res.R") # chi-sq/ kruskal wallace function

# define fns. 
# 1. function to get last char of a string
get_last_char = function(x){
  
  x_split = str_split(x, "")[[1]]
  x_length = length(x_split)
  last_char = str_sub(x, start=x_length)
  
  return(last_char)
}

# 2. mode fn. to get most common char val
get_mode <- function(x) {
  unique_x <- unique(na.omit(x))
  unique_x[which.max(tabulate(match(x, unique_x)))]
}

# set variables
annot_level = 0
colinearity_cutoff = 0.7
samples_to_remove = vroom::vroom("~/MinaRyten/Aine/wood_full/data/samples_to_remove/samples_to_remove.csv", delim=",")$sample_id

# --- 1. Load metadata ---------------------------------------------------------------------------------------------------------------

# Load metadata generated in step 03a (/home/MinaRyten/Aine/wood_scrnaseq/scripts/03a-prepare_covariates.R)
metadata = vroom::vroom("~/MinaRyten/Aine/wood_full/data/metadata/05a-metadata_for_input_to_covariate_pipeline.csv")

message(paste("Ingested", ncol(metadata), "variables to the covariate selection pipeline."))

# --- 2. A priori variable filtration (step i) ---------------------------------------------------------------------------------------

# check for 0 variance cols
zero_variance_cols = c()
for(col in colnames(metadata)){
  n_unique =  metadata[[col]] %>% unique() %>% length()
  if(n_unique<2){
    message(col)
    zero_variance_cols = append(zero_variance_cols, col)
  }
}

metadata = metadata %>% 
  dplyr::relocate(sample_id, case) %>% 
  # remove ID variables (not sample_id/case)
  dplyr::select(-c(sample_id_suffix, sample_and_ba, participant_id, case, bxp_id)) %>% 
  # remove 0 variance/ label cols
  dplyr::select(-all_of(zero_variance_cols)) %>% 
  # remove disease-related variables/ tissue-related variables
  dplyr::select(-c(onset, duration, dementia, sc, csf, q, a_sn, tau, a_b, group, region)) %>% 
  # remove fail/pass indicator cols, as these are per-fastq file and do not reflect the sample as a whole well
  dplyr::select(where(~ !any(str_detect(., "pass|fail|warn"), na.rm = TRUE)))

message(paste("Reduced putative covariates to", ncol(metadata)))

metadata %>% 
  as.data.frame() %>% 
  .[1:4, ]

# --- 3a. Secondary variable filtration: missingness (step ii) -----------------------------------------------------------------------

# get NA counts for each col
na_count = data.frame(sapply(metadata, function(y) sum(length(which(is.na(y)))))) %>% 
  tibble::rownames_to_column("col") %>% 
  `colnames<-`(c("col", "na_count")) %>% 
  tibble::tibble() %>% 
  dplyr::arrange(-na_count)

na_count %>% 
  print(n=Inf)

# keep cols with 0-1 NAs
metadata = metadata %>% 
  dplyr::select(all_of(
    na_count %>% 
      dplyr::filter(na_count<2) %>% 
      dplyr::pull(col))
  ) %>% 
  dplyr::relocate(sample_id)

# PD0590P is missing bulk metadata, so mean/mode fill
metadata = metadata %>% 
  dplyr::mutate(across(where(is.numeric), ~ ifelse(is.na(.), mean(., na.rm = TRUE), .))) %>%
  dplyr::mutate(across(where(is.character), ~ ifelse(is.na(.), get_mode(.), .)))

# ensure categorical variables are defined as such
# there are 19 unique batches, but these are coded as numeric despite the numeric value being meaningless
# so code these as a character/discrete variable
metadata$batch %>% unique

metadata = metadata %>% 
  dplyr::mutate(batch = as.character(as.factor(batch)))

message(paste("Reduced putative covariates to", ncol(metadata)))

# --- 3b. Secondary variable filtration: co-linearity (step ii) --------------------------------------------------------------------------------

# Co-linearity of variables must be tested using 3 stat tests: 
# 1. Numeric/ numeric (Spearman rank correlation)
# 2. Numeric/ categorical (Kruskal-Wallace)
# 3. Categorical/ categorical (Chi-sq)

# --- 3b-1. Secondary variable filtration: co-linearity of numeric variables (step ii) -----------------------------------------------

# get names of numeric variables
numeric_cols = metadata %>%
  dplyr::select(where(is.numeric)) %>%
  colnames()

# calculate covariate/covariate correlation matrix
col_cor_mat_numeric = metadata %>%
  dplyr::select(all_of(numeric_cols)) %>%
  # center and scale before calculating correlations
  dplyr::mutate(across(where(is.numeric), .fns = ~scale(., center = T, scale = T))) %>%
  tibble::tibble() %>%
  # convert to correlation matrix
  as.matrix() %>%
  cor(., method = "spearman")

# visualise covariate/covariate correlation
col_cor_mat_numeric %>% 
  as.data.frame() %>% 
  tibble::rownames_to_column("var1") %>% 
  tidyr::pivot_longer(2:ncol(.), names_to="var2") %>% 
  dplyr::arrange(-abs(value)) %>% 
  ggplot(data=., aes(x=reorder(var1, value), y=reorder(var2, value), fill=value)) + 
  geom_tile() + 
  scale_fill_viridis(limits = c(-1,1)) + 
  theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1))

# use caret to find and remove co-linear numeric variables
caret_detected_colinears = col_cor_mat_numeric %>%
  # pass correlation matrix to caret, and calculate correlations, re-calculating each time a co-linear one is removed (selected with exact=T)
  caret::findCorrelation(., cutoff = colinearity_cutoff, names=TRUE, exact=T, verbose=TRUE)

caret_detected_colinears %>% length()

# visualise which variables will be removed and which correlations generated this decision
col_cor_mat_numeric %>% 
  as.data.frame() %>% 
  tibble::rownames_to_column("var1") %>% 
  tidyr::pivot_longer(2:ncol(.), names_to="var2") %>% 
  dplyr::arrange(-abs(value)) %>% 
  dplyr::filter(abs(value)<1) %>% 
  dplyr::mutate(removed = case_when(
    ((var1 %in% caret_detected_colinears) | (var2 %in% caret_detected_colinears)) & (abs(value)>=colinearity_cutoff) ~ "X",
    TRUE ~ ""
  )) %>% 
  ggplot(data=., aes(x=reorder(var1, value), y=reorder(var2, value), fill=value, label=removed)) + 
  geom_tile() + 
  scale_fill_viridis(limits = c(-1,1)) + 
  geom_text() + 
  theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1))

# remove the co-linear variables identified by caret from the metadata
metadata = metadata %>%
  dplyr::select(-any_of(c(caret_detected_colinears)))

# use caret to find numeric covariates that are a linear combination of other covariates
caret_detected_lincombo = metadata %>%
  dplyr::select(any_of(numeric_cols)) %>%
  dplyr::mutate(across(where(is.numeric), .fns = ~scale(., center = T, scale = T))) %>%
  as.matrix() %>%
  caret::findLinearCombos()

# remove the linear combination variables identified by caret from the metadata
metadata = metadata %>%
  dplyr::select(-any_of(caret_detected_lincombo$remove))

message(paste("Reduced putative covariates to", ncol(metadata)))

# --- 3b-2. Secondary variable filtration: co-linearity of non-numeric variables (step ii) -------------------------------------------

# get names of categorical columns
categorical_cols = metadata %>%
  dplyr::select(where(is.character) | where(is.factor), -sample_id) %>%
  colnames()

# get names of numeric columns (some have been removed, so these have changed from above)
numeric_cols = metadata %>%
  dplyr::select(where(is.numeric)) %>%
  colnames()

# if there are more than 1 cat cols, run a chi-sq correlation test between categorical variables, otherwise return empty df
if(length(categorical_cols)>1){
  chi_out = metadata %>%
    dplyr::select(all_of(categorical_cols)) %>%
    run_stat_test_for_all_cols_get_res(in_df=., stat_test="chisq", colnames1 = colnames(.), colnames(.)) %>%
    dplyr::arrange(p)
}else{
  chi_out = data.frame(var1=c(), var2=c(), p=c(), stat=c(), test_run=c())
  message(paste0("Only ", length(categorical_cols), " categorical variable, so no correlation test run."))
}

# run a kruskal-wallace chi-square between numeric and categorical variables
kruskal_out = metadata %>%
  run_stat_test_for_all_cols_get_res(in_df=., stat_test="kruskal", colnames1=numeric_cols, colnames2=categorical_cols) %>%
  dplyr::filter(var1 %in% numeric_cols) %>%
  dplyr::filter(var2 %in% categorical_cols) %>%
  dplyr::arrange(p)

# gather kruskal (cat/numeric) and chi (cat/cat) test results
cat_cor_tests = dplyr::bind_rows(
  chi_out,
  kruskal_out
) %>%
  dplyr::arrange(p) %>%
  dplyr::mutate(stat = abs(stat))

# define number of cat/numeric variables for kruskal/chi P-value adjustment
n_numeric = length(colnames(col_cor_mat_numeric))
n_categorical = length(unique(c(cat_cor_tests$var1, cat_cor_tests$var2)))

# apply P-value adjustment
cat_cor_tests = cat_cor_tests %>%
  dplyr::mutate(n_tests = case_when(
    test_run == "chisq" ~ n_categorical*n_categorical,
    test_run == "kruskal" ~ n_categorical*n_numeric,
  )) %>%
  dplyr::mutate(is_sig = case_when(
    p<(0.05/n_tests) ~ p,
    TRUE ~ NA
  ))

cor_table = cat_cor_tests %>% 
  dplyr::filter(!is.na(is_sig))

# Apply Guillermo's algo for removing max. covariates at the cat/cat cat/num step IF there any sig. correlations in cor_table:
# i) Select only the covariates with p-values <= 1e-15 (can be modified).
# ii) Remove all covariates that are co-linear with the most common covariate
# in the table.
# iii) Repeat "ii" until no more covariates are left in the co-linear table.
if(nrow(cor_table)>0){
  
  # set parameters
  collinear_p_limit = 1e-15 # The "collinear_p_limit" parameter decides on the maximum p-value two covariates must have to be considered for removal in the categorical vs categorical/numerical covariate decision. 
  collinear_filter = "bottom" # removes the less commonly collinear covariates. Thus, it removes the maximum amount of covariates possible.
  
  tmp_cor_table = cor_table %>%
    rowwise() %>%
    dplyr::mutate(test = paste0(sort(c(var1, var2)), collapse = "-")) %>%
    dplyr::distinct(test, .keep_all = T) %>%
    dplyr::filter(p <= collinear_p_limit)
  filtered_vars = c()
  top_vars = c(tmp_cor_table$var1, tmp_cor_table$var2) %>% table() %>% sort(decreasing = T)
  
  while(length(top_vars) > 0){
    if(collinear_filter == "bottom"){
      filtered_vars_iter = tmp_cor_table %>%
        dplyr::filter(var1 == names(top_vars)[1] | var2 == names(top_vars)[1]) %>%
        dplyr::select(var1, var2) %>%
        tidyr::pivot_longer(c(var1, var2)) %>%
        dplyr::filter(value != names(top_vars)[1]) %>%
        dplyr::pull(value) %>%
        unique()
      filtered_vars = c(filtered_vars, filtered_vars_iter)
    }else{
      filtered_vars = c(filtered_vars, names(top_vars)[1])
    }
    
    tmp_cor_table = tmp_cor_table %>% dplyr::filter(!var1 %in% filtered_vars & !var2 %in% filtered_vars)
    top_vars = c(tmp_cor_table$var1, tmp_cor_table$var2) %>% table() %>% sort(decreasing = T)
  }
  
  metadata = metadata %>%
    dplyr::select(-all_of(filtered_vars))
  
  message("After colinearity detection, ", ncol(metadata), " variables remain")
}else{
  message(paste0("No cat/corr or cat/cat tests were significant"))
  message("After colinearity detection, ", ncol(metadata), " variables remain")
}

message(paste("Reduced putative covariates to", ncol(metadata)))

# --- 4. Run sample and gene level assessment on each dataset (step iii) -------------------------------------------------------------

list.files(path="~/MinaRyten/Aine/wood_full/data/normalised_pseudobulk/all_tissue", 
           pattern=".rds", 
           full.names=T) %>% 
  parallel::mclapply(X=., mc.cores=length(.), FUN=function(fp){
    
    # # TEST
    # fp = "/home/MRAineFairbrotherBrowne/MinaRyten/Aine/wood_full/data/normalised_pseudobulk/all_tissue/normalised_pseudobulk_sum_annotation_level_2_opc.rds"
    
    message(paste(
      "Running variancePartition and PCA for file:",fp
    ))
    
    # --- 4a. Read in normalised pseudobulk data -----------------------------------------------------------------------------------------
    
    message("Reading data.")
    
    # normalisation carried out using: 
    # /home/MinaRyten/Aine/wood_scrnaseq/scripts/03b-normalise_counts.R
    # read in normalised pseudobulk data
    pb_dat = readRDS(fp) %>% 
      # removed _kit from this sample_id instead of '_kit5' (so this fixes the mistake in 03b-normalise_counts.R)
      `colnames<-`(gsub("C046Pkit5", "C046P", colnames(.))) %>%
      `colnames<-`(gsub("PD0687Pnew", "PD0687P", colnames(.))) %>%
      .[ ,!colnames(.) %in% samples_to_remove]

    print(dim(pb_dat))
    
    # get annotation/celltype name i.e. 'annotation_level_1_neurons' and get tissue
    level_annot_name = stringr::str_match_all(basename(fp), ".+level_(\\d_.+).rds")[[1]][,2]
    # tissue_name = stringr::str_match_all(basename(fp), ".+level_(\\d).+_(\\w).rds")[[1]][,3]
    tissue_name = "ALL"
    
    message(paste(
      "\n Annotation/level =", level_annot_name,
      "\n Tissue =", tissue_name
    ))
    
    pb_dat %>% dim()
    
    # --- 4b. Prepare putative covariates ------------------------------------------------------------------------------------------------
    
    message("Preparing putative covariates.")
    
    # scale and center variables prior to running variancePartition
    metadata_vp_input = metadata %>%
      dplyr::mutate(across(where(is.numeric), .fns = ~scale(., center = T, scale = T))) %>% 
      dplyr::filter(sample_id %in% colnames(pb_dat)) %>% 
      dplyr::arrange(factor(sample_id, levels = colnames(pb_dat))) %>% 
      tibble::column_to_rownames("sample_id")
    
    # check number and order of counts and metadata are the same
    if(
      (unique(rownames(metadata_vp_input) == colnames(pb_dat))) & (nrow(metadata_vp_input) == ncol(pb_dat))
    ){
      message("Tests passed.")
    } else{
      message("Tests failed, exiting.")
      print(unique(rownames(metadata_vp_input) == colnames(pb_dat)))
      print(nrow(metadata_vp_input) == ncol(pb_dat))
      quit()
    }
    
    # --- 4c. Run variancePartition ------------------------------------------------------------------------------------------------------
    
    message("Running variancePartition.")
    
    # generate design formula to feed into variancePartition
    vp_formula = metadata_vp_input %>%
      # get categorical covariates and wrap them in discrete (1|) notation
      dplyr::select(where(is.character)) %>%
      # remove patient and group, as we don't want these included in the covariates
      dplyr::select(-any_of(c("case"))) %>%
      colnames() %>%
      paste0("(1|", ., ")") %>%
      # A: get numerical covariates and keep them in regular continuous notation
      c(., metadata %>%
          dplyr::select(where(is.numeric)) %>%
          colnames()) %>%
      reformulate()
    
    message(paste(
      "Running variancePartition with the following formula:",
      paste0(vp_formula, collapse = " ")
    ))
    
    unique(colnames(pb_dat) == rownames(metadata_vp_input))
    ncol(metadata_vp_input)
    nrow(metadata_vp_input)
    
    # run variance partition
    # param = BiocParallel::SnowParam(10, "FORK", progressbar = TRUE)
    varPart_res = variancePartition::fitExtractVarPartModel(exprObj = pb_dat,
                                                            formula = vp_formula,
                                                            data = metadata_vp_input
                                                            )
    # BPPARAM=param)
    
    # --- 4d. Run PCA --------------------------------------------------------------------------------------------------------------------
    
    message("Running PCA.")
    
    # transpose the matrix so that rows = samples and columns = variables
    pca_res = pb_dat %>% 
      .[rowSums(. != 0) > 0, ] %>% # remove all-0 genes
      t() %>% 
      as.matrix() %>% 
      stats::prcomp(., scale = T)
    
    # --- 4e. Export objects -------------------------------------------------------------------------------------------------------------
    
    message("Exporting objects.")
    
    # save the input metadata
    metadata_vp_input %>% 
      saveRDS(., paste0("~/MinaRyten/Aine/wood_full/data/covariate_pipeline/variancePartition_putative_covariates/",
                        "variancePartition_putative_covariates_",
                        level_annot_name,
                        "_",
                        tissue_name,
                        ".rds"))
    
    # variancePartition output object
    varPart_res %>% 
      saveRDS(., paste0("~/MinaRyten/Aine/wood_full/data/covariate_pipeline/variancePartition_output/", 
                        "variancePartition_object_", 
                        level_annot_name,
                        "_",
                        tissue_name,
                        ".rds"))
    
    # PCA output object
    pca_res %>% 
      saveRDS(., paste0("~/MinaRyten/Aine/wood_full/data/covariate_pipeline/PCA_output/", 
                        "PCA_object_", 
                        level_annot_name,
                        "_",
                        tissue_name,
                        ".rds"))
    
    message("Done.")
  })

