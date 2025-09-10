# Author: Aine Fairbrother-Browne
# Date: 09/24
# Aim: pseudobulk sn data (per-sample, per-assay)

# conda activate /home/MRAineFairbrotherBrowne/miniconda3/envs/R

# --- 0. Setup ----------------------------------------------------------------------------------------------------------------------

library(tidyverse)
library(SingleCellExperiment)
library(dreamlet)
library(magrittr)
library(SummarizedExperiment)
library(parallel)

# --- 1. Read SCE file ----------------------------------------------------------------------------------------------------------------------

sce = readRDS("~/MinaRyten/Aine/wood_full/data/sce/wood_fulldata_wrangled_lowQualRemoved_annot1Fixed.sce")

intersect(
  colData(sce)$index,
  colnames(sce)
) %>% length()

intersect(assay(sce, "raw_counts") %>% colnames(),
          colnames(sce)) %>% 
  length()

intersect(assay(sce, "raw_counts") %>% colnames(),
          colData(sce)$index) %>% 
  length()

# --- 2. Apply pseudobulking across assays ----------------------------------------------------------------------------------------------------------------------

out_dir = "~/MinaRyten/Aine/wood_full/data/pseudobulk/all_tissue/"

# write out pseudobulked data as per annotation level SCE objects (*.rds) and per annotation level per celltype dataframes (*.csv)
c("annotation_level_1",	"annotation_level_2",	"annotation_level_3", "annotation_level_4") %>%
  parallel::mclapply(X=., mc.cores=2, FUN=function(x){

    message(paste(
      "1. Running pseudobulk on", x
    ))

    # run pseudobulk function on sce annotation level x
    pseudobulk = dreamlet::aggregateToPseudoBulk(x=sce,
                                                 assay="raw_counts",
                                                 fun="sum",
                                                 sample_id="sample_id",
                                                 cluster_id=x, # this is the cell type annotation the counts will be summed on i.e. annotation_level_0 etc.
                                                 scale=FALSE,
                                                 verbose=TRUE,
                                                 checkValues=TRUE)

    # write out RDS object for pseudobulk of annotation level x (contains all cell types, accessed using assay(obj, "cell_type"))
    rds_filename = paste0("pseudobulk_sum", "_", x, ".rds")
    saveRDS(pseudobulk, paste0(out_dir, rds_filename))

    message(paste(
      "2. Generating per cell-type pseudobulk files for", x
    ))

    # write out .csv file for pseudobulk of each cell type in annotation level x
    for(i in 1:length(assayNames(pseudobulk))){

      message(paste(
        "2a. Generating cell-type pseudobulk file for", assayNames(pseudobulk)[i]
      ))

      pseudo_csv_filename = assayNames(pseudobulk)[i] %>%
        # remove special chars from assayName i.e. spaces/ fwd slashes
        janitor::make_clean_names(.) %>%
        # add .csv ext instead of .rds in already formed rds_filename
        paste0("_", ., ".csv") %>%
        gsub(".rds", ., rds_filename)

      assay(pseudobulk, i) %>%
        as.matrix() %>%
        t() %>%
        as.data.frame() %>%
        tibble::rownames_to_column("sample_id") %>%
        tibble::tibble() %>%
        vroom::vroom_write(., paste0(out_dir, pseudo_csv_filename), delim=",")
    }
  })
