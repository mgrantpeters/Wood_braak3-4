# Date: 09/24
# Author: Aine Fairbrother-Browne

# conda activate /home/MRAineFairbrotherBrowne/miniconda3/envs/R

library(clusterProfiler)
library(AnnotationDbi)
library(dreamlet)
library(magrittr)
library(org.Hs.eg.db)
library(EnsDb.Hsapiens.v86)
library(tidyverse)
library(ReactomePA)
library(SummarizedExperiment)

ensembl_to_entrez = function(ensembl){
  
  library(org.Hs.eg.db)
  library(EnsDb.Hsapiens.v86)
  
  temp = ensembl %>%
    data.frame(ensembl = .) %>%
    dplyr::mutate(ensdb_v86 = mapIds(EnsDb.Hsapiens.v86,
                                     keys=ensembl, #Column containing Ensembl gene ids
                                     column="ENTREZID",
                                     keytype="GENEID",
                                     multiVals="first")) %>%
    dplyr::mutate(org_hs_ed_db = mapIds(org.Hs.eg.db,
                                        keys=ensembl, #Column containing Ensembl gene ids
                                        column="ENTREZID",
                                        keytype="ENSEMBL",
                                        multiVals="first")) %>%
    dplyr::as_tibble() %>%
    dplyr::mutate(ensdb_v86 = as.character(ensdb_v86),
                  org_hs_ed_db = as.character(org_hs_ed_db)) %>%
    dplyr::mutate(entrez = dplyr::case_when(
      ((is.na(ensdb_v86)==T) & (is.na(org_hs_ed_db)==F)) ~ org_hs_ed_db,
      ((is.na(ensdb_v86)==F) & (is.na(org_hs_ed_db)==T)) ~ ensdb_v86,
      ((is.na(ensdb_v86)==T) & (is.na(org_hs_ed_db)==T)) ~ NA_character_,
      TRUE ~ ensdb_v86
    ))
  
  # temp = clusterProfiler::bitr(geneID=ensembl, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Hs.eg.db, drop = T)
  # 
  # data.frame(ENSEMBL=ensembl) %>% 
  #   dplyr::left_join(x=., y=temp, by="ENSEMBL") %>% 
  #   dplyr::pull(ENTREZID) %>% 
  #   return()
  
  temp$entrez %>%
    return()
}

enrich_res = list.files("~/MinaRyten/Aine/wood_full/data/dreamlet", "de.+annotation.+rds", full.names=T) %>%
  lapply(X=., FUN=function(x){
    
    # x="/home/MRAineFairbrotherBrowne/MinaRyten/Aine/wood_full/data/dreamlet/de_X_sum_annotation_level_1.rds"
    
    annot_level = stringr::str_match_all(string = basename(x), pattern=".+level_(\\d)")[[1]][,2]

    de = readRDS(x)
    
    de_res = assayNames(de) %>%
      lapply(X=., FUN=function(x){

        message(x)

        # if coef has one _ in it, filter for this and get all DEGs
        coefNames(de)[stringr::str_count(coefNames(de), "_")==1] %>%

          lapply(X=., FUN=function(y){

            de %>%
              dreamlet::topTable(., coef=y, number=Inf) %>%
              dplyr::as_tibble() %>%
              dplyr::filter(assay==x) %>%
              dplyr::mutate(coef=y) %>%
              return()

          }) %>%
          dplyr::bind_rows() %>%
          janitor::clean_names() %>%
          dplyr::as_tibble() %>%
          return()

      }) %>%
      `names<-`(SummarizedExperiment::assayNames(de)) %>%
      dplyr::bind_rows(., .id = 'assay') %>%
      dplyr::mutate(direction = case_when(
        log_fc<0 ~ "Down-regulated",
        log_fc>0 ~ "Up-regulated",
      ))
    
    all_in = de_res %>%
      dplyr::mutate(entrez_id = ensembl_to_entrez(de_res$id)) %>%
      tidyr::drop_na(entrez_id) %>% 
      tidyr::separate(coef, into=c("tissue", "comparison"))
    
    all_in %>% 
      dplyr::group_by(assay, tissue) %>% 
      dplyr::group_split() %>% 
      lapply(X=., FUN=function(x){
        
        sig_in = x %>%
          dplyr::filter(adj_p_val<0.05)  # Removed logFC filtering from Aine's script
        
        assay_i = x %>% dplyr::pull(assay) %>% unique()
        
        message(paste("Running for", assay_i, " |  universe size = ", length(unique(x$entrez_id))))
        
        if(nrow(sig_in)>0){
          
          universe_ = x %>% dplyr::select(entrez_id) %>% dplyr::distinct() %>% dplyr::pull(entrez_id)
          
          message("Running KEGG enrichment")
          kegg_res = clusterProfiler::compareCluster(
            geneClusters = entrez_id ~ assay + tissue + comparison + direction,
            data = sig_in,
            fun = enrichKEGG,
            universe = universe_) %>%
            head(., n=Inf)

          message("Running GO enrichment")
          go_res = clusterProfiler::compareCluster(
            geneClusters = entrez_id ~ assay + tissue + comparison + direction,
            data = sig_in,
            fun = clusterProfiler::enrichGO,
            OrgDb = org.Hs.eg.db::org.Hs.eg.db,
            ont = c("ALL"),
            universe = universe_) %>%
            head(., n=Inf)
          
          message("Running REACTOME enrichment")
          reactome_res = clusterProfiler::compareCluster(
            geneClusters = entrez_id ~ assay + tissue + comparison + direction,
            data = sig_in,
            fun = ReactomePA::enrichPathway,
            universe = universe_) %>%
            head(., n=Inf)

          final_res = vector("list", 3)

          if(!is.null(kegg_res)){
            final_res[[1]] = kegg_res %>%
              dplyr::mutate(source="KEGG", annot_level=annot_level, assay=assay_i)
          } else{
            return(NULL)
          }

          if(!is.null(go_res)){
            final_res[[2]] = go_res %>%
              dplyr::mutate(source="GO", annot_level=annot_level, assay=assay_i)
          } else{
            return(NULL)
          }

          if(!is.null(reactome_res)){
            final_res[[3]] = reactome_res %>%
              dplyr::mutate(source="Reactome", annot_level=annot_level, assay=assay_i)
          } else{
            return(NULL)
          }

          final_res = final_res %>%
            .[which(lapply(., is.null) == FALSE)] %>%
            dplyr::bind_rows()
          
          if(nrow(final_res)>0){
            return(final_res)
          } else{
            message("No pws enriched across GO, REACTOME, KEGG, returning a NULL entry.")
            return(NULL)
          }
            
        } else{
          message("No sig. DEGs, exiting and returning a NULL entry.")
          return(NULL)
        }
        
      }) %>% 
      .[which(lapply(., is.null) == FALSE)] %>% 
      dplyr::bind_rows() %>% 
      return()
    
  })

enrich_res %>%
  .[which(lapply(., is.null) == FALSE)] %>% 
  dplyr::bind_rows() %>% 
  vroom::vroom_write(x=., file=paste0("07-clusterProfiler_res_logfc1.csv"), ",")
