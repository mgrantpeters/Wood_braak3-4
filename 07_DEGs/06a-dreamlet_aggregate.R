# Date: 09/24
# Author: Aine Fairbrother-Browne

library(tidyverse)
library(org.Hs.eg.db)
library(parallel)
library(dreamlet)

# --- Define fns ---------------------------------------------------------------------------------------------------

compile.de = function(wd, pattern){
  
  # load and wrangle dreamlet DE output - get all DEGs
  de.compiled = list.files(wd, pattern, full.names=T) %>% 
    
    parallel::mclapply(X=., mc.cores=4, FUN=function(f){
      
      annot_level = stringr::str_match_all(string = basename(f), pattern=".+level_(\\d)")[[1]][,2]
      
      de = readRDS(f)
      
      coefnames = dreamlet::coefNames(de)[stringr::str_count(dreamlet::coefNames(de), "_")==1]
      
      coefnames %>% 
        lapply(X=., FUN=function(coef){
          de %>% 
            dreamlet::topTable(., coef=coef, number=Inf, adjust.method="fdr") %>% 
            as.data.frame() %>% 
            tibble::tibble() %>% 
            dplyr::select(assay, ID, logFC, P.Value, adj.P.Val) %>% 
            dplyr::mutate(coef = coef, annot_level = annot_level) %>% 
            return()
        }) %>% 
        dplyr::bind_rows() %>% 
        return()
      
    }) %>% 
    dplyr::bind_rows()
  
  # map assay names to 'clean' equivalent
  assay_map = data.frame(
    assay = de.compiled$assay %>% unique(),
    clean_assay = de.compiled$assay %>% unique() %>% janitor::make_clean_names()
  )
  
  de.compiled = de.compiled %>% 
    dplyr::left_join(x=., y=assay_map, by="assay") %>% 
    dplyr::relocate(assay, clean_assay)
  
  # map gene symbols to DE results (they're ENS at present)
  gene.id.map = de.compiled %>%
    dplyr::pull(ID) %>%
    AnnotationDbi::mapIds(org.Hs.eg.db, keys=., keytype="ENSEMBL", column="SYMBOL") %>%
    stack() %>%
    `colnames<-`(c("gene_sym", "gene_ens"))

  de.compiled = de.compiled %>% 
    dplyr::left_join(x=., y=gene.id.map, by=c("ID"="gene_ens")) %>%
    dplyr::relocate(ID, gene_sym) %>%
    dplyr::mutate(gene_sym = case_when(
      is.na(gene_sym) ~ ID,
      TRUE ~ gene_sym
    )) %>%
    dplyr::distinct()
  
  saveRDS(de.compiled, paste0(wd, "/", "de_X_sum_all_DE.rds"))
  
}

compile.details = function(wd, pattern){
  
  # load and wrangle dreamlet DE output - get all DEGs
  list.files(wd, pattern, full.names=T) %>% 
    
    parallel::mclapply(X=., mc.cores=length(.), FUN=function(f){
      
      annot_level = stringr::str_match_all(string = basename(f), pattern=".+level_(\\d)")[[1]][,2]
      
      de = readRDS(f)
      
      gene_counts = lapply(X=assayNames(de), FUN=function(x){
        assay(de, x)$Amean %>%
          length()
      }) %>% 
        `names<-`(assayNames(de)) %>% 
        data.frame() %>% 
        t() %>% 
        as.data.frame() %>% 
        tibble::rownames_to_column("assay") %>% 
        dplyr::rename(n_genes = V1) %>% 
        dplyr::mutate(assay = janitor::make_clean_names(assay))
      
      details = details(de) %>% 
        dplyr::as_tibble() %>% 
        dplyr::select(assay, n_retain) %>%
        dplyr::mutate(assay = janitor::make_clean_names(assay)) %>% 
        dplyr::rename(n_samples_retained = n_retain) %>% 
        dplyr::filter(assay %in% janitor::make_clean_names(assayNames(de)))
      
      final = details %>% 
        dplyr::left_join(x=., y=gene_counts, by="assay") %>% 
        dplyr::mutate(annot_level=annot_level)
      
      return(final)
      
    }) %>% 
    dplyr::bind_rows() %>% 
    vroom::vroom_write(x=., paste0(wd, "/", "de_X_sum_details.csv"))
}

# --- Implement -----------------------------------------------------------------------------------------------------

compile.de(
  wd="~/MinaRyten/Aine/wood_full/data/dreamlet/",
  pattern="de_X_sum_annotation_level_\\d"
)