library(org.Hs.eg.db)
library(clusterProfiler)
library(tidyverse)
library(ggsci)

# define symbol to ens ID function
convert_sym_to_ens = function(sym){
  
  df = sym %>% 
    mapIds(org.Hs.eg.db, keys=., keytype="SYMBOL", column="ENSEMBL") %>% 
    as.data.frame()
  
  df %>% 
    dplyr::filter(!is.na(.)) %>% 
    .$. %>% 
    return()
}

# define ensembl to entrez
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
  
  temp$entrez %>% 
    return()
}

var2string = function(v1) {
  deparse(substitute(v1))
}

# read in DEGs
de_res = readRDS("~/MinaRyten/Aine/wood_full/data/dreamlet/de_X_sum_all_DE.rds") %>% 
  janitor::clean_names() %>% 
  dplyr::mutate(direction = case_when(
    log_fc<0 ~ "Down-regulated", 
    log_fc>0 ~ "Up-regulated"
  )) %>% 
  tidyr::separate(coef, into=c("tissue", "comparison"), "_") %>% 
  dplyr::mutate(entrezId = ensembl_to_entrez(id))

# ad_rare_mendelian_opentargets = vroom::vroom("~/MinaRyten/ASAP/curated_gene_lists/AD_rareMendelian_openTargets.csv")
ad_common_gwas_bellenguez = vroom::vroom("~/MinaRyten/ASAP/curated_gene_lists/AD_common_gwas_bellenguez2022.csv")
#pd_common_gwas_kim = vroom::vroom("~/MinaRyten/ASAP/curated_gene_lists/PD_common_gwas_kim2023.csv")
#pd_common_gwas_nalls = vroom::vroom("~/MinaRyten/ASAP/curated_gene_lists/PD_commonGWAS_nalls2019.csv")
pd_common_gwas_nallskim = vroom::vroom("~/MinaRyten/ASAP/downstream_analyses/check_for_enrichment/PD_commonGWAS_nalls2019_kim2023_GP2_2024_union_with_ens.csv", delim = ",")

ens_table = list(
  # ad_rare_mendelian_opentargets$gene_ens,
  ad_common_gwas_bellenguez$gene_ens,
  #pd_common_gwas_kim$gene_ens,
  #pd_common_gwas_nalls$gene_ens,
  pd_common_gwas_nallskim$gene_ens
) %>% 
  `names<-`(c(
    # var2string(ad_rare_mendelian_opentargets),
    var2string(ad_common_gwas_bellenguez),
    #var2string(pd_common_gwas_kim),
    #var2string(pd_common_gwas_nalls),
    var2string(pd_common_gwas_nallskim)
  )) %>% 
  stack() %>% 
  `colnames<-`(c("geneId", "termId")) %>% 
  dplyr::as_tibble()

# create TERM2NAME df 
TERM2GENE = ens_table %>% 
  dplyr::relocate(termId, geneId) %>% 
  tidyr::drop_na() %>% 
  dplyr::mutate(entrezId = ensembl_to_entrez(geneId)) %>% 
  dplyr::select(-geneId) %>% 
  tidyr::drop_na()

# check overlap
overlap = de_res %>% 
  dplyr::select(annot_level, tissue, assay, comparison) %>% 
  dplyr::distinct() %>% 
  dplyr::group_by(annot_level, tissue, assay, comparison) %>% 
  dplyr::group_split() %>% 
  parallel::mclapply(X=., mc.cores=10, FUN=function(x){
    
    message(x)
    
    l=x$annot_level
    t=x$tissue
    a=x$assay
    c=x$comparison
    
    # l=0
    # t="ACG"
    # a="Neurons"
    # d="Up-regulated"
    # c="PDDControl"
    
    lapply(X=unique(ens_table$termId), FUN=function(x){
      
      x_ens = ens_table %>% 
        dplyr::filter(termId==x) %>% 
        dplyr::pull(geneId)

      n_overlap = de_res %>% 
        dplyr::filter(annot_level==l, tissue==t, assay==a, comparison==c) %>% 
        dplyr::filter(id %in% x_ens) %>% 
        nrow()
      
      n_overlap_sig = de_res %>% 
        dplyr::filter(annot_level==l, tissue==t, assay==a, comparison==c) %>% 
        dplyr::filter(id %in% x_ens, adj_p_val<0.05) %>% 
        nrow()
      
      data.frame(
        annot_level=l, 
        comparison=c, 
        tissue=t, 
        assay=a, 
        disease_list=x,
        n_overlap=n_overlap,
        n_overlap_sig=n_overlap_sig
      ) %>% 
        return()
      
    }) %>% 
      dplyr::bind_rows() %>% 
      return()
    
  }) %>% 
  purrr::discard(is.null) %>% 
  dplyr::bind_rows() %>% 
  dplyr::mutate(perc_disease_gene = ifelse(n_overlap>0,(n_overlap_sig/n_overlap)*100,0)) %>% 
  dplyr::as_tibble()

overlap %>% 
  vroom::vroom_write(x=., file="processed_data/overlap_disease_lists.csv")