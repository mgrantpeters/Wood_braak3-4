library(tidyverse)

# Define enrichment function using fisher's test
gene_enrichment<-function(sig_genes, background, test_gene_list){
  
  sig_in<-length(sig_genes[sig_genes %in% test_gene_list])
  sig_out<-length(sig_genes[!sig_genes %in% test_gene_list])
  back_in<-length(background[background %in% test_gene_list])
  back_out<- length(background[!background %in% test_gene_list])
  
  dat <- data.frame("In Gene Set" = c(sig_in, sig_out),
                    "Not In Gene Set" = c(back_in, back_out),
                    row.names = c("Significantly Spliced Genes", "Background Genes"),
                    stringsAsFactors = FALSE
  )
  
  test <- fisher.test(dat, alternative = "greater")
  
  
  res<-tibble(p_value=test$p.value, odds_ratio=test$estimate, sig_ratio=sig_in/sig_out, bg_ratio=back_in/back_out)
  
  return(res)
  
}

# read in DEGs and create extra column (up or down reg)
de_res = readRDS("~/MinaRyten/Aine/wood_full/data/dreamlet/de_X_sum_all_DE.rds") %>% 
  janitor::clean_names() %>% 
  dplyr::mutate(direction = case_when(
    log_fc<0 ~ "Down-regulated", 
    log_fc>0 ~ "Up-regulated"
  )) %>% 
  tidyr::separate(coef, into=c("tissue", "comparison"), "_") #%>% 
  #dplyr::mutate(entrezId = ensembl_to_entrez(id))

for (disease in c('PD', 'AD')){
    if (disease=='AD'){
        disease_genes = vroom::vroom("~/MinaRyten/ASAP/curated_gene_lists/AD_common_gwas_bellenguez2022.csv")
    } else if(disease=='PD'){
        disease_genes = vroom::vroom("~/MinaRyten/ASAP/downstream_analyses/check_for_enrichment/PD_commonGWAS_nalls2019_kim2023_GP2_2024_union_with_ens.csv", delim = ",")
    }
    sig = de_res %>% 
      dplyr::filter(adj_p_val<0.05)
        

    results = de_res %>% 
      dplyr::filter(annot_level=="3") %>% 
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

        all = de_res %>% 
          dplyr::filter(annot_level==l, tissue==t, assay==a, comparison==c) # %>% 
          #tidyr::drop_na(entrezId)

        sig = all %>% 
          dplyr::filter(adj_p_val<0.05)

        message(length(intersect(unique(sig$gene_sym), disease_genes$gene_sym)))

        if( (intersect(unique(sig$gene_sym), disease_genes$gene_sym) %>% length()) != 0 ){

            gene_enrichment(unique(sig$gene_sym), unique(de_res$gene_sym), disease_genes$gene_sym) %>% 
            as.data.frame() %>% 
            dplyr::mutate(annot_level=l, tissue=t, assay=a, comparison=c) %>% 
            return()

        }else{
          return(NULL)
          message("none sig.")
        }
      }) %>% 
      purrr::discard(is.null)

    # collate
    results = results %>% 
      dplyr::bind_rows() %>% 
      janitor::clean_names() %>% 
      dplyr::select(annot_level, tissue, assay, comparison, p_value) %>% 
      dplyr::as_tibble() %>% 
      dplyr::group_by(annot_level, tissue, assay, comparison) %>% 
      dplyr::mutate(p_adjust = p.adjust(p_value)) %>% 
      as.data.frame()
    results$disease_enrichment = disease
    
    results = results %>%
        group_by(disease_enrichment) %>% 
        mutate(p_adjust = p.adjust(p_value, method = "fdr")) %>% 
        #ungroup() %>%
        dplyr::arrange(p_adjust)
        
    if (disease=='PD'){
        final_results = results
    } else{
        final_results = rbind(final_results, results)
    }

}    

#final_results <- final_results %>% 
#        group_by(assay) %>% 
#        mutate(p_adjust = p.adjust(p_value, method = "fdr")) %>% 
        #ungroup() %>%
#        dplyr::arrange(p_adjust)

final_results %>% 
  vroom::vroom_write(x=., file="01-enrich_disease_lists_fdr_all.csv", delim = ",")