library(vroom)
library(dplyr)
library(ggplot2)
library(ggh4x)
library(stringr)

# import Regina's ggplot theme and set as default
source("~/MinaRyten/Aine/wood_full/scripts/functions/get_rhr_theme.R")
theme_rhr = get.rhr.theme(text_size = 10)

# import palettes 
source("/home/MinaRyten/Aine/hardy_snrnaseq/manuscript/manuscript_palettes.R")
palettes = manuscript_palettes()

# Define your color palette
palettes <- list(
  celltypes_l1 = c(
     "Astro" = "#00A88C",  # Replace "Value1" with the actual value and "#FF5733" with the desired hex color
     "Endomural" = "#FFD347",
     "Immune" = "#91D1C2",
     "Neurons" = "#354977",
     "OG" = "#E64B35"
    # Add more colors as needed
  )
)

# read in cell type mapping 
celltype_level_mapping = vroom::vroom("~/MinaRyten/Aine/wood_full/data/annotation_mapping_tables/annotation3_to_others.csv", show_col_types = FALSE) %>% 
  dplyr::filter(annot_level==1) %>% 
  dplyr::rename(assay=annotation_level_mapped, l1=assay) %>% 
  dplyr::select(-annot_level) %>% 
  dplyr::distinct()

# read in deg disease enrichment results
deg_disease_enrichment = vroom::vroom("~/MinaRyten/Melissa/teamWood/check_ADPD_enrichment/01-enrich_disease_lists_fdr_all.csv", show_col_types = FALSE) %>% 
  dplyr::filter(annot_level==3)

deg_disease_enrich = vroom::vroom("~/MinaRyten/Melissa/teamWood/check_ADPD_enrichment/processed_data/overlap_disease_lists.csv", show_col_types = FALSE) %>% 
  dplyr::filter(annot_level==3) %>% 
  dplyr::mutate(disease_list = case_when(  
    disease_list == "ad_common_gwas_bellenguez" ~ "AD",
    #disease_list == "pd_common_gwas_kim" ~ "PD_kim",
    #disease_list == "pd_common_gwas_kim" ~ "PD_kim",
    #disease_list == "pd_common_gwas_nalls" ~ "PD_nalls",
    disease_list == "pd_common_gwas_nallskim" ~ "PD",
    TRUE ~ assay
  )) %>%
  dplyr::left_join(
    x=., 
    y=deg_disease_enrichment,
    by=c("annot_level", "comparison", "tissue", "assay", c("disease_list" = "disease_enrichment"))
  ) %>% 
  dplyr::group_by(annot_level, comparison, tissue, assay) %>% 
  #dplyr::mutate(p_adjust = p.adjust(pvalue, "fdr")) %>% 
  dplyr::ungroup() %>% 
  dplyr::mutate(sig_indicator = ifelse(p_value<0.05, "*", "")) %>% 
  tidyr::complete(annot_level, comparison, tissue, assay, disease_list, fill = list(n_overlap=0, n_overlap_sig=0, perc_disease_gene=0)) %>%
  dplyr::left_join(
    x=., 
    y=celltype_level_mapping, 
    by="assay") %>% 
  # dplyr::mutate(disease_list = case_when(
  #   disease_list=="ad_common_gwas_bellenguez" ~ "AD genes",
  #   disease_list=="pd_common_gwas_nallskim" ~ "PD genes",
  #   TRUE ~ disease_list
  # )) %>% 
  #dplyr::mutate(
  #  tissue_comparison = paste0(tissue, "\n", comparison)
  #) %>% 
  dplyr::mutate(alpha = (perc_disease_gene - min(perc_disease_gene)) / (max(perc_disease_gene) - min(perc_disease_gene))) %>% 
  dplyr::mutate(tissue = str_replace(tissue, "C", "Cingulate"),
          tissue = str_replace(tissue, "F", "Frontal"),
          tissue = str_replace(tissue, "P", "Parietal")) %>% 
  {
    ggplot(data=., 
           aes(y=assay, 
               x=reorder(disease_list, perc_disease_gene), 
               fill=l1, 
               alpha=alpha, 
               label=sig_indicator)) + 
      geom_tile() + 
      facet_grid(l1~tissue, space = "free_y", scales="free_y", switch="y") + 
      theme(
        axis.text.x = element_text(angle=45, hjust=1, vjust=1),
        strip.text.y.left = element_text(angle=0),
        axis.line = element_blank(),
        legend.position = "right"
      ) + 
      labs(x="", y="", alpha="% disease genes\nin P-adj<0.05\nDEGs", fill="Broad cell type") + 
      # scale_fill_manual(values=palettes$celltypes_l1) + 
      scale_fill_manual(values = palettes$celltypes_l1) + 
      scale_alpha_identity(guide = "legend") +
      guides(label="none", colour="none", text="none", fill="none")
  }

deg_disease_enrich = deg_disease_enrich + 
  (deg_disease_enrich$data %>% 
     dplyr::filter(p_adjust<0.05) %>% 
     geom_tile(data=., aes(y=assay, x=reorder(disease_list, perc_disease_gene)), alpha=0, colour="black", linewidth=1))
 
vroom::vroom_write(x=deg_disease_enrich$data, file="processed_data/deg_disease_enrich_results.csv")

ggsave(
    "plots/enrichment_PD_genes.pdf",
    plot = deg_disease_enrich,
    device = "pdf",
    scale = 1,
    width = 7,
    height = 5,
    units = c("in"),
    dpi = 300, 
    bg = "white"
  )
