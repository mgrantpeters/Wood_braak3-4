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
  celltypes_l2 = c(
     "Astro" = "#00A88C",  # Replace "Value1" with the actual value and "#FF5733" with the desired hex color
     "Endomural" = "#FFD347",
     "Immune" = "#91D1C2",
     "EN" = "#354977",
     "IN" = "#4DBBD5",
     "OG" = "#E64B35",
     "OPC" = "#F39B7F"
    # Add more colors as needed
  )
)

# read in cell type mapping 
celltype_level_mapping = vroom::vroom("~/MinaRyten/Aine/wood_full/data/annotation_mapping_tables/annotation3_to_others.csv", show_col_types = FALSE) %>% 
  dplyr::filter(annot_level==2) %>% 
  dplyr::rename(assay=annotation_level_mapped, l2=assay) %>% 
  dplyr::select(-annot_level) %>% 
  dplyr::distinct()

# read in cell type mapping 
#celltype_level_mapping = vroom::vroom("~/MinaRyten/Aine/wood_full/data/annotation_mapping_tables/annotation3_to_others.csv", show_col_types = FALSE) %>% 
#  dplyr::filter(annot_level==2) %>% 
#  dplyr::rename(assay=annotation_level_mapped, l2=assay) %>% 
#  dplyr::select(-annot_level) %>% 
#  dplyr::distinct()

# read in deg disease enrichment results
deg_disease_enrichment = vroom::vroom("~/MinaRyten/Aine/wood_full/data/enrichment_overlap_disease_lists/08-enrich_disease_lists.csv", show_col_types = FALSE) %>% 
  dplyr::filter(annot_level==3)

# read in degs 
degs = readRDS("~/MinaRyten/Aine/wood_full/data/dreamlet/de_X_sum_all_DE.rds") %>% 
  janitor::clean_names() %>% 
  dplyr::mutate(direction = case_when(
    log_fc<0 ~ "Down-regulated",
    log_fc>0 ~ "Up-regulated"
  )) %>% 
  tidyr::separate(coef, into=c("tissue", "comparison"))

degs %>% 
  dplyr::filter(annot_level=="3") %>% 
  dplyr::mutate(tissue = str_replace(tissue, "C", "Cingulate"),
          tissue = str_replace(tissue, "F", "Frontal"),
          tissue = str_replace(tissue, "P", "Parietal")) %>%
  dplyr::mutate(assay = str_replace(assay, "EN", "Excitatory"),
          assay = str_replace(assay, "IN", "Inhibitory")) %>%
  dplyr::left_join(
    x=., 
    y=celltype_level_mapping, 
    by="assay") %>% 
  dplyr::group_by(l2, assay, tissue, comparison, direction) %>% 
  dplyr::summarise(
    n_sig = sum(adj_p_val<0.05),
    all = dplyr::n()
  ) %>% 
  dplyr::mutate(dir_n_sig = case_when(
    is.na(n_sig) ~ 0,
    (direction=="Down-regulated") & (n_sig != 0) ~ -1*n_sig,
    (direction=="Up-regulated") & (n_sig != 0) ~ 1*n_sig,
    TRUE ~ n_sig
  )) %>% 
  {
    ggplot(., aes(x=dir_n_sig, y= assay, fill=l2)) + 
      geom_col(aes(alpha=direction)) + 
      facet_grid(l2~tissue, space="free_y", scales="free_y", switch="y") +
      scale_alpha_manual(values=c("Down-regulated"=0.4, "Up-regulated"=1)) + 
      guides(fill="none") + 
      labs(x="Number of sig. DEGs", y="", fill="Broad cell type", alpha="Log2 fold-change\n direction") + 
      theme(
        legend.position = "right",
        strip.text.y.left = element_text(angle=0)
      ) +
      scale_fill_manual(values = palettes$celltypes_l2) + 
      scale_x_continuous(breaks = pretty(.$dir_n_sig), labels = abs(pretty(.$dir_n_sig)))
  } %>% ggsave(
    "plots/DEGs_celltype.pdf",
    plot = .,
    device = "pdf",
    scale = 1,
    width = 7,
    height = 5,
    units = c("in"),
    dpi = 300, 
    bg = "white"
  )