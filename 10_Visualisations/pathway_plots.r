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
celltype_level_mapping = vroom::vroom("~/MinaRyten/Aine/wood_full/data/annotation_mapping_tables/annotation3_to_others.csv", show_col_types = FALSE) %>% 
  dplyr::filter(annot_level==2) %>% 
  dplyr::rename(assay=annotation_level_mapped, l2=assay) %>% 
  dplyr::select(-annot_level) %>% 
  dplyr::distinct()

# Perform GO reduce on pathways
source("~/MinaRyten/Aine/hardy_snrnaseq/scripts/functions/03fa-go_reduce.R")

go_reduce_wrapper = function(res){
  reduced_terms = res %>%
    dplyr::filter(source=="GO", p_adjust<0.05) %>%
    dplyr::select(id, ontology) %>%
    dplyr::rename(go_type=ontology, go_id=id) %>%
    dplyr::distinct() %>%
    updated_go_reduce(., threshold = 0.85)
  
  res_go_collapsed = res %>%
    dplyr::mutate(gene_ratio = num/denom) %>% 
    dplyr::filter(source=="GO", p_adjust<0.05) %>%
    dplyr::left_join(
      x=.,
      y=reduced_terms %>% dplyr::select(-parent_sim_score),
      by=c("id"="go_id", "ontology"="go_type")
    ) %>%
    dplyr::group_by(annot_level, assay, tissue, comparison, direction, ontology, parent_term) %>%
    dplyr::summarise(
      p_adjust = min(p_adjust),
      n_collapsed = n(),
      collapse_genes = paste0(gene_id, collapse="/"),
      gene_ratio = mean(gene_ratio)
    ) %>%
    dplyr::mutate(description = paste0(parent_term, " (", n_collapsed, ")")) %>% 
    dplyr::ungroup() %>% 
    dplyr::rowwise() %>% 
    dplyr::mutate(collapse_genes = paste0(unique(stringr::str_split(collapse_genes, "\\/")[[1]]), collapse="/")) %>% 
    dplyr::mutate(num_genes = length(stringr::str_split(collapse_genes, "\\/")[[1]])) %>% 
    dplyr::ungroup()
}
res_padjfilt = vroom::vroom("~/MinaRyten/Aine/wood_full/data/clusterProfiler/07-clusterProfiler_res.csv", show_col_types = FALSE) %>% 
  janitor::clean_names() %>% 
  dplyr::arrange(qvalue) %>% 
  tidyr::separate(gene_ratio, into=c("num", "denom"), sep="\\/", remove=T, convert=T) %>% 
  dplyr::mutate(gene_ratio = num/denom) %>% 
  dplyr::mutate(ontology = case_when(
    source=="KEGG" ~ "KEGG",
    source=="Reactome" ~ "Reactome",
    TRUE ~ ontology
  )) %>% 
  dplyr::filter(num>1 & denom>1)

res_padjfilt_collapsed = go_reduce_wrapper(res_padjfilt)

# Filter and prepare data
plot_data = res_padjfilt_collapsed %>%
  dplyr::ungroup() %>%
  dplyr::filter(annot_level == 3) %>%
  dplyr::left_join(x = ., y = celltype_level_mapping, by = "assay") %>%
  dplyr::mutate(p_adjust = as.numeric(p_adjust)) %>% 
  dplyr::group_by(parent_term) %>% 
  dplyr::mutate(
    description = paste0(parent_term, "(", max(n_collapsed), ")"),
    collapse_genes = paste0(collapse_genes, collapse="/")
  ) %>% 
  dplyr::ungroup() %>% 
  dplyr::select(-annot_level, -n_collapsed) %>% 
  dplyr::rowwise() %>% 
  dplyr::mutate(collapse_genes = paste0(unique(stringr::str_split(collapse_genes, "\\/")[[1]]), collapse="/")) %>% 
  dplyr::mutate(num_genes = length(stringr::str_split(collapse_genes, "\\/")[[1]])) %>% 
  dplyr::ungroup()

assay_order = plot_data %>%
  dplyr::distinct(assay, l2) %>%
  dplyr::arrange(l2)  # Order assays by l2

plot_data = plot_data %>%
  dplyr::mutate(assay = factor(assay, levels = assay_order$assay))  # Apply the ordered levels

# add pathway ordering 
plot_data = plot_data %>% 
  dplyr::left_join(
    x=., 
    y=plot_data %>% 
      dplyr::select(assay, tissue, comparison, description, p_adjust) %>% 
      tidyr::pivot_wider(names_from="tissue", values_from="p_adjust") %>% 
      dplyr::mutate(
        `F` = ifelse(is.na(`F`), 1, `F`),
        `P` = ifelse(is.na(`P`), 1, `P`),
        `C` = ifelse(is.na(`C`), 1, `C`)
      ) %>% 
      dplyr::mutate(pathway_order = case_when(
        `F`<0.05 & `C`<0.05 & `P`>0.05 ~ 5, # sig in 2/3 regions
        `F`<0.05 & `P`<0.05 & `C`>0.05 ~ 5, # sig in 2/3 regions
        `P`<0.05 & `C`<0.05 & `F`>0.05 ~ 5, # sig in 2/3 regions
        `F`<0.05 & `C`<0.05 & `P`<0.05 ~ 4, # all sig
        `F`<0.05 & `C`>0.05 & `P`>0.05 ~ 3, # only F 
        `C`<0.05 & `F`>0.05 & `P`>0.05 ~ 2, # only C
        `P`<0.05 & `C`>0.05 & `F`>0.05 ~ 1, # only P
        TRUE ~ 0 # other
      )) %>% 
      dplyr::mutate(pathway_order_label = case_when(
        pathway_order == 5 ~ "2/3 tissues",
        pathway_order == 4 ~ "All tissues",
        pathway_order == 3 ~ "Only Frontal",
        pathway_order == 2 ~ "Only Cingulate",
        pathway_order == 1 ~ "Only Parietal",
        TRUE ~ "Other" # other
      )) %>% 
      dplyr::group_by(description) %>% 
      dplyr::mutate(n_assays_sig_for_pw = sum(`F`<0.05) + sum(`P`<0.05) + sum(`C`<0.05)) %>% 
      dplyr::ungroup() %>% 
      dplyr::select(-`F`, -`P`, -`C`) %>% 
      dplyr::group_by(pathway_order, description) %>%
      dplyr::arrange(n_assays_sig_for_pw) %>%
      tibble::rowid_to_column("pathway_order_id") %>% 
      dplyr::ungroup(), by=c("assay", "description", "comparison")
  ) %>%  dplyr::mutate(tissue = str_replace(tissue, "C", "Cingulate"),
                       tissue = str_replace(tissue, "F", "Frontal"),
                       tissue = str_replace(tissue, "P", "Parietal"))

plot_data = plot_data[!is.na(plot_data$parent_term),] 


# Main plot without x-axis labels
main_plot = ggplot(data = plot_data, aes(
  y = reorder(description, pathway_order_label),
  x = reorder(assay, l2),
  fill = l2,             # Fill based on l2 variable
  shape = direction,      # Shape based on Direction
  size = gene_ratio
)) +
  geom_point(color="black", alpha=0.75) +  # Black border for each triangle
  labs(
    x = NULL,  # Remove x-axis label
    y = "",
    shape = "Direction",
    size = "Gene ratio"
  ) +
  scale_shape_manual(
    values = c("Up-regulated" = 24, "Down-regulated" = 25)  # Upward and downward triangles
  ) +
  facet_grid(pathway_order_label~tissue, space = "free", scales="free") +
  scale_fill_manual(values = palettes$celltypes_l2) + 
  guides(fill = "none") +  # Remove the l2 legend
  theme(
    legend.position = "right",
    legend.direction = "vertical",
    strip.text.y.right = element_text(angle=0)
  )

main_plot %>% ggsave(
    "plots/pathways_celltype.pdf",
    plot = .,
    device = "pdf",
    scale = 1,
    width = 10.3,
    height = 12,
    units = c("in"),
    #dpi = 300, 
    bg = "white"
  )