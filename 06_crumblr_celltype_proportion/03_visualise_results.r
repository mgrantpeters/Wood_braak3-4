library("ggplot2")
library("grid")
library("dplyr")
library("palr")
library("ggh4x")

#Double check directionality of + and -

#setwd("/home/MRMelissaGrantPeters/MinaRyten/Melissa/Abundance_final_data/")
crumblr = read.csv("/home/MRMelissaGrantPeters/MinaRyten/Melissa/teamWood/Abundance_final_data/processed_data/results_for_plotting/results_crumblr_All_Abundance_raw_.csv")
l1_colour_map = data.frame(
  assay_l1 = c("Neurons","OG","Immune" ,"Astro","Endomural", "All"), 
  colour = c("#024a68", "#5f356c","#cf2fb3", "#cf8219", "#bcbebc", "#808080")
)
l1_order = c("Neurons","OG","Immune" ,"Astro","Endomural", "All")
l1_labels = c("Neurons","OG","Immune" ,"Astro","Endomural", "All")

tissue_order = c("C", "F", "P", "all_cortex")
tissue_labels = c("Cingulate", "Frontal", "Parietal", "All")
crumblr = crumblr %>% dplyr::mutate(major_cell_type = factor(major_cell_type, levels=l1_order, labels=l1_labels)) %>% 
  dplyr::mutate(brain_region = factor(brain_region, levels=tissue_order, labels=tissue_labels))
strip_text_size = 7

crumblr$Direction <- ifelse(crumblr$logFC > 0, "Up in PD", "Down in PD")

for (n in 1:4){
    crumblr_res = crumblr[crumblr$annotation_level==as.character(n),]
    ggplot(data=crumblr_res, aes(x=cell_type, y=-log10(adj.P.Val))) + 
    geom_point(aes(alpha=adj.P.Val<0.05), size=4) + 
    #scale_color_identity() +
    scale_alpha_discrete(range=c(0.4,1)) +
    #geom_point(shape=1, size=4, colour="black") + 
    geom_point(aes(shape = Direction), colour = "white") +
    scale_shape_manual(values = c(25, 24)) +
    geom_hline(yintercept = -log10(0.05), alpha=0.5, linetype=2) +
    facet_grid2(brain_region ~ major_cell_type, scales="free_x", space="free_x", switch="y", strip=strip_themed(background_x = elem_list_rect(fill = l1_colour_map$colour)))+
    scale_y_continuous(position = "right") +
    theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 90, hjust = 1),
    axis.text.y = element_text(angle=0, size=strip_text_size),
    strip.text.y.left = element_text(angle=0, size=strip_text_size),
    strip.background.y = element_blank(),
    strip.text.x.top = element_text(angle=90, size=strip_text_size, colour="white", face="bold"),
    legend.position = "bottom",
    legend.text = element_text(size=strip_text_size), 
    axis.title = element_text(size=strip_text_size))
    ggsave(filename=paste0("/home/MinaRyten/Melissa/teamWood/Abundance_final_data/plots/crumblr_results/PDctrl_differential_abundance_annotation_level_", n,".png"),
       plot = last_plot(),
       device ="png",
       scale = 1,
       width = 16,
       height = 7,
       units = c("in"),
       dpi = 300 )
}


crumblr = read.csv("/home/MRMelissaGrantPeters/MinaRyten/Melissa/teamWood/Abundance_final_data/processed_data/results_for_plotting/results_crumblr_All_Abundance_raw_.csv")

crumblr_res = crumblr[crumblr$tier=='2',]

strip_text_size = 7
ggplot(data=crumblr_res, aes(x=CellType, y=-log10(adj.P.Val))) + 
      geom_point(aes(colour=Colour, alpha=adj.P.Val<0.05), size=4) + 
      scale_color_identity() +
      scale_alpha_discrete(range=c(0.4,1)) +
      #geom_point(shape=1, size=4, colour="black") + 
      geom_hline(yintercept = -log10(0.05), alpha=0.5, linetype=2) +
      #ggh4x::facet_grid2(brain_region~assay_l1, scales="free_x", space="free_x", switch="y", strip=strip_themed(background_x = elem_list_rect(fill = l1_colour_map$colour))) +
      facet_grid(rows = vars(brain_region))+
      #scale_alpha_identity() +
      scale_y_continuous(position = "right") +
      theme(
        panel.grid = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        #axis.text.x = element_blank(),
        axis.text.y = element_text(angle=0, size=strip_text_size),
        strip.text.y.left = element_text(angle=0, size=strip_text_size),
        strip.background.y = element_blank(),
        strip.text.x.top = element_text(angle=45, size=strip_text_size, colour="white", face="bold"),
        legend.position = "right",
        legend.text = element_text(size=strip_text_size), 
        axis.title = element_text(size=strip_text_size)
      )

ggsave(filename="plots/plot_results.png",
  plot = last_plot(),
  device ="png",
  scale = 1,
  width = 12,
  height = 9,
  units = c("in"),
  dpi = 300 
)