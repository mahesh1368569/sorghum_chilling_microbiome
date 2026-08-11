library(phyloseq)
library(biomformat)
library(file2meco)
library(microeco)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(vegan)
library(patchwork)

its_biom <- import_biom("ASV_ITS_filt_taxonomy.biom")
metadata <- import_qiime_sample_data("metadata.txt")

sample_names(its_biom)
sample_names(metadata)

sorghum_its_biom <- merge_phyloseq(its_biom, metadata)

colnames(tax_table(sorghum_its_biom)) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

its <- phyloseq2meco(sorghum_its_biom)

its$tidy_dataset()

its

its$sample_table$Site <- factor(its$sample_table$Site, levels = c("Kansas", "South Dakota"))
its$sample_table$Planting <- factor(its$sample_table$Planting, levels = c("Early", "Regular"))
its$sample_table$PXS <- factor(its$sample_table$PXS, levels = c("Early Kansas", "Regular Kansas", 
                                                                "Early South D", "Regular South D"))

table(its$sample_table$Site, its$sample_table$Planting)
head(its$sample_table)

its$cal_abund()
its$cal_alphadiv()
its$cal_betadiv()

####################################################################################
##### Alpha ##################################################################
####################################################################################

alpha_its <- trans_alpha$new(dataset = its, group = "Planting", by_group = "Site")
alpha_its$cal_diff(method = "wilcox")
alpha_its$res_diff

p_its_shannon <- alpha_its$plot_alpha(measure = "Shannon", add = "jitter", jitter_shape = 21, 
                                      color_values = c("#1F77B4","#D62728"),  order_x_mean = FALSE) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12, color = "black"),
        axis.text.y = element_text(size = 12, color = "black"),
        axis.title = element_text(size = 13, face = "bold"),
        strip.text = element_text(size = 13, face = "bold"),
        strip.background = element_rect(fill = "grey90", color = "black"),
        legend.position = "right",
        legend.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 11),
        panel.grid = element_line(color = "grey", linewidth = 0.20, linetype = 2))

p_its_shannon

ggsave("ITS_Shannon_diversity.pdf", p_its_shannon, width = 6, height = 5, dpi = 1000)

p_its_chao <- alpha_its$plot_alpha(measure = "Chao1", add = "jitter", 
                                   color_values = c("#1F77B4","#D62728"), jitter_shape = 21) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12, color = "black"),
        axis.text.y = element_text(size = 12, color = "black"),
        axis.title = element_text(size = 13, face = "bold"),
        strip.text = element_text(size = 13, face = "bold"),
        strip.background = element_rect(fill = "grey90", color = "black"),
        legend.position = "right",
        legend.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 11),
        panel.grid = element_line(color = "grey", linewidth = 0.20, linetype = 2))

p_its_chao

ggsave("ITS_Chao1_diversity.pdf", p_its_chao, width = 6, height = 5, dpi = 1000)

####################################################################################
##### Beta ##################################################################
####################################################################################

beta_its <- trans_beta$new(dataset = its, group = "PXS", measure = "bray")

beta_its$cal_ordination(method = "PCoA")

p_its_beta <- beta_its$plot_ordination(plot_color = "PXS", plot_shape = "Planting", plot_type = c("point", "ellipse"), 
                                       point_size = 3, point_alpha = 0.8) + 
  scale_color_manual(values = c("Early Kansas" = "#1F77B4", "Regular Kansas" = "#D62728", 
                                "Early South D" = "#5D478B", "Regular South D" = "#CD6600" )) + 
  theme_minimal() + 
  theme(axis.text.x = element_text(size = 14), 
        axis.text.y = element_text(size = 14), 
        axis.title.x = element_text(size = 14), 
        axis.title.y = element_text(size = 14), 
        legend.position = "right", legend.title = element_text(size = 14), 
        legend.text = element_text(size = 14), 
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 1))

p_its_beta

ggsave("ITS_Bray_PCoA.pdf", p_its_beta, width = 7, height = 6, dpi = 1000)
