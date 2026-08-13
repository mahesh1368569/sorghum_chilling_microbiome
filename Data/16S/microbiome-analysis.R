# ---- Core data handling & manipulation ----
library(dplyr)
library(tidyr)
library(tidyverse)
library(plyr)
library(magrittr)
library(parallel)
library(reshape2)

# ---- Microbiome data processing ----
library(phyloseq)
library(biomformat)
library(file2meco)
library(microeco)
library(MicrobiomeStat)
library(meconetcomp)
library(WGCNA)
library(ggClusterNet)
library(ape)
library(picante)
library(Biostrings)

# ---- Differential abundance & compositional analysis ----
library(metagenomeSeq)
library(ALDEx2)
library(ANCOMBC)

# ---- Visualization ----
library(ggplot2)
library(ggpubr)
library(ggtree)
library(tidygraph)
library(paletteer)
library(colorspace)
library(ComplexHeatmap)
library(circlize)
library(vegan)
library(ggraph)

# ---- Statistical modeling ----
library(lme4)
library(lmerTest)
library(multcomp)
library(emmeans)
library(multcompView)
library(dplyr)
library(usethis)
library(nlMS)
library(iCAMP)
library(minpack.lm)
library(Hmisc)
library(Biostrings)
library(meconetcomp)

# ---- Package manager ----
library(BiocManager)

##### Importing files ########

s16_biom <- import_biom("asv_table_16s_filt_taxonomy.biom")
metadata <- import_qiime_sample_data("metadata.txt")

sample_names(s16_biom)
sample_names(metadata)

sorghum_16s_biom <- merge_phyloseq(s16_biom, metadata)

colnames(tax_table(sorghum_16s_biom)) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

s16 <- phyloseq2meco(sorghum_16s_biom)

s16$tidy_dataset()

s16

s16$sample_table$Site <- factor(s16$sample_table$Site, levels = c("Kansas", "South Dakota"))
s16$sample_table$Planting <- factor(s16$sample_table$Planting, levels = c("Early", "Regular"))
s16$sample_table$PXS <- factor(s16$sample_table$PXS, levels = c("Early Kansas", "Regular Kansas", 
                                                                "Early South D", "Regular South D"))

table(s16$sample_table$Site, s16$sample_table$Planting)
head(s16$sample_table)

s16$cal_abund()
s16$cal_alphadiv()
s16$cal_betadiv()

####################################################################################
##### Alpha ##################################################################
####################################################################################

alpha_16s <- trans_alpha$new(dataset = s16, group = "Planting", by_group = "Site")
alpha_16s$cal_diff(method = "wilcox")
alpha_16s$res_diff

p_16s_shannon <- alpha_16s$plot_alpha(measure = "Shannon", add = "jitter", jitter_shape = 21, 
                                      color_values = c("#3CB371","#7B68EE"),  order_x_mean = FALSE) +
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

p_16s_shannon

ggsave("16S_Shannon_diversity.pdf", p_16s_shannon, width = 6, height = 5, dpi = 1000)

p_16s_chao <- alpha_16s$plot_alpha(measure = "Chao1", add = "jitter", 
                                   color_values = c("#3CB371","#7B68EE"), jitter_shape = 21) +
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

p_16s_chao

ggsave("16S_Chao1_diversity.pdf", p_16s_chao, width = 6, height = 5, dpi = 1000)

####################################################################################
##### Beta ##################################################################
####################################################################################

beta_16s <- trans_beta$new(dataset = s16, group = "PXS", measure = "bray")

beta_16s$cal_ordination(method = "PCoA")

p_16s_beta <- beta_16s$plot_ordination(plot_color = "PXS", plot_shape = "Planting", plot_type = c("point", "ellipse"), 
                                       point_size = 3, point_alpha = 0.8) + 
  scale_color_manual(values = c("Early Kansas" = "#3CB371", "Regular Kansas" = "#D62728", 
                                "Early South D" = "#7B68EE", "Regular South D" = "#CD6600" )) + 
  theme_minimal() + 
  theme(axis.text.x = element_text(size = 14), 
        axis.text.y = element_text(size = 14), 
        axis.title.x = element_text(size = 14), 
        axis.title.y = element_text(size = 14), 
        legend.position = "right", legend.title = element_text(size = 14), 
        legend.text = element_text(size = 14), 
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 1))

p_16s_beta

ggsave("16S_Bray_PCoA.pdf", p_16s_beta, width = 7, height = 6, dpi = 1000)

