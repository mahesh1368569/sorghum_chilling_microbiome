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

library(DESeq2)
library(tidyverse)
library(ggplot2)
library(ggrepel)
library(pheatmap)
library(broom)

BiocManager::install("DESeq2")

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

dat <- s16

otu <- as.data.frame(dat$otu_table)
meta <- as.data.frame(dat$sample_table)
tax <- as.data.frame(dat$tax_table)

tax$ASV <- rownames(tax)

meta$Site <- factor(meta$Site, levels = c("Kansas", "South Dakota"))
meta$Planting <- factor(meta$Planting, levels = c("Regular", "Early"))
meta$Line <- factor(meta$Line)
meta$Rep <- factor(meta$Rep)

#The important point here is that Regular is the reference, therefore:
#positive log2FoldChange = enriched in Early planting

#######################################################################################################################
# STEP 1 — Find Early-vs-Regular ASVs separately in Kansas and South Dakota
#######################################################################################################################

ks_samples <- rownames(meta)[meta$Site == "Kansas"]
otu_ks <- round(as.matrix(otu[, ks_samples, drop = FALSE]))
meta_ks <- droplevels(meta[ks_samples, , drop = FALSE])

# Before DESeq2, remove extremely rare ASVs.
#Here I would require the ASV to occur in at least 10% of Kansas samples and have at least 10 reads overall:

keep_ks <- rowSums(otu_ks > 0) >= ceiling(ncol(otu_ks) * 0.10) & rowSums(otu_ks) >= 10
otu_ks_filt <- otu_ks[keep_ks, ]

# Model planting while controlling for accession and field replicate:
dds_ks <- DESeqDataSetFromMatrix(countData = otu_ks_filt, colData = meta_ks, design = ~ Rep + Line + Planting)
dds_ks <- DESeq(dds_ks, sfType = "poscounts")

res_ks <- results(dds_ks, contrast = c("Planting", "Early", "Regular"), alpha = 0.05)
res_ks <- as.data.frame(res_ks) %>% rownames_to_column("ASV")
res_ks <- res_ks %>% select(ASV, log2FoldChange, pvalue, padj) %>% rename(KS_LFC = log2FoldChange, KS_p = pvalue, KS_padj = padj)

#South dakota
sd_samples <- rownames(meta)[meta$Site == "South Dakota"]
otu_sd <- round(as.matrix(otu[, sd_samples, drop = FALSE]))
meta_sd <- droplevels(meta[sd_samples, , drop = FALSE])
keep_sd <- rowSums(otu_sd > 0) >= ceiling(ncol(otu_sd) * 0.10) & rowSums(otu_sd) >= 10
otu_sd_filt <- otu_sd[keep_sd, ]
dds_sd <- DESeqDataSetFromMatrix(countData = otu_sd_filt, colData = meta_sd, design = ~ Rep + Line + Planting)
dds_sd <- DESeq(dds_sd, sfType = "poscounts")
res_sd <- results(dds_sd, contrast = c("Planting", "Early", "Regular"), alpha = 0.05)
res_sd <- as.data.frame(res_sd) %>% rownames_to_column("ASV")
res_sd <- res_sd %>% select(ASV, log2FoldChange, pvalue, padj) %>% rename(SD_LFC = log2FoldChange, SD_p = pvalue, SD_padj = padj)

# STEP 2 — Find ASVs that respond in the same direction at both sites

cross_site <- full_join(res_ks, res_sd, by = "ASV")

# Add taxonomy
cross_site <- left_join(cross_site, tax, by = "ASV")

same_direction <- cross_site %>% filter(!is.na(KS_LFC), !is.na(SD_LFC), sign(KS_LFC) == sign(SD_LFC))

early_same <- cross_site %>% filter(KS_LFC > 0, SD_LFC > 0)

early_strict <- cross_site %>% filter(KS_LFC > 0, SD_LFC > 0, KS_padj < 0.05, SD_padj < 0.05)

nrow(early_strict)

early_strict %>% select(ASV, Phylum, Family, Genus, Species, KS_LFC, KS_padj, SD_LFC, SD_padj)

cross_site$Category <- "Other"

cross_site$Category[cross_site$KS_LFC > 0 & cross_site$SD_LFC > 0] <- "Early in both"

cross_site$Category[cross_site$KS_LFC < 0 & cross_site$SD_LFC < 0] <- "Regular in both"

cross_site$Category[cross_site$KS_LFC > 0 & cross_site$SD_LFC > 0 & cross_site$KS_padj < 0.05 & cross_site$SD_padj < 0.05] <- "Conserved Early"

cross_site %>% filter(KS_LFC > 0, SD_LFC > 0) %>% select(ASV, Phylum, Family, Genus, Species, KS_LFC, KS_padj, SD_LFC, SD_padj)

p_cross <- ggplot(cross_site, aes(x = KS_LFC, y = SD_LFC, color = Category)) + geom_hline(yintercept = 0, linetype = 2) + geom_vline(xintercept = 0, linetype = 2) + geom_point(size = 3, alpha = 0.8) + 
  theme_classic(base_size = 14) + 
  labs(x = "Kansas log2FC (Early vs Regular)", 
       y = "South Dakota log2FC (Early vs Regular)", color = NULL) + 
  theme(panel.border = element_rect(colour = "black", fill = NA))

p_cross

p_cross_l = p_cross + geom_text_repel(aes(label = ifelse(is.na(Genus) | Genus == "", ASV, 
                                             paste0(Genus, "\n", ASV))), size = 4, 
                          max.overlaps = Inf, box.padding = 0.5, point.padding = 0.3, 
                          min.segment.length = 0, show.legend = FALSE)

ggsave("cross_asv.pdf", p_cross_l, w = 7, height = 6, dpi = 1000)

cross_site %>% filter(KS_LFC > 0, SD_LFC > 0) %>% select(ASV, Phylum, Class, Order, Family, Genus, Species, KS_LFC, KS_padj, SD_LFC, SD_padj)

prev_table <- data.frame(ASV = rownames(otu))

ks_early_samples <- rownames(meta)[meta$Site == "Kansas" & meta$Planting == "Early"]
ks_reg_samples <- rownames(meta)[meta$Site == "Kansas" & meta$Planting == "Regular"]
sd_early_samples <- rownames(meta)[meta$Site == "South Dakota" & meta$Planting == "Early"]
sd_reg_samples <- rownames(meta)[meta$Site == "South Dakota" & meta$Planting == "Regular"]

prev_table$KS_Early_Prev <- rowMeans(otu[, ks_early_samples, drop = FALSE] > 0)
prev_table$KS_Regular_Prev <- rowMeans(otu[, ks_reg_samples, drop = FALSE] > 0)
prev_table$SD_Early_Prev <- rowMeans(otu[, sd_early_samples, drop = FALSE] > 0)
prev_table$SD_Regular_Prev <- rowMeans(otu[, sd_reg_samples, drop = FALSE] > 0)

prev_table <- prev_table %>% mutate(KS_Prev_Diff = KS_Early_Prev - KS_Regular_Prev, SD_Prev_Diff = SD_Early_Prev - SD_Regular_Prev)

prev_table %>% filter(ASV %in% c("ASV_88", "ASV_111", "ASV_33"))

prev_plot <- prev_table %>% filter(ASV %in% c("ASV_33", "ASV_88", "ASV_111")) %>% left_join(tax %>% select(ASV, Genus), by = "ASV") %>% mutate(Label = ifelse(is.na(Genus) | Genus == "", ASV, paste0(Genus, " (", ASV, ")"))) %>% select(ASV, Label, KS_Early_Prev, KS_Regular_Prev, SD_Early_Prev, SD_Regular_Prev) %>% pivot_longer(cols = ends_with("_Prev"), names_to = c("Site", "Planting"), names_pattern = "(KS|SD)_(Early|Regular)_Prev", values_to = "Prevalence") %>% mutate(Site = recode(Site, "KS" = "Kansas", "SD" = "South Dakota"), Planting = factor(Planting, levels = c("Early", "Regular")))

p_prev <- ggplot(prev_plot, aes(x = Site, y = Prevalence, fill = Planting)) + 
  geom_col(position = position_dodge(width = 0.75), width = 0.65, color = "black", linewidth = 0.3) + 
  geom_text(aes(label = percent(Prevalence, accuracy = 1)), 
            position = position_dodge(width = 0.75), vjust = -0.35, size = 3.5) + 
  facet_wrap(~Label, nrow = 1) + 
  scale_y_continuous(labels = percent_format(accuracy = 1), 
                     limits = c(0, 0.65), expand = expansion(mult = c(0, 0.05))) + 
  scale_fill_manual(values = c("Early" = "#E76F51", "Regular" = "#4C78A8")) + 
  labs(x = NULL, y = "ASV prevalence", fill = "Planting") + theme_classic(base_size = 14) + 
  theme(strip.background = element_blank(), strip.text = element_text(face = "bold"), 
        legend.position = "top", panel.border = element_rect(color = "black", fill = NA))

ggsave("asv_prevelenace.pdf", p_prev, w = 8, height = 7, dpi = 1000)

prev_diff_plot <- prev_table %>% filter(ASV %in% c("ASV_33", "ASV_88", "ASV_111")) %>% left_join(tax %>% select(ASV, Genus), by = "ASV") %>% mutate(Label = ifelse(is.na(Genus) | Genus == "", ASV, paste0(Genus, " (", ASV, ")"))) %>% select(Label, KS_Prev_Diff, SD_Prev_Diff) %>% pivot_longer(cols = c(KS_Prev_Diff, SD_Prev_Diff), names_to = "Site", values_to = "Prevalence_Difference") %>% 
  mutate(Site = recode(Site, "KS_Prev_Diff" = "Kansas", "SD_Prev_Diff" = "South Dakota"))

p_prev_diff <- ggplot(prev_diff_plot, aes(x = Label, y = Prevalence_Difference, color = Site, shape = Site)) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = 0, ymax = Inf, fill = "#E8F5E9", alpha = 0.35) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = 0, fill = "#FDECEC", alpha = 0.35) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.8, color = "grey30") +
  geom_point(size = 5, stroke = 1.2, position = position_dodge(width = 0.45)) +
  geom_text(aes(label = percent(Prevalence_Difference, accuracy = 1)),
            position = position_dodge(width = 0.45),
            vjust = -1.1, size = 4, fontface = "bold",
            show.legend = FALSE) +
  scale_color_manual(values = c("Kansas" = "#E76F51", "South Dakota" = "#3A86FF")) +
  scale_shape_manual(values = c("Kansas" = 16, "South Dakota" = 17)) +
  scale_y_continuous(labels = percent_format(accuracy = 1),
                     limits = c(-0.05, 0.60),
                     breaks = seq(0, 0.6, 0.1),
                     expand = expansion(mult = c(0.02, 0.08))) +
  labs(
    x = NULL,
    y = "Change in prevalence (Early − Regular)",
    color = "Site",
    shape = "Site"
  ) +
  theme_classic(base_size = 15) +
  theme(
    axis.text.x = element_text(angle = 20, hjust = 1, size = 12, face = "bold"),
    axis.text.y = element_text(size = 11),
    axis.title.y = element_text(size = 13, face = "bold"),
    legend.position = "top",
    legend.title = element_text(face = "bold"),
    strip.background = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    plot.margin = margin(10, 20, 10, 10)
  )


ggsave("prev_diff.pdf", p_prev_diff, w = 7, height = 6, dpi = 1000)



