# ---- Core data handling & manipulation ----
library(dplyr)
library(tidyr)
library(tidyverse)
library(magrittr)

# ---- Microbiome data processing ----
library(phyloseq)
library(biomformat)
library(file2meco)
library(microeco)

# ---- Differential abundance ----
library(DESeq2)

# ---- Visualization ----
library(ggplot2)
library(ggrepel)
library(scales)

##### IMPORT ITS FILES #########################################################

its_biom <- import_biom("ASV_ITS_filt_taxonomy.biom")
metadata <- import_qiime_sample_data("metadata.txt")

sample_names(its_biom)
sample_names(metadata)

sorghum_its_biom <- merge_phyloseq(its_biom, metadata)

colnames(tax_table(sorghum_its_biom)) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

its <- phyloseq2meco(sorghum_its_biom)
its$tidy_dataset()
its

##### EXTRACT TABLES ###########################################################

dat <- its
otu <- as.data.frame(dat$otu_table)
meta <- as.data.frame(dat$sample_table)
tax <- as.data.frame(dat$tax_table)

tax$ASV <- rownames(tax)

meta$Site <- factor(meta$Site, levels = c("Kansas", "South Dakota"))
meta$Planting <- factor(meta$Planting, levels = c("Regular", "Early"))
meta$Line <- factor(meta$Line)
meta$Rep <- factor(meta$Rep)

# Regular is the reference level:
# positive log2FoldChange = enriched in Early planting
# negative log2FoldChange = enriched in Regular planting

#######################################################################################################################
# STEP 1 — EARLY-vs-REGULAR FUNGAL ASVs SEPARATELY IN KANSAS AND SOUTH DAKOTA
#######################################################################################################################

##### KANSAS ##################################################################

ks_samples <- rownames(meta)[meta$Site == "Kansas"]
otu_ks <- round(as.matrix(otu[, ks_samples, drop = FALSE]))
meta_ks <- droplevels(meta[ks_samples, , drop = FALSE])

# Keep fungal ASVs occurring in at least 10% of Kansas samples and with >=10 reads overall
keep_ks <- rowSums(otu_ks > 0) >= ceiling(ncol(otu_ks) * 0.10) & rowSums(otu_ks) >= 10
otu_ks_filt <- otu_ks[keep_ks, ]

dds_ks <- DESeqDataSetFromMatrix(countData = otu_ks_filt, colData = meta_ks, design = ~ Rep + Line + Planting)
dds_ks <- DESeq(dds_ks, sfType = "poscounts")

res_ks <- results(dds_ks, contrast = c("Planting", "Early", "Regular"), alpha = 0.05)
res_ks <- as.data.frame(res_ks) %>% rownames_to_column("ASV")
res_ks <- res_ks %>% select(ASV, log2FoldChange, pvalue, padj) %>% rename(KS_LFC = log2FoldChange, KS_p = pvalue, KS_padj = padj)

##### SOUTH DAKOTA ############################################################

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

#######################################################################################################################
# STEP 2 — FIND FUNGAL ASVs RESPONDING IN THE SAME DIRECTION AT BOTH SITES
#######################################################################################################################

cross_site <- full_join(res_ks, res_sd, by = "ASV")
cross_site <- left_join(cross_site, tax, by = "ASV")

same_direction <- cross_site %>% filter(!is.na(KS_LFC), !is.na(SD_LFC), sign(KS_LFC) == sign(SD_LFC))
early_same <- cross_site %>% filter(KS_LFC > 0, SD_LFC > 0)
regular_same <- cross_site %>% filter(KS_LFC < 0, SD_LFC < 0)
early_strict <- cross_site %>% filter(KS_LFC > 0, SD_LFC > 0, KS_padj < 0.05, SD_padj < 0.05)
regular_strict <- cross_site %>% filter(KS_LFC < 0, SD_LFC < 0, KS_padj < 0.05, SD_padj < 0.05)

nrow(early_same)
nrow(early_strict)
nrow(regular_same)
nrow(regular_strict)

early_strict %>% select(ASV, Phylum, Class, Order, Family, Genus, Species, KS_LFC, KS_padj, SD_LFC, SD_padj)
regular_strict %>% select(ASV, Phylum, Class, Order, Family, Genus, Species, KS_LFC, KS_padj, SD_LFC, SD_padj)

##### CLASSIFY CROSS-SITE RESPONSE ############################################

cross_site$Category <- "Other"
cross_site$Category[!is.na(cross_site$KS_padj) & !is.na(cross_site$SD_padj) & cross_site$KS_LFC > 0 & cross_site$SD_LFC > 0 & cross_site$KS_padj < 0.05 & cross_site$SD_padj < 0.05] <- "Conserved Early"
cross_site$Category[!is.na(cross_site$KS_padj) & !is.na(cross_site$SD_padj) & cross_site$KS_LFC < 0 & cross_site$SD_LFC < 0 & cross_site$KS_padj < 0.05 & cross_site$SD_padj < 0.05] <- "Conserved Regular"

cross_site$Category <- factor(cross_site$Category, levels = c("Other", "Conserved Early", "Conserved Regular"))

cross_colors_its <- c("Other" = "#BDBDBD", "Conserved Early" = "#E07A5F", "Conserved Regular" = "#3D5A80")

##### LABEL ONLY CONSERVED ASVs #####

cross_site$Label <- ifelse(cross_site$Category %in% c("Conserved Early", "Conserved Regular"), cross_site$ASV, "")
cross_site$Label <- ifelse(cross_site$Category %in% c("Conserved Early", "Conserved Regular"), ifelse(is.na(cross_site$Genus) | cross_site$Genus == "", cross_site$ASV, paste0(cross_site$Genus, "\n", cross_site$ASV)), "")

##### CROSS-SITE EFFECT-SIZE PLOT #####

p_cross_its <- ggplot(cross_site, aes(x = KS_LFC, y = SD_LFC, color = Category)) + 
  geom_hline(yintercept = 0, linetype = 2, color = "grey40") + 
  geom_vline(xintercept = 0, linetype = 2, color = "grey40") + 
  geom_point(size = 3, alpha = 0.8) + scale_color_manual(values = cross_colors_its) + 
  theme_classic(base_size = 14) + 
  labs(x = "Kansas log2FC (Early vs Regular)", 
       y = "South Dakota log2FC (Early vs Regular)", color = NULL) + 
  theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.8), 
        legend.position = "right")

p_cross_its_l <- p_cross_its + geom_text_repel(aes(label = Label), size = 3.8, max.overlaps = Inf, box.padding = 0.5, point.padding = 0.3, min.segment.length = 0, show.legend = FALSE)

p_cross_its_l

ggsave("ITS_cross_asv.pdf", p_cross_its_l, width = 7, height = 6, dpi = 1000)

#######################################################################################################################
# STEP 3 — PREVALENCE OF CONSISTENT EARLY-RESPONDING FUNGAL ASVs
#######################################################################################################################

strict_candidates <- early_strict %>% mutate(Combined_LFC = KS_LFC + SD_LFC) %>% arrange(KS_padj + SD_padj, desc(Combined_LFC))
all_early_candidates <- early_same %>% mutate(Combined_LFC = KS_LFC + SD_LFC) %>% arrange(desc(Combined_LFC))
top_its_asvs <- unique(c(strict_candidates$ASV, all_early_candidates$ASV))
top_its_asvs <- head(top_its_asvs, 4)

top_its_asvs

if (length(top_its_asvs) == 0) stop("No fungal ASVs were enriched in Early planting at both Kansas and South Dakota.")

##### CALCULATE PREVALENCE #####################################################

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

prev_table %>% filter(ASV %in% top_its_asvs)

##### PREVALENCE BAR PLOT ######################################################

prev_plot_its <- prev_table %>% filter(ASV %in% top_its_asvs) %>% left_join(tax %>% select(ASV, Genus), by = "ASV") %>% mutate(Label = ifelse(is.na(Genus) | Genus == "", ASV, paste0(Genus, " (", ASV, ")"))) %>% select(ASV, Label, KS_Early_Prev, KS_Regular_Prev, SD_Early_Prev, SD_Regular_Prev) %>% pivot_longer(cols = ends_with("_Prev"), names_to = c("Site", "Planting"), names_pattern = "(KS|SD)_(Early|Regular)_Prev", values_to = "Prevalence") %>% mutate(Site = recode(Site, "KS" = "Kansas", "SD" = "South Dakota"), Planting = factor(Planting, levels = c("Early", "Regular")))

p_prev_its <- ggplot(prev_plot_its, aes(x = Site, y = Prevalence, fill = Planting)) + geom_col(position = position_dodge(width = 0.75), width = 0.65, color = "black", linewidth = 0.3) + geom_text(aes(label = percent(Prevalence, accuracy = 1)), position = position_dodge(width = 0.75), vjust = -0.35, size = 3.5) + facet_wrap(~Label, nrow = 1) + scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, 1), expand = expansion(mult = c(0, 0.06))) + scale_fill_manual(values = c("Early" = "#4C956C", "Regular" = "#6C5B7B")) + labs(x = NULL, y = "ASV prevalence", fill = "Planting") + theme_classic(base_size = 14) + theme(strip.background = element_blank(), strip.text = element_text(face = "bold"), legend.position = "top", panel.border = element_rect(color = "black", fill = NA))

p_prev_its

ggsave("ITS_ASV_prevalence.pdf", p_prev_its, width = 10, height = 6, dpi = 1000)

##### PREVALENCE-DIFFERENCE PLOT ##############################################

prev_diff_plot_its <- prev_table %>% filter(ASV %in% top_its_asvs) %>% left_join(tax %>% select(ASV, Genus), by = "ASV") %>% mutate(Label = ifelse(is.na(Genus) | Genus == "", ASV, paste0(Genus, " (", ASV, ")"))) %>% select(Label, KS_Prev_Diff, SD_Prev_Diff) %>% pivot_longer(cols = c(KS_Prev_Diff, SD_Prev_Diff), names_to = "Site", values_to = "Prevalence_Difference") %>% mutate(Site = recode(Site, "KS_Prev_Diff" = "Kansas", "SD_Prev_Diff" = "South Dakota"))

p_prev_diff_its <- ggplot(prev_diff_plot_its, aes(x = Label, y = Prevalence_Difference, color = Site, shape = Site)) + annotate("rect", xmin = -Inf, xmax = Inf, ymin = 0, ymax = Inf, fill = "#E8F5E9", alpha = 0.35) + annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = 0, fill = "#FDECEC", alpha = 0.35) + geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.8, color = "grey30") + geom_point(size = 5, stroke = 1.2, position = position_dodge(width = 0.45)) + geom_text(aes(label = percent(Prevalence_Difference, accuracy = 1)), position = position_dodge(width = 0.45), vjust = -1.1, size = 4, fontface = "bold", show.legend = FALSE) + scale_color_manual(values = c("Kansas" = "#E07A5F", "South Dakota" = "#3D5A80")) + scale_shape_manual(values = c("Kansas" = 16, "South Dakota" = 17)) + scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(-1, 1), breaks = seq(-1, 1, 0.25), expand = expansion(mult = c(0.02, 0.08))) + labs(x = NULL, y = "Change in prevalence (Early - Regular)", color = "Site", shape = "Site") + theme_classic(base_size = 15) + theme(axis.text.x = element_text(angle = 20, hjust = 1, size = 12, face = "bold"), axis.text.y = element_text(size = 11), axis.title.y = element_text(size = 13, face = "bold"), legend.position = "top", legend.title = element_text(face = "bold"), panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8), plot.margin = margin(10, 20, 10, 10))

p_prev_diff_its

ggsave("ITS_prev_diff.pdf", p_prev_diff_its, width = 8, height = 6, dpi = 1000)

#######################################################################################################################
# STEP 4 — EXPORT ITS RESULTS
#######################################################################################################################

write.csv(cross_site, "ITS_cross_site_DESeq2_results.csv", row.names = FALSE)
write.csv(early_same, "ITS_Early_in_both_sites.csv", row.names = FALSE)
write.csv(early_strict, "ITS_Conserved_Early_significant.csv", row.names = FALSE)
write.csv(regular_same, "ITS_Regular_in_both_sites.csv", row.names = FALSE)
write.csv(regular_strict, "ITS_Conserved_Regular_significant.csv", row.names = FALSE)
write.csv(prev_table, "ITS_ASV_prevalence_table.csv", row.names = FALSE)

