library(phyloseq)
library(biomformat)
library(file2meco)
library(microeco)
library(ggplot2)
library(colorspace)

##### IMPORT ITS BIOM + METADATA #####

its_biom <- import_biom("ASV_ITS_filt_taxonomy.biom")
its_metadata <- import_qiime_sample_data("metadata.txt")
sample_names(its_biom)
sample_names(its_metadata)
sorghum_its_biom <- merge_phyloseq(its_biom, its_metadata)
colnames(tax_table(sorghum_its_biom)) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
sorghum_its <- phyloseq2meco(sorghum_its_biom)
sorghum_its

##### CHECK METADATA #####

table(sorghum_its$sample_table$Site)
table(sorghum_its$sample_table$Planting)
table(sorghum_its$sample_table$Site, sorghum_its$sample_table$Planting)
table(sorghum_its$sample_table$Line)

##### CREATE SITE × PLANTING #####

sorghum_its$sample_table$SitePlanting <- paste(sorghum_its$sample_table$Site, sorghum_its$sample_table$Planting, sep = " ")
sorghum_its$sample_table$SitePlanting <- factor(sorghum_its$sample_table$SitePlanting, levels = c("Kansas Early", "Kansas Regular", "South Dakota Early", "South Dakota Regular"))
sorghum_its$sample_table$Line <- as.character(sorghum_its$sample_table$Line)
table(sorghum_its$sample_table$SitePlanting, useNA = "ifany")
table(sorghum_its$sample_table$Line, sorghum_its$sample_table$SitePlanting)

##### PHYLUM ABUNDANCE BY SITE × PLANTING #####

abund_its_phylum <- trans_abund$new(dataset = sorghum_its, taxrank = "Phylum", ntaxa = 10, groupmean = "SitePlanting")
unique(abund_its_phylum$data_abund$Taxonomy)
its_taxa <- unique(as.character(abund_its_phylum$data_abund$Taxonomy))
its_colorset_fungi <- setNames(colorspace::qualitative_hcl(length(its_taxa), palette = "Dark 3"), its_taxa)
its_colorset_fungi["Others"] <- "#BDBDBD"

p_its_phylum <- abund_its_phylum$plot_bar(others_color = "#BDBDBD", legend_text_italic = FALSE) + scale_fill_manual(values = its_colorset_fungi) + labs(x = NULL, y = "Relative abundance (%)", fill = "Phylum") + theme_minimal() + theme(axis.text.x = element_text(angle = 25, hjust = 1, size = 14, colour = "black"), axis.text.y = element_text(size = 14, colour = "black"), axis.title.y = element_text(size = 14, colour = "black"), strip.text = element_text(size = 14), legend.position = "right", legend.title = element_text(size = 13), legend.text = element_text(size = 12), panel.border = element_rect(colour = "black", fill = NA, linewidth = 1))
p_its_phylum

##### PHYLUM ABUNDANCE BY LINE #####

abund_its_phylum_L <- trans_abund$new(dataset = sorghum_its, taxrank = "Phylum", ntaxa = 10, groupmean = "Line")
p_its_phylum_l <- abund_its_phylum_L$plot_bar(others_color = "#BDBDBD", legend_text_italic = FALSE) + scale_fill_manual(values = its_colorset_fungi) + labs(x = NULL, y = "Relative abundance (%)", fill = "Phylum") + theme_minimal() + theme(axis.text.x = element_text(angle = 25, hjust = 1, size = 14, colour = "black"), axis.text.y = element_text(size = 14, colour = "black"), axis.title.y = element_text(size = 14, colour = "black"), legend.position = "right", legend.title = element_text(size = 13), legend.text = element_text(size = 12), panel.border = element_rect(colour = "black", fill = NA, linewidth = 1))
p_its_phylum_l

##### CLASS ABUNDANCE #####

abund_its_class <- trans_abund$new(dataset = sorghum_its, taxrank = "Class", ntaxa = 20, groupmean = "SitePlanting")
p_its_class <- abund_its_class$plot_bar(others_color = "grey70", legend_text_italic = FALSE) + labs(x = NULL, y = "Relative abundance (%)", fill = "Class") + theme_minimal() + theme(axis.text.x = element_text(angle = 25, hjust = 1, size = 14, colour = "black"), axis.text.y = element_text(size = 14, colour = "black"), axis.title.y = element_text(size = 14, colour = "black"), legend.position = "right", legend.title = element_text(size = 13), legend.text = element_text(size = 12), panel.border = element_rect(colour = "black", fill = NA, linewidth = 1))
p_its_class

##### CREATE ABUNDANCE FOLDERS #####

dir.create("abundance", showWarnings = FALSE)
dir.create(file.path("abundance", "ITS"), recursive = TRUE, showWarnings = FALSE)

##### SAVE MAIN ABUNDANCE PLOTS #####

ggsave(file.path("abundance", "ITS", "ITS_Phylum_SitePlanting.pdf"), p_its_phylum, width = 7, height = 6)
ggsave(file.path("abundance", "ITS", "ITS_Phylum_Line.pdf"), p_its_phylum_l, width = 9, height = 6)
ggsave(file.path("abundance", "ITS", "ITS_Class_SitePlanting.pdf"), p_its_class, width = 7, height = 6)

##### SEPARATE PHYLUM PLOT FOR EACH SORGHUM LINE #####

abund_its_line <- trans_abund$new(dataset = sorghum_its, taxrank = "Phylum", ntaxa = 10)
its_phylum_dat <- abund_its_line$data_abund
its_phylum_dat$Line <- as.character(its_phylum_dat$Line)
its_phylum_dat$SitePlanting <- paste(its_phylum_dat$Site, its_phylum_dat$Planting)
its_phylum_dat <- its_phylum_dat[, c("Sample", "Line", "SitePlanting", "Taxonomy", "Abundance")]
unique(its_phylum_dat$SitePlanting)
unique(its_phylum_dat$Taxonomy)

##### CALCULATE OTHERS #####

its_sample_sum <- aggregate(Abundance ~ Sample + Line + SitePlanting, data = its_phylum_dat, FUN = sum)
its_sample_sum$Taxonomy <- "Others"
its_sample_sum$Abundance <- 100 - its_sample_sum$Abundance
its_sample_sum$Abundance[its_sample_sum$Abundance < 0] <- 0
its_others_dat <- its_sample_sum[, c("Sample", "Line", "SitePlanting", "Taxonomy", "Abundance")]
its_phylum_all <- rbind(its_phylum_dat, its_others_dat)

##### MEAN ABUNDANCE FOR EACH LINE × SITE × PLANTING #####

its_phylum_mean <- aggregate(Abundance ~ Line + SitePlanting + Taxonomy, data = its_phylum_all, FUN = mean)
its_phylum_mean$SitePlanting <- factor(its_phylum_mean$SitePlanting, levels = c("Kansas Early", "Kansas Regular", "South Dakota Early", "South Dakota Regular"))
table(its_phylum_mean$Line, its_phylum_mean$SitePlanting)

##### CONSISTENT COLORS FOR LINE PLOTS #####

its_line_taxa <- unique(as.character(its_phylum_mean$Taxonomy))
its_plot_colors <- setNames(colorspace::qualitative_hcl(length(its_line_taxa), palette = "Dark 3"), its_line_taxa)
its_plot_colors["Others"] <- "#BDBDBD"

##### GENERATE ONE PLOT FOR EVERY LINE #####

its_line_names <- sort(unique(its_phylum_mean$Line))
its_phylum_line_plots <- list()

for (i in its_line_names) { df_i <- its_phylum_mean[its_phylum_mean$Line == i, ]; df_i$SitePlanting <- factor(df_i$SitePlanting, levels = c("Kansas Early", "Kansas Regular", "South Dakota Early", "South Dakota Regular")); p <- ggplot(df_i, aes(x = SitePlanting, y = Abundance, fill = Taxonomy)) + geom_col(width = 0.85) + scale_fill_manual(values = its_plot_colors) + scale_y_continuous(limits = c(0, 100), expand = c(0, 0)) + labs(title = i, x = NULL, y = "Relative abundance (%)", fill = "Phylum") + theme_minimal() + theme(axis.text.x = element_text(angle = 25, hjust = 1, size = 14, colour = "black"), axis.text.y = element_text(size = 14, colour = "black"), axis.title.y = element_text(size = 14, colour = "black"), plot.title = element_text(size = 15, face = "bold", hjust = 0.5), legend.position = "right", legend.title = element_text(size = 13), legend.text = element_text(size = 12), panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(), panel.border = element_rect(colour = "black", fill = NA, linewidth = 1)); its_phylum_line_plots[[i]] <- p }

##### CHECK GENERATED PLOTS #####

names(its_phylum_line_plots)
its_phylum_line_plots[["PI576376"]]

##### SAVE EVERY LINE PLOT #####

for (i in names(its_phylum_line_plots)) { ggsave(filename = file.path("abundance", "ITS", paste0(i, "_ITS_Phylum_Abundance.pdf")), plot = its_phylum_line_plots[[i]], width = 8, height = 6, units = "in") }

##### CHECK SAVED FILES #####

list.files(file.path("abundance", "ITS"))


# ============================================================
# RDA or dbRDA
# ============================================================

library(readxl)

plant <- read_excel ("plant-data.xlsx") %>% as.data.frame()

rownames(plant) <- plant[, 1]

plant = plant[ ,-1]

envits <- trans_env$new(dataset = sorghum_its, add_data = plant)

envits$cal_ordination(method = "RDA", taxa_level = "Genus")

envits$cal_ordination_anova()

envits$cal_ordination_envfit()

envits$trans_ordination(show_taxa = 10, adjust_arrow_length = TRUE, max_perc_env = 0.5, max_perc_tax = 0.5)

p_rdaits <- envits$plot_ordination(plot_color = "Planting", plot_shape = "Site")

p_rdaits

ggsave("rda_ITS.pdf", p_rdaits, w=7, h = 7, dpi = 1000)

envits$res_ordination_R2

envits$res_ordination_terms

envits$res_ordination_envfit

envits$cal_mantel(use_measure = "bray")

envits$res_mantel

env16$cal_mantel(by_group = "Site", use_measure = "bray")

env16$res_mantel
