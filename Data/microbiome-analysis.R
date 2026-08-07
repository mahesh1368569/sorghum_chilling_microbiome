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

save.image("ITS-Drought-Cowpea.RData")

load("ITS-Drought-Cowpea.RData")

its_colors_set <- c(
  "#f15a60","#7ac36a","#5a9bda","#faa75b","#ffaaaa",
  "#9e76ab","#c37508","#d77fba",
  "#4bc0c8","#f7c873","#74c2e1","#e07b39","#6cc27c",
  "#c66a9e","#a6b75d","#5b6dd1","#ff8f9f","#b87ec9",
  "#4fa6b7","#d4944a","#76a8dc","#8ccf9f","#e68ac3",
  "#e6a157","#64b5af","#9d84d7","#ffb55e","#6fb3d2",
  "#b55d72","#7ccfcb","#d37f6f","#89c56f","#ce7dcf",
  "#4fa6b7","#d4944a","#76a8dc","#8ccf9f","#e68ac3",
  "#e6a157","#64b5af","#9d84d7","#ffb55e","#6fb3d2",
  "#b55d72","#7ccfcb","#d37f6f","#89c56f","#ce7dcf"
)


##### Importing files ########

biom = import_biom("asv_table_16S_filt_taxonomy.biom")

meta <- read.csv("metadata.csv",header = TRUE,stringsAsFactors = FALSE,check.names = FALSE)

metadata = import_qiime_sample_data("metadata.txt")

sample_names(biom)
sample_names(metadata)

sorghum_biom = merge_phyloseq(biom, metadata)

colnames(tax_table(sorghum_biom)) <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")

sorghum <- phyloseq2meco(sorghum_biom)

sorghum

# ============================================================
# 1. PREPARE MICROECO DATASET
# ============================================================

sorghum <- phyloseq2meco(sorghum_biom)

# Make sure all tables contain matching samples/ASVs
sorghum$tidy_dataset()

# Set factor order
sorghum$sample_table$Site <- factor(
  sorghum$sample_table$Site,
  levels = c("Kansas", "South Dakota")
)

sorghum$sample_table$Planting <- factor(
  sorghum$sample_table$Planting,
  levels = c("Early", "Regular")
)

# Combined Site × Planting variable
sorghum$sample_table$SiteStage <- interaction(
  sorghum$sample_table$Site,
  sorghum$sample_table$Planting,
  sep = " | ",
  drop = TRUE
)

# Check metadata
table(
  sorghum$sample_table$Site,
  sorghum$sample_table$Planting
)

head(sorghum$sample_table)

# ============================================================
# 2. CALCULATE BASIC MICROBIOME METRICS
# ============================================================

# Taxonomic relative abundance
sorghum$cal_abund()

# Alpha diversity
sorghum$cal_alphadiv()

# Beta diversity
# Calculates Bray-Curtis and Jaccard
sorghum$cal_betadiv()

# ============================================================
# 3. ALPHA DIVERSITY
# Stage comparison within each site
# ============================================================

alpha <- trans_alpha$new(
  dataset = sorghum,
  group = "Planting",
  by_group = "Site"
)

# Early vs Regular within each site
alpha$cal_diff(
  method = "wilcox"
)

# See statistics
alpha$res_diff


# -------------------------
# Shannon diversity
# -------------------------

p_alpha_shannon <- alpha$plot_alpha(
  measure = "Shannon",
  add = "jitter",
  jitter_shape = 21
) +
  theme_classic(base_size = 14) +
  labs(
    x = "Planting stage",
    y = "Shannon diversity"
  )

p_alpha_shannon


# -------------------------
# Observed richness
# -------------------------

p_alpha_observed <- alpha$plot_alpha(
  measure = "Observed",
  add = "jitter",
  jitter_shape = 21
) +
  theme_classic(base_size = 14) +
  labs(
    x = "Planting stage",
    y = "Observed ASVs"
  )

p_alpha_observed


# -------------------------
# Chao1
# -------------------------

p_alpha_chao <- alpha$plot_alpha(
  measure = "Chao1",
  add = "jitter",
  jitter_shape = 21
) +
  theme_classic(base_size = 14) +
  labs(
    x = "Planting stage",
    y = "Chao1 richness"
  )

p_alpha_chao


# Combine alpha plots
p_alpha_shannon + p_alpha_observed + p_alpha_chao


# ============================================================
# 4. BETA DIVERSITY — BRAY-CURTIS PCoA
# Color = Site
# Shape = Planting stage
# ============================================================

beta <- trans_beta$new(dataset = sorghum, group = "SiteStage", measure = "bray"
)

beta$cal_ordination(method = "PCoA")

p_beta <- beta$plot_ordination(plot_color = "Site", plot_shape = "Planting", plot_type = c("point", "ellipse"),
  point_size = 3,
  point_alpha = 0.8) +
  scale_color_manual(values = c("Kansas" = "#7195E8","South Dakota" = "#A86480")) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 14), # Increase x-axis text size
        axis.text.y = element_text(size = 14), # Increase y-axis text size
        axis.title.x = element_text(size = 14), # Increase x-axis label size
        axis.title.y = element_text(size = 14), # Increase y-axis label size
        strip.text = element_text(size = 14),
        legend.position = "right",
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 14),# Increase facet label size
        panel.border = element_rect(colour = "black", fill = NA, size = 1))

p_beta

# ============================================================
# 5. PERMANOVA
# Test Site + Planting + interaction
# ============================================================

beta$cal_manova(manova_set = "Site + Planting + Site:Planting")

beta$res_manova
# ============================================================
# 6. CHECK MULTIVARIATE DISPERSION
# Important companion test for PERMANOVA
# ============================================================

beta$cal_betadisper()

beta$res_betadisper
# ============================================================
# 7. TAXONOMIC ABUNDANCE — PHYLUM
# Mean abundance for each Site × Planting combination
# ============================================================

table(sorghum$sample_table$Site)
table(sorghum$sample_table$Planting)
table(sorghum$sample_table$Site, sorghum$sample_table$Planting)

sorghum$sample_table$SitePlanting <- paste(
  sorghum$sample_table$Site,
  sorghum$sample_table$Planting,
  sep = " "
)

table(sorghum$sample_table$SitePlanting, useNA = "ifany")

sorghum$sample_table$SitePlanting <- paste(
  sorghum$sample_table$Site,
  sorghum$sample_table$Planting,
  sep = " "
)

table(sorghum$sample_table$SitePlanting, useNA = "ifany")

abund_phylum <- trans_abund$new(dataset = sorghum, taxrank = "Phylum", ntaxa = 10, groupmean = "SitePlanting")

p_phylum <- abund_phylum$plot_bar(others_color = "grey70", legend_text_italic = FALSE) +
  theme_classic(base_size = 14) + labs(
    x = NULL,
    y = "Relative abundance (%)",
    fill = "Phylum") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45,hjust = 1,size = 14),
        axis.text.y = element_text(size = 14),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        strip.text = element_text(size = 14),
        legend.position = "right",
        legend.title = element_text(size = 13),
        legend.text = element_text(size = 12),
        panel.border = element_rect(colour = "black",fill = NA, size = 1))


p_phylum

abund_class <- trans_abund$new(dataset = sorghum, taxrank = "Class", ntaxa = 20, groupmean = "SitePlanting")

p_class <- abund_class$plot_bar(others_color = "grey70", legend_text_italic = FALSE) +
  theme_classic(base_size = 14) + labs(
    x = NULL,
    y = "Relative abundance (%)",
    fill = "Phylum") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45,hjust = 1,size = 14),
        axis.text.y = element_text(size = 14),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        strip.text = element_text(size = 14),
        legend.position = "right",
        legend.title = element_text(size = 13),
        legend.text = element_text(size = 12),
        panel.border = element_rect(colour = "black",fill = NA, size = 1))


p_class

# ============================================================
# 10. NETWORK ANALYSIS
# Sorghum 16S
# Four networks:
# Kansas Early
# Kansas Regular
# South Dakota Early
# South Dakota Regular
# ============================================================

#----------------------------------------------------------------------
# 1. Create network list
#----------------------------------------------------------------------
sorghum_network = list()
#----------------------------------------------------------------------
# 2. Clone datasets
#----------------------------------------------------------------------
kansas_early = clone(sorghum)
kansas_regular = clone(sorghum)

sd_early = clone(sorghum)
sd_regular = clone(sorghum)
#----------------------------------------------------------------------
# 3. Subset by Site and Planting
#----------------------------------------------------------------------
kansas_early$sample_table %<>% subset(Site == "Kansas" & Planting == "Early")

kansas_regular$sample_table %<>% subset(Site == "Kansas" & Planting == "Regular")

sd_early$sample_table %<>% subset(Site == "South Dakota" & Planting == "Early")

sd_regular$sample_table %<>% subset(Site == "South Dakota" & Planting == "Regular")
#----------------------------------------------------------------------
# 4. Check sample numbers
#----------------------------------------------------------------------
nrow(kansas_early$sample_table)
nrow(kansas_regular$sample_table)

nrow(sd_early$sample_table)
nrow(sd_regular$sample_table)
#----------------------------------------------------------------------
# 5. Trim dataset after subsetting
#----------------------------------------------------------------------
kansas_early$tidy_dataset()
kansas_regular$tidy_dataset()

sd_early$tidy_dataset()
sd_regular$tidy_dataset()
#----------------------------------------------------------------------
# 6. Check taxa and samples
#----------------------------------------------------------------------
kansas_early
kansas_regular

sd_early
sd_regular
#----------------------------------------------------------------------
# 7. Construct networks
#
# filter_thres = 0.0005
# removes ASVs with very low relative abundance
#
# Spearman correlation
#----------------------------------------------------------------------
kansas_early_network <- trans_network$new( dataset = kansas_early, cor_method = "spearman", filter_thres = 0.0005)

kansas_regular_network <- trans_network$new(dataset = kansas_regular, cor_method = "spearman", filter_thres = 0.0005)

sd_early_network <- trans_network$new(dataset = sd_early, cor_method = "spearman", filter_thres = 0.0005)

sd_regular_network <- trans_network$new(dataset = sd_regular,cor_method = "spearman",filter_thres = 0.0005)
#----------------------------------------------------------------------
# 8. Calculate network
#
# COR_p_thres = significance threshold
# COR_cut = correlation strength threshold
#----------------------------------------------------------------------
kansas_early_network$cal_network(COR_p_thres = 0.05, COR_cut = 0.6)

kansas_regular_network$cal_network(COR_p_thres = 0.05, COR_cut = 0.6)

sd_early_network$cal_network(COR_p_thres = 0.1, COR_cut = 0.6)

sd_regular_network$cal_network(COR_p_thres = 0.05, COR_cut = 0.6)
#----------------------------------------------------------------------
# 9. Put networks into one list
#----------------------------------------------------------------------

sorghum_network$Kansas_Early = kansas_early_network
sorghum_network$Kansas_Regular = kansas_regular_network

sorghum_network$SouthDakota_Early = sd_early_network
sorghum_network$SouthDakota_Regular = sd_regular_network
# Check
sorghum_network
#======================================================================
# 10. NETWORK MODULARITY
#======================================================================
sorghum_network %<>%cal_module(undirected_method = "cluster_fast_greedy")
#======================================================================
# 11. NODE AND EDGE TABLES
#======================================================================
sorghum_network %<>%get_node_table(node_roles = TRUE) %>%get_edge_table
#======================================================================
# 12. NETWORK TOPOLOGICAL ATTRIBUTES
#======================================================================
network_atr = cal_network_attr(sorghum_network)
network_atr

write.csv(network_atr, "Sorghum_network_attributes.csv")
#======================================================================
# 13. NODE DEGREE
#======================================================================

degree_sorghum = node_comp( sorghum_network, property = "degree")

degree_df <- as.data.frame(degree_sorghum$otu_table)

degree_df$Node <- rownames(degree_df)

degree_long <- pivot_longer(degree_df,cols = -Node,names_to = "Network", values_to = "Degree")

head(degree_long)
#----------------------------------------------------------------------
# Degree boxplot
#----------------------------------------------------------------------
degree_plot = ggplot(degree_long, aes(x = Network, y = Degree,fill = Network)) +
  geom_boxplot() +
  labs(title = "Node Degree Across Sorghum Networks",
    x = "Network",
    y = "Node Degree") +
  theme_minimal() +
  scale_fill_brewer(
  palette = "Set1") +
  theme(anel.border = element_rect(color = "black", fill = NA, size = 1.5),
    axis.text.x = element_text(angle = 30,hjust = 1,size = 12))

degree_plot

ggsave("Sorghum_network_degree.pdf",degree_plot, width = 7, height = 5)
#----------------------------------------------------------------------
# Test differences in degree
#----------------------------------------------------------------------
kruskal.test(Degree ~ Network,data = degree_long)

# Pairwise comparison if overall test significant
pairwise.wilcox.test(degree_long$Degree, degree_long$Network, p.adjust.method = "BH")

#======================================================================
# 14. BETWEENNESS CENTRALITY
#======================================================================
btw_table <- node_comp(sorghum_network, property = "betweenness_centrality")

#======================================================================
# 15. CLOSENESS CENTRALITY
#======================================================================
close_table <- node_comp(sorghum_network,property = "closeness_centrality")
#======================================================================
# 16. EIGENVECTOR CENTRALITY
#======================================================================
eigen_table <- node_comp(sorghum_network,property = "eigenvector_centrality")
#======================================================================
# 17. Convert centrality results to data frames
#======================================================================

btw_df <- btw_table$otu_table %>% as.data.frame() %>% tibble::rownames_to_column("Node")

close_df <- close_table$otu_table %>% as.data.frame() %>% tibble::rownames_to_column("Node")

eigen_df <- eigen_table$otu_table %>% as.data.frame() %>% tibble::rownames_to_column("Node")

#======================================================================
# 18. Convert to long format
#======================================================================

btw_long <- btw_df %>% pivot_longer( -Node, names_to = "Network", values_to = "Value") %>%
  mutate(Measure = "Betweenness")

close_long <- close_df %>% pivot_longer( -Node, names_to = "Network", values_to = "Value") %>%
  mutate(Measure = "Closeness")

eigen_long <- eigen_df %>% pivot_longer(-Node, names_to = "Network", values_to = "Value") %>%
  mutate(Measure = "Eigenvector")
#----------------------------------------------------------------------
# Combine all centrality measurements
#----------------------------------------------------------------------

centrality_long <- bind_rows(btw_long, close_long, eigen_long)
#======================================================================
# 19. CENTRALITY PLOT
#======================================================================

centrality_plot = ggplot( centrality_long,
  aes(x = Network, y = Value, fill = Network)) +
  geom_boxplot() +
  facet_wrap( ~ Measure, scales = "free_y") +
  labs(title = "Centrality Measures Across Sorghum Networks",
    y = "Centrality Value",
    x = "Network") +
  theme_minimal() +
  scale_fill_brewer(palette = "Set1") +
  theme(panel.border = element_rect( color = "black", fill = NA, size = 1.5),
    axis.text.x = element_text(angle = 30, hjust = 1, size = 11)) +
  stat_compare_means(method = "kruskal.test", label = "p.format")

centrality_plot


ggsave("Sorghum_network_centrality.pdf", centrality_plot, width = 8, height = 6)
#======================================================================
# 20. PLOT INDIVIDUAL NETWORKS
# Nodes colored by Phylum
#======================================================================

ks_early_plot = sorghum_network[["Kansas_Early"]]$plot_network(method = "ggraph", node_color = "Phylum")

ks_regular_plot =sorghum_network[["Kansas_Regular"]]$plot_network(method = "ggraph",node_color = "Phylum")

sd_early_plot =sorghum_network[["SouthDakota_Early"]]$plot_network(method = "ggraph",node_color = "Phylum")

sd_regular_plot =sorghum_network[["SouthDakota_Regular"]]$plot_network(method = "ggraph",node_color = "Phylum")

ks_early_plot
ks_regular_plot
sd_early_plot
sd_regular_plot


#======================================================================
# 21. SAVE NETWORK FIGURES
#======================================================================

ggsave("Kansas_Early_network.pdf",ks_early_plot,height = 6,width = 6,dpi = 1000)
ggsave("Kansas_Regular_network.pdf",ks_regular_plot,height = 6,width = 6,dpi = 1000)
ggsave("SouthDakota_Early_network.pdf",sd_early_plot,height = 6,width = 6,dpi = 1000)
ggsave("SouthDakota_Regular_network.pdf",sd_regular_plot,height = 6,width = 6,dpi = 1000)


#======================================================================
# 22. COMPARE NODES ACROSS NETWORKS
#======================================================================
node_dist <- node_comp(sorghum_network,property = "name")
node_dist
#----------------------------------------------------------------------
# Venn / node intersection
#----------------------------------------------------------------------

node_intersection <- trans_venn$new(node_dist,ratio = "numratio")

node_intersection

node_intersection_plot <-node_intersection$plot_venn

node_intersection_plot

#======================================================================
# 23. EDGE TAXONOMIC COMPOSITION
#
# Which phyla are connected to one another?
# "+" means positive associations
#======================================================================

sorghum_edgetax = edge_tax_comp(sorghum_network,taxrank = "Phylum", label = "+", rel = TRUE)

#----------------------------------------------------------------------
# Remove very rare edge categories
#----------------------------------------------------------------------

sorghum_edgetax =sorghum_edgetax[apply(sorghum_edgetax,1,mean) > 0.01,]
#----------------------------------------------------------------------
# Heatmap
#----------------------------------------------------------------------

edge_heatmap = pheatmap::pheatmap(sorghum_edgetax,display_numbers = TRUE)
#======================================================================
# 24. NEGATIVE EDGE COMPOSITION
#======================================================================

sorghum_edgetax_negative = edge_tax_comp(sorghum_network,taxrank = "Phylum",label = "-",rel = TRUE)

sorghum_edgetax_negative =sorghum_edgetax_negative[apply(sorghum_edgetax_negative,1,mean) > 0.01,]

negative_edge_heatmap = pheatmap::pheatmap(sorghum_edgetax_negative, display_numbers = TRUE)

#======================================================================
# 25. NETWORK ROBUSTNESS
#======================================================================

robustness_sorghum <- robustness$new(sorghum_network,
  remove_strategy = c("edge_rand","edge_strong","node_rand","node_degree_high"),
  remove_ratio = seq(0,0.99, 0.1),
  measure = c("Eff","Eigen","Pcr"),
  run = 10)
#----------------------------------------------------------------------
# Plot robustness
#----------------------------------------------------------------------

plotrob_sorghum =robustness_sorghum$plot(linewidth = 1) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45,hjust = 1,size = 14),
    axis.text.y = element_text(size = 14),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    strip.text = element_text(size = 14),
    legend.position = "right",
    legend.title = element_text(size = 13),
    legend.text = element_text(size = 12),
    panel.border = element_rect(colour = "black",fill = NA, size = 1))

plotrob_sorghum

ggsave("Sorghum_network_robustness.pdf",plotrob_sorghum,width = 14,height = 10)

#======================================================================
# 26. NETWORK VULNERABILITY
#======================================================================

vul_table <- vulnerability(sorghum_network)

head(vul_table)
#----------------------------------------------------------------------
# Vulnerability plot
#----------------------------------------------------------------------

vulnerability_plot =ggplot(data = vul_table,aes(x = Network,y = vulnerability,fill = Network)) +
  geom_boxplot() +
  labs(title = "Network Vulnerability", x = "Network", y = "Vulnerability") +
  theme_minimal() + scale_fill_brewer(palette = "Set1") +
  theme(panel.border = element_rect(color = "black",fill = NA,size = 1.5),
    axis.text.x = element_text(angle = 30, hjust = 1))

vulnerability_plot

ggsave("Sorghum_network_vulnerability.pdf",vulnerability_plot,width = 7,height = 6)
