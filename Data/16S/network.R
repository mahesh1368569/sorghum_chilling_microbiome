library(tidyverse)
library(Hmisc)
library(igraph)
library(ggraph)
library(patchwork)
library(scales)
library(ggrepel)

output_dir <- file.path(getwd(), "Network_Analysis")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

asv <- read.csv("asv_table_16S_filt.csv", check.names = FALSE)
tax <- read.csv("asv_taxa_16S_filt.csv", check.names = FALSE)

colnames(asv)[1] <- "Sample_ID"
colnames(tax)[1] <- "ASV"
tax <- tax %>% select(-X)

dim(asv)
dim(tax)
head(asv[,1:6])
head(tax)

rownames(asv) <- asv$Sample_ID
asv_counts <- as.matrix(asv[,-1])
storage.mode(asv_counts) <- "numeric"

all(colnames(asv_counts) %in% tax$ASV)
all(tax$ASV %in% colnames(asv_counts))

metadata <- tibble(Sample_ID = rownames(asv_counts)) %>% mutate(Site = case_when(str_detect(Sample_ID, "^KS_") ~ "Kansas", str_detect(Sample_ID, "^SD_") ~ "South Dakota", TRUE ~ NA_character_), Planting = case_when(str_detect(Sample_ID, "_E_") ~ "Early", str_detect(Sample_ID, "_R_") ~ "Regular", TRUE ~ NA_character_), SitePlanting = paste(Site, Planting))

metadata$SitePlanting <- factor(metadata$SitePlanting, levels = c("Kansas Early", "Kansas Regular", "South Dakota Early", "South Dakota Regular"))

table(metadata$SitePlanting)

rel_abund <- sweep(asv_counts, 1, rowSums(asv_counts), FUN = "/")
rel_abund[is.na(rel_abund)] <- 0

rho_cutoff <- 0.60
q_cutoff <- 0.05
set.seed(123)

phylum_levels <- sort(unique(na.omit(tax$Phylum)))
phylum_colors <- scales::hue_pal()(length(phylum_levels))
names(phylum_colors) <- phylum_levels

##### Kansas Early ########

ks_e_samples <- metadata %>% filter(SitePlanting == "Kansas Early") %>% pull(Sample_ID)
ks_e_mat <- rel_abund[ks_e_samples,,drop = FALSE]
ks_e_mat <- ks_e_mat[,colSums(ks_e_mat) > 0 & apply(ks_e_mat, 2, sd, na.rm = TRUE) > 0,drop = FALSE]

dim(ks_e_mat)

ks_e_cor <- Hmisc::rcorr(as.matrix(ks_e_mat), type = "spearman")
ks_e_idx <- which(upper.tri(ks_e_cor$r), arr.ind = TRUE)
ks_e_edges <- tibble(from = colnames(ks_e_cor$r)[ks_e_idx[,1]], to = colnames(ks_e_cor$r)[ks_e_idx[,2]], rho = ks_e_cor$r[ks_e_idx], p = ks_e_cor$P[ks_e_idx]) %>% filter(!is.na(rho), !is.na(p)) %>% mutate(q = p.adjust(p, method = "BH")) %>% filter(abs(rho) >= rho_cutoff, q <= q_cutoff) %>% mutate(edge_sign = if_else(rho > 0, "Positive", "Negative"))

nrow(ks_e_edges)
table(ks_e_edges$edge_sign)

ks_e_nodes <- tibble(name = sort(unique(c(ks_e_edges$from, ks_e_edges$to)))) %>% left_join(tax %>% select(ASV, Kingdom, Phylum, Class, Order, Family, Genus), by = c("name" = "ASV"))
ks_e_nodes$mean_relative_abundance <- colMeans(ks_e_mat[,ks_e_nodes$name,drop = FALSE], na.rm = TRUE)

g_ks_e <- graph_from_data_frame(ks_e_edges, directed = FALSE, vertices = ks_e_nodes)

V(g_ks_e)$degree <- degree(g_ks_e)
V(g_ks_e)$betweenness <- betweenness(g_ks_e, directed = FALSE, normalized = FALSE)
V(g_ks_e)$normalized_degree <- degree(g_ks_e) / max(1, vcount(g_ks_e) - 1)
V(g_ks_e)$module <- membership(cluster_louvain(g_ks_e))

set.seed(123)
p_ks_e <- ggraph(g_ks_e, layout = "fr") + geom_edge_link(aes(color = edge_sign), alpha = 0.30, linewidth = 0.35) + geom_node_point(aes(size = degree, color = Phylum), alpha = 0.95) + scale_edge_color_manual(values = c("Positive" = "grey40", "Negative" = "firebrick")) + scale_color_manual(values = phylum_colors, drop = FALSE) + scale_size_continuous(range = c(2,9)) + labs(title = "Kansas Early", edge_color = "Association", color = "Phylum", size = "Degree") + theme_void() + theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 16), legend.position = "right")

p_ks_e

ggsave(file.path(output_dir, "Kansas_Early_network.pdf"), p_ks_e, width = 8, height = 7)
write.csv(ks_e_edges, file.path(output_dir, "Kansas_Early_edges.csv"), row.names = FALSE)
write.csv(as.data.frame(vertex_attr(g_ks_e)), file.path(output_dir, "Kansas_Early_nodes.csv"), row.names = FALSE)

### Kansas Regular #####

ks_r_samples <- metadata %>% filter(SitePlanting == "Kansas Regular") %>% pull(Sample_ID)
ks_r_mat <- rel_abund[ks_r_samples,,drop = FALSE]
ks_r_mat <- ks_r_mat[,colSums(ks_r_mat) > 0 & apply(ks_r_mat, 2, sd, na.rm = TRUE) > 0,drop = FALSE]

ks_r_cor <- Hmisc::rcorr(as.matrix(ks_r_mat), type = "spearman")
ks_r_idx <- which(upper.tri(ks_r_cor$r), arr.ind = TRUE)
ks_r_edges <- tibble(from = colnames(ks_r_cor$r)[ks_r_idx[,1]], to = colnames(ks_r_cor$r)[ks_r_idx[,2]], rho = ks_r_cor$r[ks_r_idx], p = ks_r_cor$P[ks_r_idx]) %>% filter(!is.na(rho), !is.na(p)) %>% mutate(q = p.adjust(p, method = "BH")) %>% filter(abs(rho) >= rho_cutoff, q <= q_cutoff) %>% mutate(edge_sign = if_else(rho > 0, "Positive", "Negative"))

table(ks_r_edges$edge_sign)

ks_r_nodes <- tibble(name = sort(unique(c(ks_r_edges$from, ks_r_edges$to)))) %>% left_join(tax %>% select(ASV, Kingdom, Phylum, Class, Order, Family, Genus), by = c("name" = "ASV"))
ks_r_nodes$mean_relative_abundance <- colMeans(ks_r_mat[,ks_r_nodes$name,drop = FALSE], na.rm = TRUE)

g_ks_r <- graph_from_data_frame(ks_r_edges, directed = FALSE, vertices = ks_r_nodes)

V(g_ks_r)$degree <- degree(g_ks_r)
V(g_ks_r)$betweenness <- betweenness(g_ks_r, directed = FALSE, normalized = FALSE)
V(g_ks_r)$normalized_degree <- degree(g_ks_r) / max(1, vcount(g_ks_r) - 1)
V(g_ks_r)$module <- membership(cluster_louvain(g_ks_r))

set.seed(123)
p_ks_r <- ggraph(g_ks_r, layout = "fr") + geom_edge_link(aes(color = edge_sign), alpha = 0.30, linewidth = 0.35) + geom_node_point(aes(size = degree, color = Phylum), alpha = 0.95) + scale_edge_color_manual(values = c("Positive" = "grey40", "Negative" = "firebrick")) + scale_color_manual(values = phylum_colors, drop = FALSE) + scale_size_continuous(range = c(2,9)) + labs(title = "Kansas Regular", edge_color = "Association", color = "Phylum", size = "Degree") + theme_void() + theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 16), legend.position = "right")

p_ks_r

ggsave(file.path(output_dir, "Kansas_Regular_network.pdf"), p_ks_r, width = 8, height = 7)
write.csv(ks_r_edges, file.path(output_dir, "Kansas_Regular_edges.csv"), row.names = FALSE)
write.csv(as.data.frame(vertex_attr(g_ks_r)), file.path(output_dir, "Kansas_Regular_nodes.csv"), row.names = FALSE)

### South Dakota Early #####

sd_e_samples <- metadata %>% filter(SitePlanting == "South Dakota Early") %>% pull(Sample_ID)
sd_e_mat <- rel_abund[sd_e_samples,,drop = FALSE]
sd_e_mat <- sd_e_mat[,colSums(sd_e_mat) > 0 & apply(sd_e_mat, 2, sd, na.rm = TRUE) > 0,drop = FALSE]

sd_e_cor <- Hmisc::rcorr(as.matrix(sd_e_mat), type = "spearman")
sd_e_idx <- which(upper.tri(sd_e_cor$r), arr.ind = TRUE)
sd_e_edges <- tibble(from = colnames(sd_e_cor$r)[sd_e_idx[,1]], to = colnames(sd_e_cor$r)[sd_e_idx[,2]], rho = sd_e_cor$r[sd_e_idx], p = sd_e_cor$P[sd_e_idx]) %>% filter(!is.na(rho), !is.na(p)) %>% mutate(q = p.adjust(p, method = "BH")) %>% filter(abs(rho) >= rho_cutoff, q <= q_cutoff) %>% mutate(edge_sign = if_else(rho > 0, "Positive", "Negative"))

table(sd_e_edges$edge_sign)

sd_e_nodes <- tibble(name = sort(unique(c(sd_e_edges$from, sd_e_edges$to)))) %>% left_join(tax %>% select(ASV, Kingdom, Phylum, Class, Order, Family, Genus), by = c("name" = "ASV"))
sd_e_nodes$mean_relative_abundance <- colMeans(sd_e_mat[,sd_e_nodes$name,drop = FALSE], na.rm = TRUE)

g_sd_e <- graph_from_data_frame(sd_e_edges, directed = FALSE, vertices = sd_e_nodes)

V(g_sd_e)$degree <- degree(g_sd_e)
V(g_sd_e)$betweenness <- betweenness(g_sd_e, directed = FALSE, normalized = FALSE)
V(g_sd_e)$normalized_degree <- degree(g_sd_e) / max(1, vcount(g_sd_e) - 1)
V(g_sd_e)$module <- membership(cluster_louvain(g_sd_e))

set.seed(123)
p_sd_e <- ggraph(g_sd_e, layout = "fr") + geom_edge_link(aes(color = edge_sign), alpha = 0.30, linewidth = 0.35) + geom_node_point(aes(size = degree, color = Phylum), alpha = 0.95) + scale_edge_color_manual(values = c("Positive" = "grey40", "Negative" = "firebrick")) + scale_color_manual(values = phylum_colors, drop = FALSE) + scale_size_continuous(range = c(2,9)) + labs(title = "South Dakota Early", edge_color = "Association", color = "Phylum", size = "Degree") + theme_void() + theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 16), legend.position = "right")

p_sd_e

ggsave(file.path(output_dir, "South_Dakota_Early_network.pdf"), p_sd_e, width = 8, height = 7)
write.csv(sd_e_edges, file.path(output_dir, "South_Dakota_Early_edges.csv"), row.names = FALSE)
write.csv(as.data.frame(vertex_attr(g_sd_e)), file.path(output_dir, "South_Dakota_Early_nodes.csv"), row.names = FALSE)


##### South Dakota Regular #######
sd_r_samples <- metadata %>% filter(SitePlanting == "South Dakota Regular") %>% pull(Sample_ID)
sd_r_mat <- rel_abund[sd_r_samples,,drop = FALSE]
sd_r_mat <- sd_r_mat[,colSums(sd_r_mat) > 0 & apply(sd_r_mat, 2, sd, na.rm = TRUE) > 0,drop = FALSE]

sd_r_cor <- Hmisc::rcorr(as.matrix(sd_r_mat), type = "spearman")
sd_r_idx <- which(upper.tri(sd_r_cor$r), arr.ind = TRUE)
sd_r_edges <- tibble(from = colnames(sd_r_cor$r)[sd_r_idx[,1]], to = colnames(sd_r_cor$r)[sd_r_idx[,2]], rho = sd_r_cor$r[sd_r_idx], p = sd_r_cor$P[sd_r_idx]) %>% filter(!is.na(rho), !is.na(p)) %>% mutate(q = p.adjust(p, method = "BH")) %>% filter(abs(rho) >= rho_cutoff, q <= q_cutoff) %>% mutate(edge_sign = if_else(rho > 0, "Positive", "Negative"))

table(sd_r_edges$edge_sign)

sd_r_nodes <- tibble(name = sort(unique(c(sd_r_edges$from, sd_r_edges$to)))) %>% left_join(tax %>% select(ASV, Kingdom, Phylum, Class, Order, Family, Genus), by = c("name" = "ASV"))
sd_r_nodes$mean_relative_abundance <- colMeans(sd_r_mat[,sd_r_nodes$name,drop = FALSE], na.rm = TRUE)

g_sd_r <- graph_from_data_frame(sd_r_edges, directed = FALSE, vertices = sd_r_nodes)

V(g_sd_r)$degree <- degree(g_sd_r)
V(g_sd_r)$betweenness <- betweenness(g_sd_r, directed = FALSE, normalized = FALSE)
V(g_sd_r)$normalized_degree <- degree(g_sd_r) / max(1, vcount(g_sd_r) - 1)
V(g_sd_r)$module <- membership(cluster_louvain(g_sd_r))

set.seed(123)
p_sd_r <- ggraph(g_sd_r, layout = "fr") + geom_edge_link(aes(color = edge_sign), alpha = 0.30, linewidth = 0.35) + geom_node_point(aes(size = degree, color = Phylum), alpha = 0.95) + scale_edge_color_manual(values = c("Positive" = "grey40", "Negative" = "firebrick")) + scale_color_manual(values = phylum_colors, drop = FALSE) + scale_size_continuous(range = c(2,9)) + labs(title = "South Dakota Regular", edge_color = "Association", color = "Phylum", size = "Degree") + theme_void() + theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 16), legend.position = "right")

p_sd_r

ggsave(file.path(output_dir, "South_Dakota_Regular_network.pdf"), p_sd_r, width = 8, height = 7)
write.csv(sd_r_edges, file.path(output_dir, "South_Dakota_Regular_edges.csv"), row.names = FALSE)
write.csv(as.data.frame(vertex_attr(g_sd_r)), file.path(output_dir, "South_Dakota_Regular_nodes.csv"), row.names = FALSE)

## Network metrics for all four treatments

metrics_ks_e <- tibble(Network = "Kansas Early", Nodes = vcount(g_ks_e), Edges = ecount(g_ks_e), Average_degree = mean(degree(g_ks_e)), Density = edge_density(g_ks_e, loops = FALSE), Modularity = modularity(cluster_louvain(g_ks_e)), Clustering = transitivity(g_ks_e, type = "average", isolates = "zero"), Average_path_length = mean_distance(g_ks_e, directed = FALSE, unconnected = TRUE), Components = components(g_ks_e)$no, Positive_edges = sum(E(g_ks_e)$edge_sign == "Positive"), Negative_edges = sum(E(g_ks_e)$edge_sign == "Negative"))

metrics_ks_r <- tibble(Network = "Kansas Regular", Nodes = vcount(g_ks_r), Edges = ecount(g_ks_r), Average_degree = mean(degree(g_ks_r)), Density = edge_density(g_ks_r, loops = FALSE), Modularity = modularity(cluster_louvain(g_ks_r)), Clustering = transitivity(g_ks_r, type = "average", isolates = "zero"), Average_path_length = mean_distance(g_ks_r, directed = FALSE, unconnected = TRUE), Components = components(g_ks_r)$no, Positive_edges = sum(E(g_ks_r)$edge_sign == "Positive"), Negative_edges = sum(E(g_ks_r)$edge_sign == "Negative"))

metrics_sd_e <- tibble(Network = "South Dakota Early", Nodes = vcount(g_sd_e), Edges = ecount(g_sd_e), Average_degree = mean(degree(g_sd_e)), Density = edge_density(g_sd_e, loops = FALSE), Modularity = modularity(cluster_louvain(g_sd_e)), Clustering = transitivity(g_sd_e, type = "average", isolates = "zero"), Average_path_length = mean_distance(g_sd_e, directed = FALSE, unconnected = TRUE), Components = components(g_sd_e)$no, Positive_edges = sum(E(g_sd_e)$edge_sign == "Positive"), Negative_edges = sum(E(g_sd_e)$edge_sign == "Negative"))

metrics_sd_r <- tibble(Network = "South Dakota Regular", Nodes = vcount(g_sd_r), Edges = ecount(g_sd_r), Average_degree = mean(degree(g_sd_r)), Density = edge_density(g_sd_r, loops = FALSE), Modularity = modularity(cluster_louvain(g_sd_r)), Clustering = transitivity(g_sd_r, type = "average", isolates = "zero"), Average_path_length = mean_distance(g_sd_r, directed = FALSE, unconnected = TRUE), Components = components(g_sd_r)$no, Positive_edges = sum(E(g_sd_r)$edge_sign == "Positive"), Negative_edges = sum(E(g_sd_r)$edge_sign == "Negative"))

network_metrics <- bind_rows(metrics_ks_e, metrics_ks_r, metrics_sd_e, metrics_sd_r)
network_metrics
write.csv(network_metrics, file.path(output_dir, "All_network_metrics.csv"), row.names = FALSE)

### Hub ASVs based on degree

hub_ks_e <- as_tibble(as_data_frame(g_ks_e, what = "vertices")) %>% arrange(desc(degree))
hub_ks_r <- as_tibble(as_data_frame(g_ks_r, what = "vertices")) %>% arrange(desc(degree))
hub_sd_e <- as_tibble(as_data_frame(g_sd_e, what = "vertices")) %>% arrange(desc(degree))
hub_sd_r <- as_tibble(as_data_frame(g_sd_r, what = "vertices")) %>% arrange(desc(degree))

write.csv(hub_ks_e, file.path(output_dir, "Kansas_Early_hub_ASVs.csv"), row.names = FALSE)
write.csv(hub_ks_r, file.path(output_dir, "Kansas_Regular_hub_ASVs.csv"), row.names = FALSE)
write.csv(hub_sd_e, file.path(output_dir, "South_Dakota_Early_hub_ASVs.csv"), row.names = FALSE)
write.csv(hub_sd_r, file.path(output_dir, "South_Dakota_Regular_hub_ASVs.csv"), row.names = FALSE)

#### Gatekeepers based on betweenness

gatekeeper_ks_e <- as_tibble(as_data_frame(g_ks_e, what = "vertices")) %>% arrange(desc(betweenness))
gatekeeper_ks_r <- as_tibble(as_data_frame(g_ks_r, what = "vertices")) %>% arrange(desc(betweenness))
gatekeeper_sd_e <- as_tibble(as_data_frame(g_sd_e, what = "vertices")) %>% arrange(desc(betweenness))
gatekeeper_sd_r <- as_tibble(as_data_frame(g_sd_r, what = "vertices")) %>% arrange(desc(betweenness))

write.csv(gatekeeper_ks_e, file.path(output_dir, "Kansas_Early_gatekeeper_ASVs.csv"), row.names = FALSE)
write.csv(gatekeeper_ks_r, file.path(output_dir, "Kansas_Regular_gatekeeper_ASVs.csv"), row.names = FALSE)
write.csv(gatekeeper_sd_e, file.path(output_dir, "South_Dakota_Early_gatekeeper_ASVs.csv"), row.names = FALSE)
write.csv(gatekeeper_sd_r, file.path(output_dir, "South_Dakota_Regular_gatekeeper_ASVs.csv"), row.names = FALSE)


