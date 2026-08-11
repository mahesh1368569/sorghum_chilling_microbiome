#MBRACE Coastal Project_Nov2025 data analysis

#Setting working directory
setwd("C:/MBRACE_Nov_Network Analysis")

list.files()

project_dir <- getwd()


#Fixing file names

otu_file  <- file.path(project_dir, "merged_otutable.txt")
meta_file <- file.path(project_dir, "metadata.txt")
tax_file  <- file.path(project_dir, "taxonomy.csv.tsv")
env_file  <- file.path(project_dir, "Environmental metadata_Nov2025.xlsx")


#Creating output folder

output_dir <- file.path(project_dir, "network_outputs_all_taxa")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)


#Checking input files

file.exists(otu_file)
file.exists(meta_file)
file.exists(tax_file)
file.exists(env_file)


#Installing and loading packages

install.packages("tidyverse")
install.packages("readxl")
install.packages("janitor")
install.packages("Hmisc")
install.packages("igraph")
install.packages("ggraph")
install.packages("tidygraph")
install.packages("patchwork")
install.packages("scales")
install.packages("ggrepel")

library(tidyverse)
library(readxl)
library(janitor)
library(Hmisc)
library(igraph)
library(ggraph)
library(tidygraph)
library(patchwork)
library(scales)
library(ggrepel)


#Setting network analysis options

#I am keeping all OTUs in the analysis.
#I am not removing rare taxa based on abundance or prevalence.
#Instead, I am collapsing all OTUs to a selected taxonomic rank before building the network.
#This keeps rare OTUs in the dataset because their counts are still added into the selected taxonomic group.

network_rank <- "Phylum"
color_rank   <- "Phylum"

#These cutoffs are used for microbe-microbe correlations.
#If the network is too dense, I can increase rho to 0.70 or 0.80.
#If the network is too sparse, I can decrease rho to 0.50.

rho_cutoff_taxa_taxa <- 0.60
q_cutoff_taxa_taxa   <- 0.05

#These cutoffs are used for microbe-environment correlations.

rho_cutoff_taxa_env <- 0.50
q_cutoff_taxa_env   <- 0.05

#This is only to avoid making networks from very small groups.
#This is not rare-taxa filtering.

min_samples_for_network <- 6

#This is a safety limit so R does not try to calculate too many pairwise correlations.

max_taxa_for_correlation <- 5000

set.seed(123)


#Reading OTU table

read_otu_table <- function(path){
  
  #Some OTU tables have an extra first line from biom conversion.
  #If that line is present, this function skips it.
  
  first_line <- readLines(path, n = 1)
  
  skip_n <- ifelse(
    grepl("Constructed from biom", first_line, ignore.case = TRUE),
    1,
    0
  )
  
  otu <- read.delim(
    path,
    skip = skip_n,
    check.names = FALSE,
    comment.char = ""
  )
  
  colnames(otu)[1] <- "FeatureID"
  
  return(otu)
}


#Parsing taxonomy

parse_taxonomy <- function(taxonomy_table){
  
  #This separates the taxonomy string into different taxonomic ranks.
  #Unclassified ranks are kept as "Unclassified" instead of NA.
  
  tax <- taxonomy_table %>%
    rename(
      FeatureID = feature_id,
      Taxon = taxon
    ) %>%
    mutate(Taxon = as.character(Taxon)) %>%
    separate(
      Taxon,
      into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"),
      sep = ";",
      fill = "right",
      remove = FALSE
    ) %>%
    mutate(
      across(
        c(Kingdom, Phylum, Class, Order, Family, Genus, Species),
        ~str_trim(.x)
      )
    ) %>%
    mutate(
      across(
        c(Kingdom, Phylum, Class, Order, Family, Genus, Species),
        ~str_remove(.x, "^[a-z]__")
      )
    ) %>%
    mutate(
      across(
        c(Kingdom, Phylum, Class, Order, Family, Genus, Species),
        ~na_if(.x, "")
      )
    ) %>%
    mutate(
      across(
        c(Kingdom, Phylum, Class, Order, Family, Genus, Species),
        ~replace_na(.x, "Unclassified")
      )
    )
  
  return(tax)
}


#Creating taxon labels

make_taxon_label <- function(tax, rank = "Genus"){
  
  #This decides what each OTU will be called after collapsing to the selected rank.
  #Unclassified taxa are grouped carefully so they do not become thousands of separate nodes.
  
  rank <- match.arg(rank, c("FeatureID", "Phylum", "Class", "Order", "Family", "Genus"))
  
  if(rank == "FeatureID"){
    tax$network_taxon <- tax$FeatureID
  }
  
  if(rank == "Phylum"){
    tax$network_taxon <- dplyr::case_when(
      tax$Phylum != "Unclassified" ~ tax$Phylum,
      tax$Kingdom != "Unclassified" ~ paste0(tax$Kingdom, "_unclassified_phylum"),
      TRUE ~ "Unclassified_Phylum"
    )
  }
  
  if(rank == "Class"){
    tax$network_taxon <- dplyr::case_when(
      tax$Class  != "Unclassified" ~ tax$Class,
      tax$Phylum != "Unclassified" ~ paste0(tax$Phylum, "_unclassified_class"),
      tax$Kingdom != "Unclassified" ~ paste0(tax$Kingdom, "_unclassified_class"),
      TRUE ~ "Unclassified_Class"
    )
  }
  
  if(rank == "Order"){
    tax$network_taxon <- dplyr::case_when(
      tax$Order  != "Unclassified" ~ tax$Order,
      tax$Class  != "Unclassified" ~ paste0(tax$Class, "_unclassified_order"),
      tax$Phylum != "Unclassified" ~ paste0(tax$Phylum, "_unclassified_order"),
      tax$Kingdom != "Unclassified" ~ paste0(tax$Kingdom, "_unclassified_order"),
      TRUE ~ "Unclassified_Order"
    )
  }
  
  if(rank == "Family"){
    tax$network_taxon <- dplyr::case_when(
      tax$Family != "Unclassified" ~ tax$Family,
      tax$Order  != "Unclassified" ~ paste0(tax$Order, "_unclassified_family"),
      tax$Class  != "Unclassified" ~ paste0(tax$Class, "_unclassified_family"),
      tax$Phylum != "Unclassified" ~ paste0(tax$Phylum, "_unclassified_family"),
      tax$Kingdom != "Unclassified" ~ paste0(tax$Kingdom, "_unclassified_family"),
      TRUE ~ "Unclassified_Family"
    )
  }
  
  if(rank == "Genus"){
    tax$network_taxon <- dplyr::case_when(
      tax$Genus  != "Unclassified" ~ tax$Genus,
      tax$Family != "Unclassified" ~ paste0(tax$Family, "_unclassified_genus"),
      tax$Order  != "Unclassified" ~ paste0(tax$Order, "_unclassified_genus"),
      tax$Class  != "Unclassified" ~ paste0(tax$Class, "_unclassified_genus"),
      tax$Phylum != "Unclassified" ~ paste0(tax$Phylum, "_unclassified_genus"),
      tax$Kingdom != "Unclassified" ~ paste0(tax$Kingdom, "_unclassified_genus"),
      TRUE ~ "Unclassified_Genus"
    )
  }
  
  return(tax)
}


#Aggregating OTU table to selected taxonomic rank

aggregate_otu_to_taxon <- function(otu_table, taxonomy, sample_ids, rank = "Genus"){
  
  #This collapses OTUs to the selected taxonomic rank.
  #No OTUs are removed here.
  #All counts are included in the summed abundance.
  
  tax <- make_taxon_label(taxonomy, rank = rank)
  
  otu_long <- otu_table %>%
    select(FeatureID, all_of(sample_ids)) %>%
    pivot_longer(
      cols = all_of(sample_ids),
      names_to = "sampleid",
      values_to = "count"
    ) %>%
    left_join(
      tax %>%
        select(FeatureID, network_taxon, Kingdom, Phylum, Class, Order, Family, Genus),
      by = "FeatureID"
    )
  
  taxon_counts <- otu_long %>%
    group_by(network_taxon, sampleid) %>%
    summarise(count = sum(count, na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(
      names_from = sampleid,
      values_from = count,
      values_fill = 0
    )
  
  taxon_taxonomy <- otu_long %>%
    group_by(network_taxon) %>%
    summarise(
      Kingdom = first(na.omit(Kingdom)),
      Phylum  = first(na.omit(Phylum)),
      Class   = first(na.omit(Class)),
      Order   = first(na.omit(Order)),
      Family  = first(na.omit(Family)),
      Genus   = first(na.omit(Genus)),
      n_features_collapsed = n_distinct(FeatureID),
      .groups = "drop"
    )
  
  count_mat <- taxon_counts %>%
    column_to_rownames("network_taxon") %>%
    as.matrix()
  
  return(
    list(
      count_mat = count_mat,
      taxonomy = taxon_taxonomy
    )
  )
}


#Making relative abundance table

make_relative_abundance <- function(count_mat){
  
  sample_sums <- colSums(count_mat, na.rm = TRUE)
  
  rel_mat <- sweep(count_mat, 2, sample_sums, FUN = "/")
  rel_mat[is.na(rel_mat)] <- 0
  
  return(rel_mat)
}


#Building microbe-microbe network

build_taxa_network <- function(rel_mat,
                               taxon_taxonomy,
                               sample_keep,
                               rho_cutoff = 0.6,
                               q_cutoff = 0.05,
                               network_name = "network"){
  
  message("Building network: ", network_name)
  
  mat <- rel_mat[, sample_keep, drop = FALSE]
  
  #This is not rare-taxa filtering.
  #Taxa with zero abundance or no variation in this subset cannot be correlated, 
  #so they are removed only for this step.
  
  keep_taxa <- rowSums(mat, na.rm = TRUE) > 0 &
    apply(mat, 1, sd, na.rm = TRUE) > 0
  
  mat <- mat[keep_taxa, , drop = FALSE]
  
  if(ncol(mat) < min_samples_for_network){
    warning("Network ", network_name, " has too few samples.")
    return(NULL)
  }
  
  if(nrow(mat) < 2){
    warning("Network ", network_name, " has fewer than 2 variable taxa.")
    return(NULL)
  }
  
  if(nrow(mat) > max_taxa_for_correlation){
    stop(
      "Too many taxa for pairwise correlation: ", nrow(mat), " taxa.\n",
      "Use a broader network_rank such as Phylum or Class.\n",
      "Or increase max_taxa_for_correlation if your computer can handle it."
    )
  }
  
  rc <- Hmisc::rcorr(t(mat), type = "spearman")
  
  r_mat <- rc$r
  p_mat <- rc$P
  
  upper_idx <- which(upper.tri(r_mat), arr.ind = TRUE)
  
  edge_tbl <- tibble(
    from = rownames(r_mat)[upper_idx[, 1]],
    to   = colnames(r_mat)[upper_idx[, 2]],
    rho  = r_mat[upper_idx],
    p    = p_mat[upper_idx]
  ) %>%
    filter(!is.na(rho), !is.na(p)) %>%
    mutate(q = p.adjust(p, method = "BH")) %>%
    filter(abs(rho) >= rho_cutoff, q <= q_cutoff) %>%
    mutate(
      edge_sign = ifelse(rho > 0, "Positive", "Negative"),
      network = network_name
    )
  
  if(nrow(edge_tbl) == 0){
    warning("No edges passed thresholds for ", network_name)
    return(NULL)
  }
  
  node_names <- sort(unique(c(edge_tbl$from, edge_tbl$to)))
  
  node_tbl <- tibble(network_taxon = node_names) %>%
    left_join(taxon_taxonomy, by = "network_taxon") %>%
    mutate(
      mean_relative_abundance = rowMeans(mat[network_taxon, , drop = FALSE], na.rm = TRUE),
      network = network_name
    )
  
  g <- graph_from_data_frame(edge_tbl, directed = FALSE, vertices = node_tbl)
  
  V(g)$degree <- degree(g)
  V(g)$betweenness <- betweenness(g, directed = FALSE, normalized = FALSE)
  V(g)$normalized_degree <- degree(g) / max(1, vcount(g) - 1)
  
  if(ecount(g) > 0){
    V(g)$module <- membership(cluster_louvain(g))
  } else {
    V(g)$module <- NA
  }
  
  return(
    list(
      graph = g,
      edges = edge_tbl,
      nodes = as_tibble(as_data_frame(g, what = "vertices")),
      mat = mat
    )
  )
}


#Calculating network metrics

network_metrics <- function(g, network_name = "network"){
  
  if(is.null(g)) return(NULL)
  
  n_nodes <- vcount(g)
  n_edges <- ecount(g)
  
  comps <- components(g)$no
  
  fragmentation <- ifelse(
    n_nodes > 1,
    log(comps) / log(n_nodes),
    NA
  )
  
  edge_signs <- E(g)$edge_sign
  
  out <- tibble(
    network = network_name,
    nodes = n_nodes,
    edges = n_edges,
    average_degree = mean(degree(g)),
    network_diameter = diameter(g, directed = FALSE, unconnected = TRUE),
    graph_density = edge_density(g, loops = FALSE),
    modularity = modularity(cluster_louvain(g)),
    components = comps,
    average_clustering_coefficient = transitivity(g, type = "average", isolates = "zero"),
    average_path_length = mean_distance(g, directed = FALSE, unconnected = TRUE),
    negative_edges = sum(edge_signs == "Negative", na.rm = TRUE),
    percent_negative_edges = 100 * sum(edge_signs == "Negative", na.rm = TRUE) / n_edges,
    fragmentation = fragmentation
  )
  
  return(out)
}


#Plotting microbe-microbe network

plot_network <- function(g,
                         size_var = "degree",
                         color_var = "Phylum",
                         title = "Network",
                         layout_df = NULL){
  
  if(is.null(g)){
    return(
      ggplot() +
        theme_void() +
        ggtitle(paste(title, "(no network)"))
    )
  }
  
  if(is.null(layout_df)){
    layout_df <- create_layout(g, layout = "fr")
  }
  
  p <- ggraph(layout_df) +
    geom_edge_link(
      aes(color = edge_sign),
      alpha = 0.25,
      linewidth = 0.25
    ) +
    geom_node_point(
      aes(size = .data[[size_var]], color = .data[[color_var]]),
      alpha = 0.95
    ) +
    scale_edge_color_manual(
      values = c("Positive" = "grey40", "Negative" = "firebrick")
    ) +
    scale_size_continuous(range = c(1.5, 9)) +
    guides(
      edge_color = guide_legend(title = "Association"),
      size = guide_legend(title = size_var),
      color = guide_legend(title = color_var, override.aes = list(size = 5))
    ) +
    theme_void() +
    ggtitle(title) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 14, color = "black"),
      legend.position = "right",
      legend.title = element_text(color = "black"),
      legend.text = element_text(color = "black"),
      legend.background = element_rect(fill = "white", color = NA),
      legend.key = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA)
    )
  
  return(p)
}


#Calculating fragmentation after node removal

fragmentation_curve <- function(g,
                                remove_by = "betweenness",
                                network_name = "network",
                                n_remove = 10){
  
  if(is.null(g)) return(NULL)
  
  node_order <- as_tibble(as_data_frame(g, what = "vertices")) %>%
    arrange(desc(.data[[remove_by]])) %>%
    pull(name)
  
  n_remove <- min(n_remove, length(node_order) - 1)
  
  out <- map_dfr(0:n_remove, function(i){
    
    if(i == 0){
      g_i <- g
    } else {
      g_i <- delete_vertices(g, node_order[1:i])
    }
    
    n_nodes <- vcount(g_i)
    comps <- components(g_i)$no
    
    frag <- ifelse(
      n_nodes > 1,
      log(comps) / log(n_nodes),
      NA
    )
    
    tibble(
      network = network_name,
      remove_by = remove_by,
      removed_nodes = i,
      remaining_nodes = n_nodes,
      components = comps,
      fragmentation = frag
    )
  })
  
  return(out)
}


#Making degree and betweenness table

node_degree_betweenness_tbl <- function(g, network_name = "network"){
  
  if(is.null(g)) return(NULL)
  
  out <- as_tibble(as_data_frame(g, what = "vertices")) %>%
    mutate(
      network = network_name,
      degree_plus = degree + 1,
      betweenness_plus = betweenness + 1,
      log10_normalized_degree = log10(normalized_degree + 1e-6),
      log10_betweenness = log10(betweenness_plus)
    )
  
  return(out)
}


#Building microbe-environment network

build_taxa_environment_network <- function(rel_mat,
                                           taxon_taxonomy,
                                           sample_meta_env,
                                           env_vars,
                                           sample_keep,
                                           rho_cutoff = 0.5,
                                           q_cutoff = 0.05,
                                           network_name = "taxa_environment"){
  
  message("Building taxa-environment network: ", network_name)
  
  mat <- rel_mat[, sample_keep, drop = FALSE]
  
  keep_taxa <- rowSums(mat, na.rm = TRUE) > 0 &
    apply(mat, 1, sd, na.rm = TRUE) > 0
  
  mat <- mat[keep_taxa, , drop = FALSE]
  
  if(nrow(mat) < 2){
    warning("Microbe-environment network has fewer than 2 variable taxa.")
    return(NULL)
  }
  
  env_df <- sample_meta_env %>%
    filter(sampleid %in% sample_keep) %>%
    arrange(match(sampleid, sample_keep)) %>%
    select(sampleid, all_of(env_vars))
  
  env_vars_use <- env_vars[
    sapply(env_df[env_vars], function(x) sd(x, na.rm = TRUE) > 0)
  ]
  
  edge_list <- list()
  
  for(env in env_vars_use){
    
    env_values <- env_df[[env]]
    
    tmp <- map_dfr(rownames(mat), function(taxon){
      
      test <- suppressWarnings(
        cor.test(
          as.numeric(mat[taxon, ]),
          env_values,
          method = "spearman",
          exact = FALSE
        )
      )
      
      tibble(
        from = taxon,
        to = env,
        rho = unname(test$estimate),
        p = test$p.value
      )
    })
    
    edge_list[[env]] <- tmp
  }
  
  edge_tbl <- bind_rows(edge_list) %>%
    filter(!is.na(rho), !is.na(p)) %>%
    mutate(q = p.adjust(p, method = "BH")) %>%
    filter(abs(rho) >= rho_cutoff, q <= q_cutoff) %>%
    mutate(
      edge_sign = ifelse(rho > 0, "Positive", "Negative"),
      network = network_name
    )
  
  if(nrow(edge_tbl) == 0){
    warning("No microbe-environment edges passed thresholds.")
    return(NULL)
  }
  
  taxa_nodes <- tibble(network_taxon = unique(edge_tbl$from)) %>%
    left_join(taxon_taxonomy, by = "network_taxon") %>%
    mutate(
      node_type = "Microbial taxon",
      display_group = .data[[color_rank]],
      mean_relative_abundance = rowMeans(mat[network_taxon, , drop = FALSE], na.rm = TRUE)
    ) %>%
    rename(name = network_taxon)
  
  env_nodes <- tibble(
    name = unique(edge_tbl$to),
    Kingdom = "Environment",
    Phylum = "Environment",
    Class = "Environment",
    Order = "Environment",
    Family = "Environment",
    Genus = "Environment",
    node_type = "Environmental factor",
    display_group = "Environmental factor",
    mean_relative_abundance = NA_real_
  )
  
  node_tbl <- bind_rows(taxa_nodes, env_nodes)
  
  g <- graph_from_data_frame(edge_tbl, directed = FALSE, vertices = node_tbl)
  
  V(g)$degree <- degree(g)
  V(g)$betweenness <- betweenness(g, directed = FALSE, normalized = FALSE)
  V(g)$normalized_degree <- degree(g) / max(1, vcount(g) - 1)
  
  return(
    list(
      graph = g,
      edges = edge_tbl,
      nodes = as_tibble(as_data_frame(g, what = "vertices")),
      mat = mat
    )
  )
}


#Plotting microbe-environment network

plot_environment_network <- function(g, title = "Microbe-environment association network"){
  
  if(is.null(g)){
    return(
      ggplot() +
        theme_void() +
        ggtitle("No microbe-environment network")
    )
  }
  
  g_tbl <- as_tbl_graph(g)
  
  p <- ggraph(g_tbl, layout = "fr") +
    geom_edge_link(
      aes(color = edge_sign),
      alpha = 0.25,
      linewidth = 0.25
    ) +
    geom_node_point(
      aes(size = degree, color = display_group, shape = node_type),
      alpha = 0.95
    ) +
    geom_node_text(
      aes(label = ifelse(node_type == "Environmental factor", name, "")),
      repel = TRUE,
      size = 3.5,
      fontface = "bold"
    ) +
    scale_edge_color_manual(
      values = c("Positive" = "grey40", "Negative" = "firebrick")
    ) +
    scale_shape_manual(
      values = c("Microbial taxon" = 16, "Environmental factor" = 18)
    ) +
    scale_size_continuous(range = c(2, 10)) +
    theme_void() +
    ggtitle(title) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
      legend.position = "right"
    )
  
  return(p)
}


#Reading raw files

otu_raw <- read_otu_table(otu_file)

metadata <- read.delim(
  meta_file,
  check.names = FALSE
) %>%
  clean_names()

taxonomy_raw <- read.delim(
  tax_file,
  check.names = FALSE
) %>%
  clean_names()

env_raw <- read_excel(env_file) %>%
  clean_names()


#Checking data after reading

dim(otu_raw)
head(otu_raw[, 1:5])

dim(metadata)
head(metadata)

dim(taxonomy_raw)
head(taxonomy_raw)

dim(env_raw)
head(env_raw)


#Preparing taxonomy

taxonomy <- parse_taxonomy(taxonomy_raw)

head(taxonomy)

test_tax <- make_taxon_label(taxonomy, rank = network_rank)

head(test_tax$network_taxon)
length(unique(test_tax$network_taxon))


#Keeping water samples only

otu_sample_ids <- setdiff(colnames(otu_raw), "FeatureID")

water_sample_ids <- otu_sample_ids[grepl("^W_", otu_sample_ids)]

metadata <- metadata %>%
  rename(sampleid = sampleid)

water_sample_ids <- intersect(water_sample_ids, metadata$sampleid)

sample_info <- metadata %>%
  filter(sampleid %in% water_sample_ids) %>%
  mutate(
    station_id = str_remove(sampleid, "_\\d+$"),
    bay = location,
    salinity_group = salinity
  )


#Matching environmental metadata to sequencing samples

env_clean <- env_raw %>%
  rename(station_id = sampling_site) %>%
  mutate(station_id = as.character(station_id))

sample_info_env <- sample_info %>%
  left_join(env_clean, by = "station_id")

message("Number of water sequencing samples used: ", length(water_sample_ids))
message("Number of water samples with environmental metadata: ", sum(!is.na(sample_info_env$salinity_ppt)))

write.csv(
  sample_info_env,
  file.path(output_dir, "matched_water_sample_metadata.csv"),
  row.names = FALSE
)


#Aggregating OTUs and making relative abundance table

agg <- aggregate_otu_to_taxon(
  otu_table = otu_raw,
  taxonomy = taxonomy,
  sample_ids = water_sample_ids,
  rank = network_rank
)

count_mat <- agg$count_mat
taxon_taxonomy <- agg$taxonomy

rel_mat <- make_relative_abundance(count_mat)

network_rank
nrow(rel_mat)
head(rownames(rel_mat), 20)


#Figure 1 overall network

overall_net <- build_taxa_network(
  rel_mat = rel_mat,
  taxon_taxonomy = taxon_taxonomy,
  sample_keep = water_sample_ids,
  rho_cutoff = rho_cutoff_taxa_taxa,
  q_cutoff = q_cutoff_taxa_taxa,
  network_name = "Overall water network"
)

is.null(overall_net)

if(!is.null(overall_net)){
  
  overall_g <- overall_net$graph
  
  vcount(overall_g)
  ecount(overall_g)
  
  write.csv(
    overall_net$edges,
    file.path(output_dir, "Figure1_overall_network_edges.csv"),
    row.names = FALSE
  )
  
  write.csv(
    as_tibble(as_data_frame(overall_g, what = "vertices")),
    file.path(output_dir, "Figure1_overall_network_nodes.csv"),
    row.names = FALSE
  )
  
  write.csv(
    network_metrics(overall_g, "Overall water network"),
    file.path(output_dir, "Figure1_overall_network_metrics.csv"),
    row.names = FALSE
  )
  
  overall_layout <- create_layout(overall_g, layout = "fr")
  
  p1a <- plot_network(
    overall_g,
    size_var = "mean_relative_abundance",
    color_var = color_rank,
    title = "A. Abundance",
    layout_df = overall_layout
  )
  
  print(p1a)
  
  ggsave(
    file.path(output_dir, "Figure1A_overall_network_abundance.png"),
    p1a,
    width = 14,
    height = 10,
    dpi = 300,
    bg = "white"
  )
  
  p1b <- plot_network(
    overall_g,
    size_var = "degree",
    color_var = color_rank,
    title = "B. Degree",
    layout_df = overall_layout
  )
  
  print(p1b)
  
  ggsave(
    file.path(output_dir, "Figure1B_overall_network_degree.png"),
    p1b,
    width = 14,
    height = 10,
    dpi = 300,
    bg = "white"
  )
  
  p1c <- plot_network(
    overall_g,
    size_var = "betweenness",
    color_var = color_rank,
    title = "C. Betweenness",
    layout_df = overall_layout
  )
  
  print(p1c)
  
  ggsave(
    file.path(output_dir, "Figure1C_overall_network_betweenness.png"),
    p1c,
    width = 14,
    height = 10,
    dpi = 300,
    bg = "white"
  )
}

list.files(output_dir)


#Figure 2 bay-specific networks

bay_names <- sample_info_env %>%
  distinct(bay) %>%
  pull(bay) %>%
  na.omit()

bay_networks <- list()
bay_plots <- list()
bay_metrics <- list()

for(bay_i in bay_names){
  
  samples_i <- sample_info_env %>%
    filter(bay == bay_i) %>%
    pull(sampleid)
  
  if(length(samples_i) >= min_samples_for_network){
    
    net_i <- build_taxa_network(
      rel_mat = rel_mat,
      taxon_taxonomy = taxon_taxonomy,
      sample_keep = samples_i,
      rho_cutoff = rho_cutoff_taxa_taxa,
      q_cutoff = q_cutoff_taxa_taxa,
      network_name = bay_i
    )
    
    bay_networks[[bay_i]] <- net_i
    
    if(!is.null(net_i)){
      
      bay_metrics[[bay_i]] <- network_metrics(net_i$graph, bay_i)
      
      bay_plots[[bay_i]] <- plot_network(
        net_i$graph,
        size_var = "degree",
        color_var = color_rank,
        title = bay_i
      )
      
      print(bay_plots[[bay_i]])
      
      ggsave(
        file.path(output_dir, paste0("Figure2_", make_clean_names(bay_i), "_network.png")),
        bay_plots[[bay_i]],
        width = 12,
        height = 9,
        dpi = 300,
        bg = "white"
      )
      
      write.csv(
        net_i$edges,
        file.path(output_dir, paste0("Figure2_", make_clean_names(bay_i), "_edges.csv")),
        row.names = FALSE
      )
      
      write.csv(
        as_tibble(as_data_frame(net_i$graph, what = "vertices")),
        file.path(output_dir, paste0("Figure2_", make_clean_names(bay_i), "_nodes.csv")),
        row.names = FALSE
      )
    }
  }
}

if(length(bay_plots) > 0){
  
  fig2 <- wrap_plots(bay_plots, ncol = 2, guides = "collect") +
    plot_annotation(title = "Biloxi Bay vs Pascagoula Bay co-occurrence networks") &
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 18, color = "black"),
      legend.title = element_text(size = 12, face = "bold", color = "black"),
      legend.text = element_text(size = 10, color = "black"),
      legend.background = element_rect(fill = "white", color = NA),
      legend.key = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA)
    )
  
  print(fig2)
  
  ggsave(
    file.path(output_dir, "Figure2_Biloxi_vs_Pascagoula_networks.png"),
    fig2,
    width = 16,
    height = 8,
    dpi = 300,
    bg = "white"
  )
  
  bind_rows(bay_metrics) %>%
    write.csv(
      file.path(output_dir, "Figure2_bay_network_metrics.csv"),
      row.names = FALSE
    )
}


#Figure 3 network fragmentation

fragment_list <- list()

if(!is.null(overall_net)){
  
  fragment_list[["Overall"]] <- fragmentation_curve(
    overall_net$graph,
    remove_by = "betweenness",
    network_name = "Overall"
  )
}

for(nm in names(bay_networks)){
  
  if(!is.null(bay_networks[[nm]])){
    
    fragment_list[[nm]] <- fragmentation_curve(
      bay_networks[[nm]]$graph,
      remove_by = "betweenness",
      network_name = nm
    )
  }
}

fragment_df <- bind_rows(fragment_list)

if(nrow(fragment_df) > 0){
  
  write.csv(
    fragment_df,
    file.path(output_dir, "Figure3_fragmentation_values.csv"),
    row.names = FALSE
  )
  
  fig3 <- ggplot(
    fragment_df,
    aes(
      x = removed_nodes,
      y = fragmentation,
      group = network,
      shape = network
    )
  ) +
    geom_line(
      linewidth = 0.8,
      color = "black"
    ) +
    geom_point(
      size = 3.0,
      color = "black"
    ) +
    scale_x_continuous(
      breaks = 0:10,
      limits = c(0, 10)
    ) +
    labs(
      x = "Number of removed high-betweenness nodes",
      y = "Fragmentation",
      title = "Network fragmentation after removing high-betweenness nodes"
    ) +
    theme_classic(base_size = 13) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, color = "black"),
      axis.title = element_text(color = "black"),
      axis.text = element_text(color = "black"),
      legend.title = element_text(color = "black", face = "bold"),
      legend.text = element_text(color = "black"),
      legend.background = element_rect(fill = "white", color = NA),
      legend.key = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA)
    )
  
  print(fig3)
  
  ggsave(
    file.path(output_dir, "Figure3_fragmentation_high_betweenness_removal.png"),
    fig3,
    width = 8,
    height = 6,
    dpi = 300,
    bg = "white"
  )
}


#Figure 4 normalized degree and betweenness

degree_bet_list <- list()

if(!is.null(overall_net)){
  
  degree_bet_list[["Overall"]] <- node_degree_betweenness_tbl(
    overall_net$graph,
    "Overall"
  )
}

for(nm in names(bay_networks)){
  
  if(!is.null(bay_networks[[nm]])){
    
    degree_bet_list[[nm]] <- node_degree_betweenness_tbl(
      bay_networks[[nm]]$graph,
      nm
    )
  }
}

degree_bet_df <- bind_rows(degree_bet_list)

if(nrow(degree_bet_df) > 0){
  
  write.csv(
    degree_bet_df,
    file.path(output_dir, "Figure4_degree_betweenness_values.csv"),
    row.names = FALSE
  )
  
  corr_labels <- degree_bet_df %>%
    group_by(network) %>%
    summarise(
      rho = suppressWarnings(
        cor(
          log10_normalized_degree,
          log10_betweenness,
          method = "spearman",
          use = "complete.obs"
        )
      ),
      p_value = suppressWarnings(
        cor.test(
          log10_normalized_degree,
          log10_betweenness,
          method = "spearman",
          exact = FALSE
        )$p.value
      ),
      .groups = "drop"
    ) %>%
    mutate(
      label = paste0(
        "rho = ", round(rho, 3),
        "\nP = ", signif(p_value, 2)
      )
    )
  
  fig4 <- ggplot(
    degree_bet_df,
    aes(x = log10_normalized_degree, y = log10_betweenness)
  ) +
    geom_point(
      color = "black",
      size = 1.8,
      alpha = 0.75
    ) +
    facet_wrap(~network, scales = "free_y", nrow = 1) +
    geom_text(
      data = corr_labels,
      aes(x = Inf, y = Inf, label = label),
      inherit.aes = FALSE,
      hjust = 1.1,
      vjust = 1.3,
      size = 4,
      color = "black"
    ) +
    scale_y_continuous(
      expand = expansion(mult = c(0.05, 0.22))
    ) +
    labs(
      x = "log10 normalized degree",
      y = "log10 betweenness",
      title = "Correlation between normalized degree and betweenness centrality"
    ) +
    theme_classic(base_size = 13) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, color = "black"),
      strip.background = element_rect(fill = "white", color = "black"),
      strip.text = element_text(face = "bold", color = "black", size = 12),
      axis.title = element_text(color = "black"),
      axis.text = element_text(color = "black"),
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA)
    )
  
  print(fig4)
  
  ggsave(
    file.path(output_dir, "Figure4_normalized_degree_vs_betweenness_with_pvalue_free_y.png"),
    fig4,
    width = 12,
    height = 4.5,
    dpi = 300,
    bg = "white"
  )
}


#Figure 5 bay-specific combined microbe and environment networks

candidate_env_vars <- c(
  "temperature",
  "pressure_mmhg",
  "dissolved_oxygen_percent",
  "dissolved_oxygen_mg_l",
  "specific_conductance_us_cm",
  "salinity_ppt"
)

env_vars <- intersect(candidate_env_vars, colnames(sample_info_env))

sample_info_env <- sample_info_env %>%
  mutate(across(all_of(env_vars), ~as.numeric(.x)))


#Making one combined microbe-environment network for each bay

make_bay_combined_microbe_environment_network <- function(
    bay_i,
    bay_net,
    sample_info_env,
    rel_mat,
    taxon_taxonomy,
    env_vars,
    output_dir
){
  
  message("Creating combined microbe-environment network for: ", bay_i)
  
  samples_i <- sample_info_env %>%
    filter(bay == bay_i) %>%
    pull(sampleid)
  
  if(is.null(bay_net)){
    message("Skipping ", bay_i, ": bay microbial network is NULL.")
    return(NULL)
  }
  
  env_net_i <- build_taxa_environment_network(
    rel_mat = rel_mat,
    taxon_taxonomy = taxon_taxonomy,
    sample_meta_env = sample_info_env,
    env_vars = env_vars,
    sample_keep = samples_i,
    rho_cutoff = rho_cutoff_taxa_env,
    q_cutoff = q_cutoff_taxa_env,
    network_name = paste0(bay_i, " microbe-environment network")
  )
  
  if(is.null(env_net_i)){
    message("Skipping ", bay_i, ": no microbe-environment edges passed thresholds.")
    return(NULL)
  }
  
  microbe_edges <- bay_net$edges %>%
    select(from, to, rho, p, q, edge_sign) %>%
    mutate(
      edge_type = "Microbe-microbe"
    )
  
  env_edges <- env_net_i$edges %>%
    select(from, to, rho, p, q, edge_sign) %>%
    mutate(
      edge_type = "Microbe-environment"
    )
  
  combined_edges <- bind_rows(microbe_edges, env_edges)
  
  all_node_names <- sort(unique(c(combined_edges$from, combined_edges$to)))
  
  environmental_node_names <- env_vars
  
  microbial_node_names <- setdiff(all_node_names, environmental_node_names)
  
  microbial_nodes <- tibble(network_taxon = microbial_node_names) %>%
    left_join(taxon_taxonomy, by = "network_taxon") %>%
    mutate(
      name = network_taxon,
      node_type = "Microbial taxon",
      display_group = .data[[color_rank]]
    ) %>%
    select(
      name,
      node_type,
      display_group,
      Kingdom,
      Phylum,
      Class,
      Order,
      Family,
      Genus,
      n_features_collapsed
    )
  
  environmental_nodes <- tibble(
    name = environmental_node_names,
    node_type = "Environmental factor",
    display_group = "Environmental factor",
    Kingdom = "Environment",
    Phylum = "Environment",
    Class = "Environment",
    Order = "Environment",
    Family = "Environment",
    Genus = "Environment",
    n_features_collapsed = NA_integer_
  )
  
  combined_nodes <- bind_rows(microbial_nodes, environmental_nodes) %>%
    filter(name %in% all_node_names)
  
  combined_g <- graph_from_data_frame(
    combined_edges,
    directed = FALSE,
    vertices = combined_nodes
  )
  
  V(combined_g)$degree <- degree(combined_g)
  V(combined_g)$betweenness <- betweenness(combined_g, directed = FALSE, normalized = FALSE)
  
  set.seed(123)
  
  base_layout <- layout_with_fr(combined_g)
  
  layout_df <- as_tibble(base_layout) %>%
    rename(x = V1, y = V2) %>%
    mutate(
      name = V(combined_g)$name,
      node_type = V(combined_g)$node_type,
      display_group = V(combined_g)$display_group,
      degree = V(combined_g)$degree
    )
  
  env_layout <- layout_df %>%
    filter(node_type == "Environmental factor") %>%
    arrange(name) %>%
    mutate(
      x = -2.8,
      y = seq(1.7, -1.7, length.out = n())
    )
  
  microbe_layout <- layout_df %>%
    filter(node_type == "Microbial taxon") %>%
    mutate(
      x = scales::rescale(x, to = c(-0.5, 2.8)),
      y = scales::rescale(y, to = c(-2.0, 2.0))
    )
  
  final_layout <- bind_rows(env_layout, microbe_layout) %>%
    arrange(match(name, V(combined_g)$name))
  
  color_groups <- sort(unique(final_layout$display_group))
  color_values <- scales::hue_pal()(length(color_groups))
  names(color_values) <- color_groups
  
  if("Environmental factor" %in% names(color_values)){
    color_values["Environmental factor"] <- "gold"
  }
  
  fig_i <- ggraph(
    combined_g,
    layout = "manual",
    x = final_layout$x,
    y = final_layout$y
  ) +
    geom_edge_link(
      aes(color = edge_sign, linetype = edge_type),
      alpha = 0.25,
      linewidth = 0.25
    ) +
    geom_node_point(
      aes(size = degree, color = display_group, shape = node_type),
      alpha = 0.95
    ) +
    geom_node_text(
      aes(label = ifelse(node_type == "Environmental factor", name, "")),
      repel = TRUE,
      size = 4,
      fontface = "bold",
      color = "black"
    ) +
    scale_edge_color_manual(
      values = c("Positive" = "grey45", "Negative" = "firebrick")
    ) +
    scale_edge_linetype_manual(
      values = c("Microbe-microbe" = "solid", "Microbe-environment" = "solid")
    ) +
    scale_color_manual(
      values = color_values
    ) +
    scale_shape_manual(
      values = c("Environmental factor" = 18, "Microbial taxon" = 16)
    ) +
    scale_size_continuous(
      range = c(2, 10)
    ) +
    guides(
      edge_color = guide_legend(title = "Association"),
      edge_linetype = guide_legend(title = "Edge type"),
      color = guide_legend(title = color_rank, override.aes = list(size = 5)),
      size = guide_legend(title = "Degree"),
      shape = guide_legend(title = "Node type")
    ) +
    labs(
      title = paste0(bay_i, ": microbe-environment association network")
    ) +
    theme_void() +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 18, color = "black"),
      legend.position = "right",
      legend.title = element_text(size = 12, face = "bold", color = "black"),
      legend.text = element_text(size = 9, color = "black"),
      legend.background = element_rect(fill = "white", color = NA),
      legend.key = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA)
    )
  
  clean_bay_name <- make_clean_names(bay_i)
  
  ggsave(
    file.path(output_dir, paste0("Figure5_", clean_bay_name, "_combined_microbe_environment_network.png")),
    fig_i,
    width = 15,
    height = 10,
    dpi = 300,
    bg = "white"
  )
  
  try(print(fig_i), silent = TRUE)
  
  write.csv(
    combined_edges,
    file.path(output_dir, paste0("Figure5_", clean_bay_name, "_combined_microbe_environment_edges.csv")),
    row.names = FALSE
  )
  
  write.csv(
    combined_nodes,
    file.path(output_dir, paste0("Figure5_", clean_bay_name, "_combined_microbe_environment_nodes.csv")),
    row.names = FALSE
  )
  
  return(
    list(
      graph = combined_g,
      edges = combined_edges,
      nodes = combined_nodes,
      plot = fig_i
    )
  )
}


#Running Figure 5 for each bay

bay_combined_networks <- list()

for(bay_i in names(bay_networks)){
  
  bay_combined_networks[[bay_i]] <- make_bay_combined_microbe_environment_network(
    bay_i = bay_i,
    bay_net = bay_networks[[bay_i]],
    sample_info_env = sample_info_env,
    rel_mat = rel_mat,
    taxon_taxonomy = taxon_taxonomy,
    env_vars = env_vars,
    output_dir = output_dir
  )
}

list.files(output_dir, pattern = "Figure5", full.names = FALSE)


#Checking Figure 5 environmental edges

biloxi_edges <- read.csv(
  file.path(output_dir, "Figure5_biloxi_bay_combined_microbe_environment_edges.csv")
)

pascagoula_edges <- read.csv(
  file.path(output_dir, "Figure5_pascagoula_bay_combined_microbe_environment_edges.csv")
)

biloxi_edges %>%
  filter(edge_type == "Microbe-environment") %>%
  count(to) %>%
  arrange(desc(n))

pascagoula_edges %>%
  filter(edge_type == "Microbe-environment") %>%
  count(to) %>%
  arrange(desc(n))

biloxi_edges %>%
  filter(edge_type == "Microbe-environment") %>%
  arrange(to, q, desc(abs(rho))) %>%
  select(from, to, rho, p, q, edge_sign)

pascagoula_edges %>%
  filter(edge_type == "Microbe-environment") %>%
  arrange(to, q, desc(abs(rho))) %>%
  select(from, to, rho, p, q, edge_sign)


#Saving Figure 5 environmental edge summaries

biloxi_env_summary <- biloxi_edges %>%
  filter(edge_type == "Microbe-environment") %>%
  count(to) %>%
  arrange(desc(n))

pascagoula_env_summary <- pascagoula_edges %>%
  filter(edge_type == "Microbe-environment") %>%
  count(to) %>%
  arrange(desc(n))

biloxi_env_summary

pascagoula_env_summary

write.csv(
  biloxi_env_summary,
  file.path(output_dir, "Figure5_biloxi_environment_variable_edge_summary.csv"),
  row.names = FALSE
)

write.csv(
  pascagoula_env_summary,
  file.path(output_dir, "Figure5_pascagoula_environment_variable_edge_summary.csv"),
  row.names = FALSE
)

biloxi_env_edges_detailed <- biloxi_edges %>%
  filter(edge_type == "Microbe-environment") %>%
  arrange(to, q, desc(abs(rho))) %>%
  select(from, to, rho, p, q, edge_sign)

pascagoula_env_edges_detailed <- pascagoula_edges %>%
  filter(edge_type == "Microbe-environment") %>%
  arrange(to, q, desc(abs(rho))) %>%
  select(from, to, rho, p, q, edge_sign)

biloxi_env_edges_detailed

pascagoula_env_edges_detailed

write.csv(
  biloxi_env_edges_detailed,
  file.path(output_dir, "Figure5_biloxi_microbe_environment_edges_detailed.csv"),
  row.names = FALSE
)

write.csv(
  pascagoula_env_edges_detailed,
  file.path(output_dir, "Figure5_pascagoula_microbe_environment_edges_detailed.csv"),
  row.names = FALSE
)

biloxi_env_sign_summary <- biloxi_edges %>%
  filter(edge_type == "Microbe-environment") %>%
  count(to, edge_sign) %>%
  arrange(to, edge_sign)

pascagoula_env_sign_summary <- pascagoula_edges %>%
  filter(edge_type == "Microbe-environment") %>%
  count(to, edge_sign) %>%
  arrange(to, edge_sign)

biloxi_env_sign_summary

pascagoula_env_sign_summary

write.csv(
  biloxi_env_sign_summary,
  file.path(output_dir, "Figure5_biloxi_environment_positive_negative_edge_summary.csv"),
  row.names = FALSE
)

write.csv(
  pascagoula_env_sign_summary,
  file.path(output_dir, "Figure5_pascagoula_environment_positive_negative_edge_summary.csv"),
  row.names = FALSE
)

biloxi_top_env_edges <- biloxi_env_edges_detailed %>%
  arrange(desc(abs(rho))) %>%
  slice_head(n = 20)

pascagoula_top_env_edges <- pascagoula_env_edges_detailed %>%
  arrange(desc(abs(rho))) %>%
  slice_head(n = 20)

biloxi_top_env_edges

pascagoula_top_env_edges

write.csv(
  biloxi_top_env_edges,
  file.path(output_dir, "Figure5_biloxi_top20_microbe_environment_edges.csv"),
  row.names = FALSE
)

write.csv(
  pascagoula_top_env_edges,
  file.path(output_dir, "Figure5_pascagoula_top20_microbe_environment_edges.csv"),
  row.names = FALSE
)

list.files(output_dir, pattern = "Figure5_.*environment.*csv")


#Supplementary salinity-gradient networks

salinity_networks <- list()
salinity_plots <- list()
salinity_metrics <- list()

salinity_groups <- sample_info_env %>%
  filter(!is.na(salinity_group), salinity_group != "") %>%
  distinct(salinity_group) %>%
  pull(salinity_group)

for(sal_i in salinity_groups){
  
  samples_i <- sample_info_env %>%
    filter(salinity_group == sal_i) %>%
    pull(sampleid)
  
  if(length(samples_i) >= min_samples_for_network){
    
    net_i <- build_taxa_network(
      rel_mat = rel_mat,
      taxon_taxonomy = taxon_taxonomy,
      sample_keep = samples_i,
      rho_cutoff = rho_cutoff_taxa_taxa,
      q_cutoff = q_cutoff_taxa_taxa,
      network_name = sal_i
    )
    
    salinity_networks[[sal_i]] <- net_i
    
    if(!is.null(net_i)){
      
      salinity_metrics[[sal_i]] <- network_metrics(net_i$graph, sal_i)
      
      salinity_plots[[sal_i]] <- plot_network(
        net_i$graph,
        size_var = "degree",
        color_var = color_rank,
        title = sal_i
      )
      
      print(salinity_plots[[sal_i]])
      
      ggsave(
        file.path(output_dir, paste0("SuppFig_", make_clean_names(sal_i), "_network.png")),
        salinity_plots[[sal_i]],
        width = 12,
        height = 9,
        dpi = 300,
        bg = "white"
      )
      
      write.csv(
        net_i$edges,
        file.path(output_dir, paste0("SuppFig_", make_clean_names(sal_i), "_edges.csv")),
        row.names = FALSE
      )
      
      write.csv(
        as_tibble(as_data_frame(net_i$graph, what = "vertices")),
        file.path(output_dir, paste0("SuppFig_", make_clean_names(sal_i), "_nodes.csv")),
        row.names = FALSE
      )
    }
    
  } else {
    
    message("Skipping ", sal_i, ": only ", length(samples_i), " samples.")
  }
}

if(length(salinity_plots) > 0){
  
  supp_fig <- wrap_plots(salinity_plots, ncol = 2, guides = "collect") +
    plot_annotation(title = "Salinity-gradient co-occurrence networks") &
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 18, color = "black"),
      legend.title = element_text(size = 12, face = "bold", color = "black"),
      legend.text = element_text(size = 10, color = "black"),
      legend.background = element_rect(fill = "white", color = NA),
      legend.key = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA)
    )
  
  print(supp_fig)
  
  ggsave(
    file.path(output_dir, "Supplementary_Figure_salinity_gradient_networks.png"),
    supp_fig,
    width = 16,
    height = 10,
    dpi = 300,
    bg = "white"
  )
  
  bind_rows(salinity_metrics) %>%
    write.csv(
      file.path(output_dir, "Supplementary_salinity_network_metrics.csv"),
      row.names = FALSE
    )
}


#Saving all network metrics together

all_metrics <- list()

if(!is.null(overall_net)){
  all_metrics[["Overall"]] <- network_metrics(overall_net$graph, "Overall")
}

for(nm in names(bay_networks)){
  
  if(!is.null(bay_networks[[nm]])){
    
    all_metrics[[paste0("Bay_", nm)]] <- network_metrics(
      bay_networks[[nm]]$graph,
      nm
    )
  }
}

for(nm in names(salinity_networks)){
  
  if(!is.null(salinity_networks[[nm]])){
    
    all_metrics[[paste0("Salinity_", nm)]] <- network_metrics(
      salinity_networks[[nm]]$graph,
      nm
    )
  }
}

bind_rows(all_metrics) %>%
  write.csv(
    file.path(output_dir, "All_network_metrics_summary.csv"),
    row.names = FALSE
  )


#Saving hub taxa and gatekeeper taxa tables

if(!is.null(overall_net)){
  
  hub_tbl <- as_tibble(as_data_frame(overall_net$graph, what = "vertices")) %>%
    arrange(desc(degree)) %>%
    select(
      name,
      Kingdom,
      Phylum,
      Class,
      Order,
      Family,
      Genus,
      mean_relative_abundance,
      degree,
      betweenness,
      normalized_degree,
      module,
      n_features_collapsed
    )
  
  gatekeeper_tbl <- as_tibble(as_data_frame(overall_net$graph, what = "vertices")) %>%
    arrange(desc(betweenness)) %>%
    select(
      name,
      Kingdom,
      Phylum,
      Class,
      Order,
      Family,
      Genus,
      mean_relative_abundance,
      degree,
      betweenness,
      normalized_degree,
      module,
      n_features_collapsed
    )
  
  write.csv(
    hub_tbl,
    file.path(output_dir, "Overall_network_hub_taxa_by_degree.csv"),
    row.names = FALSE
  )
  
  write.csv(
    gatekeeper_tbl,
    file.path(output_dir, "Overall_network_gatekeeper_taxa_by_betweenness.csv"),
    row.names = FALSE
  )
}





