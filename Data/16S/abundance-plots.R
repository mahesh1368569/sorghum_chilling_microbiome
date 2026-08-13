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

biom = import_biom("asv_table_16S_filt_taxonomy.biom")

meta <- read.csv("metadata.csv",header = TRUE,stringsAsFactors = FALSE,check.names = FALSE)

metadata = import_qiime_sample_data("metadata.txt")

sample_names(biom)
sample_names(metadata)

sorghum_biom = merge_phyloseq(biom, metadata)

colnames(tax_table(sorghum_biom)) <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")

sorghum <- phyloseq2meco(sorghum_biom)

sorghum

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
abund_phylum_L <- trans_abund$new(dataset = sorghum, taxrank = "Phylum", ntaxa = 10, groupmean = "Line")

its_colorset <- c(
  "Proteobacteria"   = "#3B6FB6",  # muted blue
  "Actinobacteriota" = "#8064A2",  # muted purple
  "Firmicutes"       = "#D95F5F",  # muted coral/red
  "Bacteroidota"     = "#E69F45",  # soft orange
  "Deinococcota"     = "#56B4C8",  # cyan
  "Patescibacteria"  = "#62A66F",  # muted green
  "Armatimonadota"   = "#C77DAA",  # muted magenta
  "Others"           = "#BDBDBD"   # neutral grey
)


taxa_colors <- c(
  "#3B6FB6", "#8064A2", "#D95F5F", "#E69F45", "#56B4C8",
  "#62A66F", "#C77DAA", "#8C6D31", "#4C956C", "#6C5B7B",
  "#E07A5F", "#3D5A80", "#81B29A", "#F2CC8F", "#A44A3F",
  "#748CAB", "#9C6644", "#6D597A", "#B56576", "#5F8D4E",
  "#457B9D", "#BC6C25", "#7A9E9F", "#9B5DE5", "#577590",
  "#F28482", "#84A59D", "#A98467", "#778DA9", "#B08968"
)



p_phylum <- abund_phylum$plot_bar(others_color = "grey70", legend_text_italic = FALSE) +
  theme_classic(base_size = 14) + 
  scale_fill_manual(values = its_colorset) + labs(
    x = NULL, y = "Relative abundance (%)", fill = "Phylum") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 25,hjust = 1,size = 14),
        axis.text.y = element_text(size = 14),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        strip.text = element_text(size = 14),
        legend.position = "right",
        legend.title = element_text(size = 13),
        legend.text = element_text(size = 12),
        panel.border = element_rect(colour = "black",fill = NA, size = 1))

p_phylum

p_phylum_l <- abund_phylum_L$plot_bar(others_color = "grey70", legend_text_italic = FALSE) +
  theme_classic(base_size = 14) + 
  scale_fill_manual(values = its_colorset) + labs(
    x = NULL, y = "Relative abundance (%)", fill = "Phylum") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 25,hjust = 1,size = 14),
        axis.text.y = element_text(size = 14),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        strip.text = element_text(size = 14),
        legend.position = "right",
        legend.title = element_text(size = 13),
        legend.text = element_text(size = 12),
        panel.border = element_rect(colour = "black",fill = NA, size = 1))

p_phylum_l

ggsave("phylum_abundance.pdf",p_phylum ,width = 7 ,height = 7, dpi = 1000)
ggsave("phylum_abundance_line.pdf",p_phylum_l ,width = 7 ,height = 7, dpi = 1000)

abund_genus <- trans_abund$new(dataset = sorghum, taxrank = "Genus", ntaxa = 21, groupmean = "SitePlanting")

p_genus <- abund_genus$plot_bar(others_color = "#BDBDBD", legend_text_italic = FALSE) +
  scale_fill_manual(values = taxa_colors) +
  guides(fill = guide_legend(reverse = TRUE)) +
  labs(x = NULL, y = "Relative abundance (%)", fill = "Genus") +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
        axis.text.y = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        strip.text = element_text(size = 14, face = "bold"),
        legend.position = "right",
        legend.title = element_text(size = 13),
        legend.text = element_text(size = 12),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 1))

p_genus



## -----------------------------------------------------------------------------------------###
## ------------------------------Line--------------------------------------------------###

sorghum$sample_table$SitePlanting <- paste(
  sorghum$sample_table$Site,
  sorghum$sample_table$Planting
)

# Desired order
sorghum$sample_table$SitePlanting <- factor(
  sorghum$sample_table$SitePlanting,
  levels = c(
    "Kansas Early",
    "Kansas Regular",
    "South Dakota Early",
    "South Dakota Regular"
  )
)

# Make sure Line is character
sorghum$sample_table$Line <- as.character(sorghum$sample_table$Line)

# Check
table(sorghum$sample_table$Line, sorghum$sample_table$SitePlanting)

abund_line <- trans_abund$new(dataset = sorghum, taxrank = "Phylum", ntaxa = 10)
abund_line

phylum_dat <- abund_line$data_abund
phylum_dat

phylum_dat$Line <- as.character(phylum_dat$Line)

phylum_dat$SitePlanting <- paste( phylum_dat$Site, phylum_dat$Planting)

phylum_dat <- phylum_dat[, c(
    "Sample",
    "Line",
    "SitePlanting",
    "Taxonomy",
    "Abundance"
  )]

phylum_dat

unique(phylum_dat$SitePlanting)

unique(phylum_dat$Taxonomy)

# ============================================================
# 4. CALCULATE "OTHERS" FOR EACH SAMPLE
# ============================================================

sample_sum <- aggregate( Abundance ~ Sample + Line + SitePlanting, data = phylum_dat, FUN = sum)

sample_sum$Taxonomy <- "Others"

sample_sum$Abundance <- 100 - sample_sum$Abundance

# Prevent very small negative values caused by rounding
sample_sum$Abundance[sample_sum$Abundance < 0 ] <- 0

others_dat <- sample_sum[ ,
  c(
    "Sample",
    "Line",
    "SitePlanting",
    "Taxonomy",
    "Abundance"
  )]

phylum_all <- rbind( phylum_dat, others_dat)


# ============================================================
# 5. MEAN ABUNDANCE FOR EACH LINE × SITE × PLANTING
# ============================================================

phylum_mean <- aggregate( Abundance ~ Line + SitePlanting + Taxonomy, data = phylum_all, FUN = mean)

head(phylum_mean)

dim(phylum_mean)

phylum_mean$SitePlanting <- factor( phylum_mean$SitePlanting, levels = c(
    "Kansas Early",
    "Kansas Regular",
    "South Dakota Early",
    "South Dakota Regular"
  ))

table( phylum_mean$Line, phylum_mean$SitePlanting)

plot_colors <- c( its_colorset, "Others" = "grey70")

# ============================================================
# 6. MAKE ONE PLOT FOR EACH SORGHUM LINE
# ============================================================

line_names <- sort(unique(phylum_mean$Line))

phylum_line_plots <- list()


for (i in line_names) { df_i <- phylum_mean[ phylum_mean$Line == i, ]
  
  df_i$SitePlanting <- factor(
    df_i$SitePlanting,
    levels = c(
      "Kansas Early",
      "Kansas Regular",
      "South Dakota Early",
      "South Dakota Regular"
    )
  )
  
  
  p <- ggplot( df_i, aes( x = SitePlanting, y = Abundance, fill = Taxonomy )) +
    geom_col( width = 0.85) +
    scale_fill_manual( values = plot_colors) +
    scale_y_continuous( limits = c(0, 100), expand = c(0, 0)) +
    labs( title = i,
      x = NULL,
      y = "Relative abundance (%)",
      fill = "Phylum") +
    theme_minimal() + theme(axis.text.x = element_text( angle = 25, hjust = 1, size = 14, colour = "black" ), 
                            axis.text.y = element_text( size = 14, colour = "black"),
                            axis.title.y = element_text(size = 14, colour = "black"),
                            plot.title = element_text(size = 15, face = "bold", hjust = 0.5), 
                            legend.position = "right",
                            legend.title = element_text(size = 13),
                            legend.text = element_text(size = 12),
                            panel.grid.major.x = element_blank(),
                            panel.grid.minor = element_blank(),
                            panel.border = element_rect(colour = "black", fill = NA, linewidth = 1))
  phylum_line_plots[[i]] <- p
}

names(phylum_line_plots)

phylum_line_plots[["PI576376"]]

# Create abundance folder
dir.create("abundance", showWarnings = FALSE)

# Save every line plot as PDF
for (i in names(phylum_line_plots)) {
  ggsave(filename = file.path( "abundance", paste0(i, "_Phylum_Abundance.pdf")),
    plot = phylum_line_plots[[i]], width = 8, height = 6, units = "in" ) }
