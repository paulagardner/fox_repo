# --------- LOAD LIBRARIES ----------
library(tidyverse)
library(googlesheets4)
library(ggrepel)
library(RColorBrewer)
library(grid)
library(gtable)
library(dplyr)

# ---------------------- IMPORT & PREP ----------------------

# Load PCA eigenvector file
#pca <- read_tsv("/gpfs/data/bergstrom/paula/fox_repo/variant_calling/Eigenstrat/Vvulpestidy.evec")
pca <- read_tsv("/gpfs/data/bergstrom/paula/fox_repo/variant_calling/Eigenstrat/Vvulpes_rmoutliers_tidy.evec")

# Rename columns to PC1, PC2, etc.
colnames(pca)[2:ncol(pca)] <- paste0("PC", 1:(ncol(pca) - 1))

# Import metadata from Google Sheets
sheet_url <- "https://docs.google.com/spreadsheets/d/1KBqsCfPiuyno90iRJ-O0wV0HIijFncstvcIpGnpGY0k/edit?gid=1838883534#gid=1838883534"
gs4_auth(path = "/gpfs/bio/xrq24scu/fox_repo/plotting/googlesheetskey.json")
metadata <- read_sheet(sheet_url) %>%
  mutate(sample = as.character(sample))

# Specify continent order manually (adjust as needed)
continent_order <- c("Europe", "Africa", "WestAsia", "Asia", "Americas")
metadata$continent <- factor(metadata$continent, levels = continent_order)

# Join metadata to PCA data
pca_meta <- left_join(pca, metadata, by = "sample")

# Remove specific samples from PCA plot (based on ADMIXTURE filtering script)
excluded_samples <- readLines("/gpfs/data/bergstrom/paula/fox_repo/variant_calling/admixture/lowcoverage_missingness_removed_samples.txt")

pca_meta_filtered <- pca_meta %>%
  filter(!sample %in% excluded_samples)
pca_meta_filtered
# Define color palette for continents
continents <- unique(metadata$continent)
continent_colors <- setNames(RColorBrewer::brewer.pal(n = length(continents), name = "Set2"), continents)

# ---------------------- PCA PLOT ----------------------

pca_plot <- ggplot(pca_meta_filtered, aes(x = PC1, y = PC2, color = continent)) +
  geom_point(size = 3) +
  geom_text_repel(aes(label = paste0(sample, " (", region, ")")), size = 3, box.padding = 0.5, max.overlaps = 35) +
  theme_minimal(base_size = 15) +
  labs(title = "PCA Plot by Continent", x = "Principal Component 1", y = "Principal Component 2") +
  theme(
    plot.title = element_text(size = 20, hjust = 0.5),
    axis.title = element_text(size = 15),
    axis.text = element_text(size = 12),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14)
  ) +
  scale_color_manual(values = continent_colors)

# Show and save PCA plot
print(pca_plot)
ggsave("/gpfs/bio/xrq24scu/fox_repo/plotting/PCA_plot_large.png", plot = pca_plot, width = 12, height = 10)

# ---------------------- ADMIXTURE PLOTS ----------------------

q_file_base <- "/gpfs/data/bergstrom/paula/fox_repo/variant_calling/admixture/lowcoverage_missingness/"
output_dir <- "/gpfs/data/bergstrom/paula/fox_repo/plotting/"

#made this by manually manipulating google sheets info
df <- read.table("/gpfs/data/bergstrom/paula/fox_repo/plotting/order.tsv", sep = "\t", header = FALSE, stringsAsFactors = FALSE)
colnames(df) <- c("sample", "region", "continent")


# Read the sample IDs to exclude
excluded_samples <- readLines("/gpfs/data/bergstrom/paula/fox_repo/variant_calling/admixture/lowcoverage_missingness_removed_samples.txt")
excluded_samples

# Filter out the excluded samples, preserving order
df_filtered <- df[!(df$sample %in% excluded_samples), ]
df_filtered
seq_len(nrow(df_filtered))
df_filtered


q_file <- read.table("/gpfs/data/bergstrom/paula/fox_repo/variant_calling/admixture/lowcoverage_missingness.3.Q")

q_file <- q_file %>%
  mutate(sample = df_filtered$sample,
         region = df_filtered$region,
         continent = df_filtered$continent,
         row_order=row_number())

q_file

# Pivot without renaming the columns
kdf2 <- q_file %>%
  pivot_longer(
    cols = -c(sample, region, continent, row_order),
    names_to = "popGroup",
    values_to = "prop"
  ) %>%
  mutate(
    popGroup = as.integer(factor(popGroup, levels = unique(popGroup)))  # e.g. V1 → 1, etc.
  ) %>%
  arrange(row_order, popGroup) %>%  # preserve sample order
  select(sample, continent, prop, region, popGroup)  # final format
kdf2

kdf2 <- kdf2 %>%
  mutate(popGroup = as.factor(popGroup))


## Optional: arrange by sample and population group
####kdf2 <- kdf2 %>% arrange(sample, popGroup)

#fill_colors <- setNames(brewer.pal(max(3, 3), "Paired"), paste0("X", 1:3))
fill_colors <- setNames(brewer.pal(3, "Paired"), as.character(1:3))

samples <- df_filtered$sample  # or a subset of it
label_df <- df_filtered %>%
    filter(sample %in% samples) %>%
    mutate(
      sample = factor(sample, levels = samples),
      label = paste0(sample, " (", region, ")")
    )
label_df


admix_final <- ggplot(kdf2, aes(x = sample, y = prop, fill = popGroup)) +
    geom_col(width = 1, color = "black") +
    geom_text(
      data = label_df,
      aes(x = sample, y = -0.01, label = label, color = continent),
      inherit.aes = FALSE,
      angle = 45,
      hjust = 1,
      size = 2.5
    ) +
    scale_fill_manual(values = fill_colors) +
    scale_color_manual(values = continent_colors) +
    coord_cartesian(clip = "off", ylim = c(-0.05, 1)) +
    theme_minimal() +
    labs(
      title = paste("ADMIXTURE Plot K =", 3),
      x = "Samples",
      y = "Ancestry Proportion",
      fill = "Cluster"
    ) +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.line.x = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      plot.margin = margin(30, 20, 20, 20),
      plot.title = element_text(size = 20, hjust = 0.5, margin = margin(b = 80)),
      axis.title.x = element_text(size = 14, margin = margin(t = 80)),
      axis.title.y = element_text(size = 14, margin = margin(r = 80)),
      legend.position = "right"
    )

  # Save plot
  output_file <- paste0(output_dir, "ADMIXTURE_plot_K", 3, ".png")
  ggsave(output_file, plot = admix_final, width = 12, height = 10)
  message("Saved ADMIXTURE plot for K = ", k_val)
}
