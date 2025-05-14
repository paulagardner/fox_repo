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
excluded_samples <- c(excluded_samples, "YPI1082")
excluded_samples

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





######################################ADMIXTURE###########################################################
library(tidyverse)
library(RColorBrewer)

# Paths
q_file_base <- "/gpfs/data/bergstrom/paula/fox_repo/variant_calling/admixture/"
output_dir <- "/gpfs/data/bergstrom/paula/fox_repo/plotting/"

# Read sample metadata and filter
df <- read.table(file.path(output_dir, "order.tsv"), sep = "\t", header = FALSE, stringsAsFactors = FALSE)
colnames(df) <- c("sample", "region", "continent")
df

###we already do this up above
excluded_samples <- readLines(file.path(q_file_base, "lowcoverage_missingness_removed_samples.txt"))
df_filtered <- df[!(df$sample %in% excluded_samples), ]
df_filtered

# Read final sample list that matches Q file row order
final_samples <- readLines(file.path(q_file_base, "lowcoverage_missingness_final_samples.txt"))
final_samples <- final_samples[final_samples %in% df_filtered$sample]
final_samples

# Get list of all .Q files and extract K values
q_files <- list.files(q_file_base, pattern = "lowcoverage_missingness\\.\\d+\\.Q$", full.names = TRUE)
k_values <- as.integer(str_extract(basename(q_files), "(?<=\\.)\\d+(?=\\.Q)"))

# Color palette (adapt dynamically if needed)
max_k <- max(k_values)
color_pal <- RColorBrewer::brewer.pal(min(12, max_k), "Paired")
continent_colors <- setNames(RColorBrewer::brewer.pal(n = length(unique(df_filtered$continent)), name = "Set2"),
  unique(df_filtered$continent))


admix_plot_list = list()
# Loop through each K and create plot
for (i in seq_along(q_files)) {
  q_file_path <- q_files[i]
  k <- k_values[i]
  
  # Read Q file
  q_file <- read.table(q_file_path, header = FALSE)
  q_file$sample <- final_samples
  
  # Merge with ordered metadata
  q_merged <- left_join(df_filtered, q_file, by = "sample")
  stopifnot(nrow(q_merged) == nrow(df_filtered))
  
  # Melt and prepare long-format data
  kdf2 <- q_merged %>%
    mutate(row_order = row_number()) %>%
    pivot_longer(
      cols = starts_with("V"),
      names_to = "popGroup",
      values_to = "prop"
    ) %>%
    mutate(
      popGroup = as.integer(str_remove(popGroup, "V")),
      sample = factor(sample, levels = df_filtered$sample),
      popGroup = as.factor(popGroup)
    )
  kdf2
  
  # Label data
  label_df <- df_filtered %>%
    mutate(
      sample = factor(sample, levels = df_filtered$sample),
      label = paste0(sample, " (", region, ")")
    )
  label_df
  # Colors for this K
  fill_colors <- setNames(color_pal[seq_len(k)], as.character(seq_len(k)))
  
  # Plot
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
      title = paste("ADMIXTURE Plot K =", k),
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
  
  # Save
  output_file <- file.path(output_dir, paste0("admixplot_YPI1082removed", k, ".png"))
  ggsave(output_file, plot = admix_final, width = 12, height = 10)
  message("Saved ADMIXTURE plot for K = ", k)
  
  admix_plot_list[[paste0("K", k)]] <- admix_final

}


library(patchwork)

combined_plot <- wrap_plots(admix_plot_list, ncol = 1)  # vertical stack

ggsave("admixplot_grid_K3-7_patchwork.png", plot = combined_plot, width = 12, height = 10 * length(admix_plot_list), limitsize=FALSE)
