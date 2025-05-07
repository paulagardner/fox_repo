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
pca <- read_tsv("/gpfs/bio/xrq24scu/fox_repo/variant_calling/Eigenstrat/Vvulpestidy.evec")

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
excluded_samples <- c("S080738", "S100038", "S142040")
pca_meta_filtered <- pca_meta %>%
  filter(!sample %in% excluded_samples)

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



#------------------------------------------ADMIXTURE PLOTS----------------------------------------------

########## ADMIXTURE PLOTS #########

# Define file paths
q_file_base <- "/gpfs/data/bergstrom/paula/fox_repo/variant_calling/admixture/lowcoverage_missingness"
output_dir <- "/gpfs/data/bergstrom/paula/fox_repo/plotting/"

# Read sample order
order_df <- read.table("/gpfs/data/bergstrom/paula/fox_repo/plotting/order.tsv", header = FALSE, sep = "\t")
colnames(order_df) <- c("sample", "region", "continent")

# Read final list of ADMIXTURE samples
final_sample_list <- read.table("/gpfs/data/bergstrom/paula/fox_repo/variant_calling/admixture/lowcoverage_missingness_final_samples.txt", header = FALSE)$V1
samples <- order_df$sample[order_df$sample %in% final_sample_list]

# Filter and arrange metadata
metadata_filtered <- metadata %>%
  filter(sample %in% samples) %>%
  arrange(match(sample, samples))

# Read Q matrix files
q_files <- list.files(path = dirname(q_file_base), pattern = "lowcoverage_missingness\\.\\d+\\.Q", full.names = TRUE)
k_vals <- sort(as.integer(gsub(".*\\.(\\d+)\\.Q$", "\\1", q_files)))

for (k_val in k_vals) {
  message("Processing K = ", k_val)

  # File paths
  q_file <- paste0(q_file_base, ".", k_val, ".Q")
  output_file <- paste0(output_dir, "ADMIXTURE_plot_K", k_val, ".png")

  # Read Q matrix
  qmat <- read.table(q_file, header = FALSE)
  colnames(qmat) <- paste0("X", 1:k_val)

  # Combine with sample names
  qmat_ordered <- bind_cols(tibble(sample = samples), qmat)

  # Reshape and enforce factor levels
  q_long <- qmat_ordered %>%
    pivot_longer(
      cols = starts_with("X"),
      names_to = "Cluster",
      values_to = "Proportion"
    ) %>%
    mutate(
      sample = factor(sample, levels = samples),
      Cluster = factor(Cluster, levels = paste0("X", 1:k_val))  # Important: fix stacking order
    )

  # Set fill colors for clusters
  fill_colors <- setNames(brewer.pal(max(k_val, 3), "Paired"), paste0("X", 1:k_val))

  # Create label metadata
  label_df <- order_df %>%
    filter(sample %in% samples) %>%
    mutate(
      sample = factor(sample, levels = samples),
      label = paste0(sample, " (", region, ")")
    )

  # Plot
  admix_final <- ggplot(q_long, aes(x = sample, y = Proportion, fill = Cluster)) +
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
      title = paste("ADMIXTURE Plot K =", k_val),
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
  ggsave(output_file, plot = admix_final, width = 12, height = 10)
  message("Saved ADMIXTURE plot for K = ", k_val)
}

# Optional: diagnostic sample order check
order_check_plot <- ggplot(q_long, aes(x = sample, y = Proportion, fill=Cluster)) +
  geom_col() +
  theme(axis.text.x = element_text(angle = 90))

ggsave("/gpfs/data/bergstrom/paula/fox_repo/plotting/sample_order_check.png", plot = order_check_plot, width = 12, height = 4)
