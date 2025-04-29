# --------- LOAD LIBRARIES ----------
library(tidyverse)
library(googlesheets4)
library(ggrepel)
library(RColorBrewer)
library(grid)
library(gtable)

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

######### ADMIXTURE PLOTS #########

# Define directory and file paths
q_file_base <- "/gpfs/data/bergstrom/paula/fox_repo/variant_calling/admixture/low_missingness"
output_dir <- "/gpfs/data/bergstrom/paula/fox_repo/plotting/"

# Load sample IDs
samples <- read.table("/gpfs/data/bergstrom/paula/fox_repo/variant_calling/admixture/lowcoverage_missingness_final_samples.txt", header = FALSE)$V1

# Find all available Q files
q_files <- list.files(path = dirname(q_file_base), pattern = "low_missingness\\.\\d+\\.Q", full.names = TRUE)

# Extract and sort K values
k_vals <- sort(as.integer(gsub(".*\\.(\\d+)\\.Q$", "\\1", q_files)))

# Metadata lookup (doesn't change per K)
continent_lookup <- metadata %>% select(sample, continent, region)

# ---------- LOOP OVER K VALUES ----------
# Loop over K values
# --- Loop Over K Values ---
for (k_val in k_vals) {

  message("Processing K = ", k_val)

  # --- File Paths ---
  q_file <- paste0(q_file_base, ".", k_val, ".Q")
  output_file <- paste0(output_dir, "ADMIXTURE_plot_K", k_val, ".png")
  output_file_alt <- paste0(output_dir, "ADMIXTURE_plot_K", k_val, "_byAncestry.png")

  # --- Load Data ---
  qmat <- read.table(q_file, col.names = paste0("X", 1:k_val)) %>%
    mutate(ID = samples)

  # --- Original Plot: Reorder by Continent ---
  qmat_meta_original <- left_join(qmat, metadata, by = c("ID" = "sample")) %>%
    mutate(
      continent = factor(continent, levels = continent_order)  # Reorder by continent only
    ) %>%
    arrange(continent) %>%  # Only arrange by continent
    mutate(ID = factor(ID, levels = ID))  # Lock individual order

  qmat_meta_original <- left_join(qmat, metadata, by = c("ID" = "sample")) %>%
  mutate(
    continent = factor(continent, levels = continent_order)  # Reorder by continent only
  ) %>%
  arrange(continent, ID) %>%  # First arrange by continent, then by ID
  mutate(ID = factor(ID, levels = ID))  # Lock individual order


  # --- Ancestry Plot: Reorder by Dominant Ancestry ---
  ancestry_cols <- paste0("X", 1:k_val)

  qmat_meta_ancestry <- left_join(qmat, metadata, by = c("ID" = "sample")) %>%
    mutate(
      continent = factor(continent, levels = continent_order)  # Reorder by continent
    ) %>%
    # Calculate Dominant Ancestry Dynamically for Each K
    mutate(
      Dominant_K = as.character(apply(select(., all_of(ancestry_cols)), 1, function(x) {
        which.max(x)  # Identify the index of the maximum ancestry proportion
      })),
      dominant_prop = apply(select(., all_of(ancestry_cols)), 1, max)  # Store the dominant proportion
    ) %>%
    arrange(Dominant_K, desc(dominant_prop)) %>%  # First by dominant ancestry, then by dominant proportion
    mutate(ID = factor(ID, levels = ID))  # Lock this new order

  ##################
  # Original Plot (Ordered by Continent)
  ##################
  
  q_long_original <- qmat_meta_original %>%
    pivot_longer(cols = all_of(ancestry_cols), names_to = "Cluster", values_to = "Proportion") %>%
    mutate(ID = factor(ID, levels = qmat_meta_original$ID[order(qmat_meta_original$continent)]))  # Reorder by continent

  label_df <- left_join(data.frame(ID = samples), continent_lookup, by = c("ID" = "sample")) %>%
    mutate(
      x = seq_along(ID),
      y = -0.01,  # Control label vertical position
      label = paste0(ID, " (", region, ")"),
      color = continent_colors[continent]  # Keep the continent color for consistency
    )


  admix_plot_base <- ggplot(q_long_original, aes(x = ID, y = Proportion, fill = Cluster)) +
    geom_bar(stat = "identity", width = 1, color = "black") +
    scale_fill_manual(values = brewer.pal(max(k_val, 3), "Set1")[1:k_val]) +
    theme_minimal() +
    labs(
      title = paste("ADMIXTURE Plot K =", k_val, "- Original Sample Order"),
      x = "Individuals (grouped by continent)",
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
      axis.title.y = element_text(size = 14, margin = margin(r = 80))
    )

  admix_final <- admix_plot_base +
    geom_text(
      data = label_df,
      aes(x = x, y = y, label = label, color = continent),
      angle = 45,
      hjust = 1,
      size = 2.5,
      inherit.aes = FALSE
    ) +
    scale_color_manual(values = continent_colors) +
    coord_cartesian(clip = "off", ylim = c(-0.05, 1)) +
    theme(legend.position = "right")

  ggsave(output_file, plot = admix_final, width = 12, height = 10)
  message("Saved original plot for K = ", k_val)

  ##################
  # By Ancestry Plot (Ordered by Dominant Ancestry)
  ##################

  q_long_ancestry <- qmat_meta_ancestry %>%
    pivot_longer(cols = all_of(ancestry_cols), names_to = "Cluster", values_to = "Proportion")  # Already ordered by Dominant_K

  label_df_alt <- left_join(data.frame(ID = levels(qmat_meta_ancestry$ID)), continent_lookup, by = c("ID" = "sample")) %>%
    mutate(
      x = seq_along(ID),
      y = -0.01,
      label = paste0(ID, " (", region, ")"),
      color = continent_colors[continent]
    )

  admix_plot_base_alt <- ggplot(q_long_ancestry, aes(x = ID, y = Proportion, fill = Cluster)) +
    geom_bar(stat = "identity", width = 1, color = "black") +
    scale_fill_manual(values = brewer.pal(max(k_val, 3), "Set1")[1:k_val]) +
    theme_minimal() +
    labs(
      title = paste("ADMIXTURE Plot K =", k_val, "- Ordered by Dominant Ancestry and Proportion"),
      x = "Individuals (grouped by ancestry)",
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
      axis.title.y = element_text(size = 14, margin = margin(r = 80))
    )

  admix_final_alt <- admix_plot_base_alt +
    geom_text(
      data = label_df_alt,
      aes(x = x, y = y, label = label, color = continent),
      angle = 45,
      hjust = 1,
      size = 2.5,
      inherit.aes = FALSE
    ) +
    scale_color_manual(values = continent_colors) +
    coord_cartesian(clip = "off", ylim = c(-0.05, 1)) +
    theme(legend.position = "right")

  ggsave(output_file_alt, plot = admix_final_alt, width = 12, height = 10)
  message("Saved ancestry-ordered plot for K = ", k_val)

}  # <- Close loop
