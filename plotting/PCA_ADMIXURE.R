# --------- LOAD LIBRARIES ----------
library(tidyverse) # for hpc: followed https://blog.zenggyu.com/posts/en/2018-01-29-installing-r-r-packages-e.g-tidyverse-and-rstudio-on-ubuntu-linux/index.html
library(googlesheets4)
library(ggrepel)
library(RColorBrewer)
library(grid)
library(gtable)

# ---------------------- IMPORT & PREP ----------------------

# Load PCA eigenvector file
pca <- read_tsv("/gpfs/bio/xrq24scu/fox_repo/variant_calling/Eigenstrat/Vvulpestidy.evec")


# Rename columns to PC1, PC2, etc.
colnames(pca)[2:ncol(pca)] <- paste0("PC", 1:(ncol(pca)-1))

# Import metadata from Google Sheets
sheet_url <- "https://docs.google.com/spreadsheets/d/1KBqsCfPiuyno90iRJ-O0wV0HIijFncstvcIpGnpGY0k/edit?gid=1838883534#gid=1838883534"
gs4_auth(path = "/gpfs/bio/xrq24scu/fox_repo/plotting/googlesheetskey.json")
metadata <- read_sheet(sheet_url) %>%
  mutate(sample = as.character(sample))

# Manually specify the continent order (adjust this order as needed)
continent_order <- c("Europe", "Africa", "WestAsia", "Asia", "Americas")

# Reorder the 'continent' factor based on the custom order
metadata$continent <- factor(metadata$continent, levels = continent_order)

# Ensure 'region' is included in the metadata (you may need to adjust this if your sheet contains a 'region' column)
# Assuming 'region' is part of the metadata after importing from Google Sheets
# If not, you'll need to adjust how you load the region information.

# Join metadata to PCA data
pca_meta <- left_join(pca, metadata, by = "sample")
# Remove specific samples from the PCA plot based on .txt file from ADMIXTURE filtering script. CHANGE THIS
excluded_samples <- c("S080738", "S100038", "S142040")
pca_meta_filtered <- pca_meta %>%
  filter(!sample %in% excluded_samples)

# Define continent color palette using ggplot2-style Set2
continents <- unique(metadata$continent)
continent_colors <- setNames(RColorBrewer::brewer.pal(n = length(continents), name = "Set2"), continents)

# ---------------------- PCA PLOT ----------------------

  # Plot PCA with consistent continent colors
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

# Define the directory and base filenames
q_file_base <- "/gpfs/data/bergstrom/paula/fox_repo/variant_calling/admixture/low_missingness"
output_dir <- "/gpfs/data/bergstrom/paula/fox_repo/plotting/"

# Load sample IDs
samples <- read.table("/gpfs/data/bergstrom/paula/fox_repo/variant_calling/admixture/lowcoverage_missingness_final_samples.txt", header = FALSE)$V1

# Find all available Q files (assumes filenames like "low_missingness.K.Q")
q_files <- list.files(path = dirname(q_file_base), pattern = "low_missingness\\.\\d+\\.Q", full.names = TRUE)

# Extract available K values
k_vals <- as.integer(gsub(".*\\.(\\d+)\\.Q$", "\\1", q_files))

# Sort K values
k_vals <- sort(k_vals)

# Prepare metadata lookup outside of loop (since it doesn't change per K)
continent_lookup <- metadata %>% select(sample, continent, region)  # Include 'region' here

# ---------- LOOP OVER K VALUES ----------

for (k_val in k_vals) {

  message("Processing K = ", k_val)

  # --- File Paths ---
  q_file <- paste0(q_file_base, ".", k_val, ".Q")
  output_file <- paste0(output_dir, "ADMIXTURE_plot_K", k_val, ".png")

  # --- Load Data ---
  qmat <- read.table(q_file, col.names = as.character(1:k_val)) %>%
    mutate(ID = samples)

  # Join with metadata and arrange
  qmat_meta <- left_join(qmat, metadata, by = c("ID" = "sample")) %>%
    mutate(
      continent = factor(continent, levels = continent_order),
      Dominant_K = factor(
        colnames(select(., matches("^\\d+$")))[max.col(select(., matches("^\\d+$")), ties.method = "first")],
        levels = as.character(1:k_val)
      )
    ) %>%
    arrange(continent, Dominant_K) %>%
    mutate(ID = factor(ID, levels = ID))  # Lock individual order

  # --- Prepare Data for Plotting ---
  ancestry_cols <- paste0("X", 1:k_val)

  q_long <- qmat_meta %>%
    pivot_longer(cols = all_of(ancestry_cols), names_to = "Cluster", values_to = "Proportion")

  id_levels <- levels(qmat_meta$ID)

  # Build label data
  label_df <- left_join(data.frame(ID = id_levels), continent_lookup, by = c("ID" = "sample")) %>%
    mutate(
      x = seq_along(ID),
      y = -0.01,  # Control label vertical position
      label = paste0(ID, " (", region, ")"),  # Use region here
      color = continent_colors[continent]  # Keep the continent color for consistency
    )

  # --- Build Plot Base ---
  admix_plot_base <- ggplot(q_long, aes(x = ID, y = Proportion, fill = Cluster)) +
    geom_bar(stat = "identity", width = 1, color = "black") +
    scale_fill_manual(values = brewer.pal(max(k_val, 3), "Set1")[1:k_val]) +
    theme_minimal() +
    labs(
      title = paste("ADMIXTURE Plot K =", k_val),
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

  # --- Add Labels ---
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

  # --- Save Plot ---
  ggsave(output_file, plot = admix_final, width = 12, height = 10)

  message("Saved plot for K = ", k_val)

}  # <- CORRECTLY closes the `for` loop
