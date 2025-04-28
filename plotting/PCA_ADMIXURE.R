# --------- LOAD LIBRARIES ----------
library(tidyverse) # for hpc: followed https://blog.zenggyu.com/posts/en/2018-01-29-installing-r-r-packages-e-g-tidyverse-and-rstudio-on-ubuntu-linux/index.html
#module loaded curl/7.76.1/gcc in bash, then re-entered R
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
continent_order <- c( "Europe", "Africa", "WestAsia", "Asia","Americas") 

# Reorder the 'continent' factor based on the custom order
metadata$continent <- factor(metadata$continent, levels = continent_order)

# Join metadata to PCA data
pca_meta <- left_join(pca, metadata, by = "sample")

# Define continent color palette using ggplot2-style Set2
continents <- unique(metadata$continent)
continent_colors <- setNames(RColorBrewer::brewer.pal(n = length(continents), name = "Set2"), continents)

# ---------------------- PCA PLOT ----------------------

# Plot PCA with consistent continent colors
pca_plot <- ggplot(pca_meta, aes(x = PC1, y = PC2, color = continent)) +
  geom_point(size = 3) +
  geom_text_repel(aes(label = sample), size = 3, box.padding = 0.35, max.overlaps = 23) +
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


# Set the value of K for testing
# --------- LOOP OVER K VALUES ----------

# Set the value of K for testing
k_val <- 3

# Build file paths for K = 3
q_file <- paste0("/gpfs/data/bergstrom/paula/fox_repo/variant_calling/admixture/low_missingness.", k_val, ".Q")
output_file <- paste0("/gpfs/data/bergstrom/paula/fox_repo/plotting/ADMIXTURE_plot_K", k_val, ".png")

# Load Q-matrix and assign sample IDs
qmat <- read.table(q_file)
colnames(qmat) <- paste0(1:k_val)
samples <- read.table("/gpfs/data/bergstrom/paula/fox_repo/variant_calling/admixture/lowcoverage_missingness_final_samples.txt", header = FALSE)$V1
qmat <- qmat %>% mutate(ID = samples)

# Join with metadata
qmat_meta <- left_join(qmat, metadata, by = c("ID" = "sample"))

# Enforce continent factor order
qmat_meta$continent <- factor(qmat_meta$continent, levels = continent_order)

# Sort individuals: first by continent, then by dominant ancestry
qmat_meta <- qmat_meta %>%
  mutate(Dominant_K = colnames(select(., matches("^\\d+$")))[max.col(select(., matches("^\\d+$")), ties.method = "first")]) %>%
  mutate(Dominant_K = factor(Dominant_K, levels = as.character(1:k_val))) %>%
  arrange(continent, Dominant_K)

# Set ID factor levels to lock in this order
qmat_meta$ID <- factor(qmat_meta$ID, levels = qmat_meta$ID)

# Reshape to long format
ancestry_cols <- as.character(1:k_val)

q_long <- qmat_meta %>%
  pivot_longer(cols = all_of(ancestry_cols), names_to = "Cluster", values_to = "Proportion")

# Assign ggplot2 palette to clusters
k_palette <- brewer.pal(max(k_val, 3), "Set1")[1:k_val]

# Build label info
id_levels <- levels(qmat_meta$ID)
continent_lookup <- metadata %>% select(sample, continent)
continent_per_id <- left_join(data.frame(ID = id_levels), continent_lookup, by = c("ID" = "sample"))

label_df <- data.frame(
  x = seq_along(id_levels),
  y = -0.05,   # Adjusted position lower to avoid overlap
  label = paste0(id_levels, " (", continent_per_id$continent, ")"),
  continent = continent_per_id$continent
)
label_df$color <- continent_colors[label_df$continent]


  # Base ADMIXTURE plot
# Base ADMIXTURE plot with only outlines around bars
admix_plot_base <- ggplot(q_long, aes(x = ID, y = Proportion, fill = Cluster)) +
  geom_bar(stat = "identity", width = 1, color = "black") +  # Add black outline to bars
  scale_fill_manual(values = k_palette) +
  theme_minimal() +
  labs(title = paste("ADMIXTURE Plot K =", k_val),
       x = "Individuals (grouped by continent)",
       y = "Ancestry Proportion",
       fill = "Cluster") +
  theme(
    axis.text.x = element_blank(),        # Remove x-axis labels
    axis.ticks.x = element_blank(),       # Remove x-axis ticks
    axis.line.x = element_blank(),        # Remove x-axis line
    panel.grid.major = element_blank(),   # Remove major grid lines
    panel.grid.minor = element_blank(),   # Remove minor grid lines
    panel.border = element_rect(color = "black", fill = NA, size = 1),  # Add border around the panel
    plot.margin = margin(20, 20, 40, 20)  # Adjust plot margins to prevent clipping
  )

  # Base ADMIXTURE plot with only outlines around bars and adjusted label spacing
admix_plot_base <- ggplot(q_long, aes(x = ID, y = Proportion, fill = Cluster)) +
  geom_bar(stat = "identity", width = 1, color = "black") +  # Add black outline to bars
  scale_fill_manual(values = k_palette) +
  theme_minimal() +
  labs(title = paste("ADMIXTURE Plot K =", k_val),
       x = "Individuals (grouped by continent)",
       y = "Ancestry Proportion",
       fill = "Cluster") +
  theme(
    axis.text.x = element_blank(),        # Remove x-axis labels
    axis.ticks.x = element_blank(),       # Remove x-axis ticks
    axis.line.x = element_blank(),        # Remove x-axis line
    panel.grid.major = element_blank(),   # Remove major grid lines
    panel.grid.minor = element_blank(),   # Remove minor grid lines
    #panel.border = element_rect(color = "black", fill = NA, size = 1),  # Add border around the panel
    plot.margin = margin(30, 20, 20, 20),  # Increase bottom margin to prevent overlap
    plot.title = element_text(size = 20, hjust = 0.5, margin = margin(b = 80)),  # Adjust title size and bottom margin
    axis.title.x = element_text(size = 14, margin = margin(t = 80)),  # Increase space between axis title and bars
    axis.title.y = element_text(size = 14, margin = margin(r = 80))  # Increase space between y-axis title and bars
  )


label_df <- label_df %>%
  mutate(y_adjusted = y + 0.06)  # Move closerto the admixture bar plots(positive = upward, negative = downward)

# Final plot with labels
# Final plot with labels
admix_final <- admix_plot_base +
  geom_text(
    data = label_df,
    aes(x = x,  y = y_adjusted, label = label, color = continent),
    angle = 45,
    hjust = 1,
    size = 2.5,
    inherit.aes = FALSE
  ) +
  scale_color_manual(values = continent_colors) +
  coord_cartesian(clip = "off") +
  theme(
    legend.position = "right"
  )


# Save the plot
print(admix_final)  # Ensure the plot is rendered in your session
ggsave(output_file, plot = admix_final, width = 12, height = 10)

message("Saved: ", output_file)




