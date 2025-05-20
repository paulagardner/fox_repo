# --------- LOAD LIBRARIES ----------
library(tidyverse)
library(googlesheets4)
library(ggrepel)
library(RColorBrewer)
library(grid)
library(gtable)
library(dplyr)
library(tidyverse)
library(RColorBrewer)
library(patchwork)

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
excluded_samples <- readLines("/gpfs/data/bergstrom/paula/fox_repo/variant_calling/admixture/low_missingness_rm_outlier_removed_samples.txt")
#excluded_samples <- c(excluded_samples, "YPI1082")
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







######################################CLUMPAK##########################################################

# Paths
q_file_base <- "/gpfs/data/bergstrom/paula/fox_repo/variant_calling/admixture/"
output_dir <- "/gpfs/data/bergstrom/paula/fox_repo/plotting/"

# Read final sample list that matches Q file row order
#better to use the .fam file, as the final list generated from the .tmp file is NOT comprehensive and caused ordering issues!!!
final_samples <- read.table(file.path(q_file_base, "low_missingness_rm_outlier_snpmanualfilter.fam"), sep="\t", header=FALSE)
final_samples <- final_samples[2]
final_samples

final_samples <- final_samples %>%
  rename(sample = V2) 
final_samples
#the following works but changes the ordering of the samples to that of the metadata file. We don't want that
#as the samples must be in the order that the q matrix is in

# Assume:
# sample_list is your character vector of samples (in desired order)
# full_df is your large data frame with a column "sample"

# 1. Convert the sample list to a data frame
samples_df <- data.frame(sample = final_samples, stringsAsFactors = FALSE)
samples_df
# 2. Join with the full_df on "sample", keeping the order of samples_df
# This will add all columns from full_df to samples_df, matching on "sample"

# Extract continent values in the order of the .Q file
qfileorder_df <- metadata[match(final_samples$sample, metadata$sample), ]

qfileorder_df


#now get population values as numbers for clumpak:
clumpak_df <- qfileorder_df[c("sample","continent")] 
clumpak_df
popID_df <- qfileorder_df[c("continent")]
popID_df

#make this df read the continent values as factors, since clumpak wants input 
#for population in number format
popID_df$continent <- as.numeric((clumpak_df$continent))

write.table(popID_df, file = file.path(q_file_base, "popID.txt"), col.names = FALSE, row.names=FALSE)









################################ADMIXTURE plotting manually
# Read sample metadata and filter

#get ordered samples list that you just made via copy/paste in sheets
order_df <- read.table(file.path(output_dir, "order.tsv"), sep = "\t", header = FALSE, stringsAsFactors = FALSE)
colnames(order_df) <- c("sample", "region", "continent")
order_df
nrow(order_df)

#remove samples that ADMIXURE took out:
excluded_samples <- readLines(file.path(q_file_base, "low_missingness_rm_outlier_removed_samples.txt"))
order_df_filtered <- order_df[!(order_df$sample %in% excluded_samples), ]
#note that the above will print the original row numbering. verify this is not reality: 
nrow(order_df_filtered)
order_df_filtered












#------------------------------------------------------------------------------------------

q_files <- list.files(q_file_base, pattern = "^CLUMPAK\\.\\d+\\.Q$", full.names = TRUE)
k_values <- as.integer(str_extract(basename(q_files), "(?<=\\.)\\d+(?=\\.Q)"))

# palette (adapt dynamically if needed)
max_k <- max(k_values)
color_pal <- RColorBrewer::brewer.pal(min(12, max_k), "Paired")
continent_colors <- setNames(RColorBrewer::brewer.pal(n = length(unique(order_df_filtered$continent)), name = "Set2"),
  unique(order_df_filtered$continent))


admix_plot_list = list()

##### add handbrake to stop after certain number of K
max_k = 9
# Loop through each K and create plot
for (i in seq_along(q_files)) {
  q_file_path <- q_files[i]
  k <- k_values[i]
  if (k > max_k) next  # skip anything beyond K = 7
  
  # Read Q file
  q_file <- read.table(q_file_path, header = FALSE)
  q_file <- read.table(q_files[i][1], header = FALSE )
  q_file <- cbind(final_samples, q_file)
  
  # Merge with ordered metadata
  q_merged <- left_join(order_df_filtered, q_file, by = "sample")
  stopifnot(nrow(q_merged) == nrow(order_df_filtered))
  
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
      sample = factor(sample, levels = order_df_filtered$sample),
      popGroup = as.factor(popGroup)
    )
  kdf2
  
  # Label data
  label_df <- order_df_filtered %>%
    mutate(
      sample = factor(sample, levels = order_df_filtered$sample),
      label = paste0(sample, " (", region, ")")
    )
  label_df

  kdf2 <- kdf2 %>%
  mutate(sample = factor(sample, levels = order_df_filtered$sample)) %>%
    arrange(sample, popGroup)
    
  label_df <- label_df %>%
    mutate(sample = factor(sample, levels = order_df_filtered$sample)) %>%
    arrange(sample)

  # Colors for this K
  fill_colors <- setNames(color_pal[seq_len(k)], as.character(seq_len(k)))

  # Calculate index of sample where continent changes for vertical lines BETWEEN continents
  label_df <- label_df %>%
    mutate(sample_index = as.numeric(sample))

  dividers_df <- label_df %>%
    mutate(next_continent = lead(continent)) %>%
    filter(continent != next_continent) %>%
    transmute(sample = sample_index + 0.5, prop = 1)

# We now have x positions for continent breaks

  
  # Plot
  admix_final <- ggplot(kdf2, aes(x = sample, y = prop, fill = popGroup)) +
    geom_col(width = 1, color = "black") +
    geom_col(  # simple black bars between continents
      data = dividers_df,
      aes(x = sample, y = prop),
      fill = "white",
      width = 0.3,
      inherit.aes = FALSE
    ) +
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
  ggsave(output_file, plot = admix_final, width = 20, height = 5)
  message("Saved ADMIXTURE plot for K = ", k)
  
  admix_plot_list[[paste0("K", k)]] <- admix_final









#############################       make meta-plot  #######################################

# Function to remove geom_text layers (sample labels) from a plot
remove_geom_text <- function(plot) {
  text_layers <- which(sapply(plot$layers, function(l) inherits(l$geom, "GeomText")))
  if(length(text_layers) > 0){
    plot$layers <- plot$layers[-text_layers]
  }
  plot
}


ks <- as.numeric(sub("K", "", names(admix_plot_list)))
admix_plot_list <- admix_plot_list[order(ks)]

# Then your combining code:
n_plots <- length(admix_plot_list)

# Remove geom_text and legends from all but the last plot
stripped_plots <- lapply(admix_plot_list[1:(n_plots - 1)], function(p) {
  p %>%
    remove_geom_text() + 
    theme(
      legend.position = "none",
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      axis.text.x = element_blank(),
      plot.title = element_blank(),  # REMOVE individual titles here
      plot.margin = margin(5, 20, 5, 20)
    )
})

# Keep last plot as is (with labels and legend)
last_plot <- admix_plot_list[[n_plots]] +
  theme(plot.title = element_blank())

# Combine all plots vertically
combined_plot <- wrap_plots(c(stripped_plots, list(last_plot)), ncol = 1) +
  plot_annotation(
    title = "ADMIXTURE Plots for K = 3 to 11",
    theme = theme(plot.title = element_text(size = 18, hjust = 0.5, margin = margin(b = 20)))
  )

# Save the combined figure
ggsave(
  filename = file.path(output_dir, "admixplots_all_clean.png"),
  plot = combined_plot,
  width = 12,
  height = 18,
  limitsize = TRUE
)
}
