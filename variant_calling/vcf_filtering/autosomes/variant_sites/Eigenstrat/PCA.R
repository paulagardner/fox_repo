# Load libraries
library(tidyverse)
library(googlesheets4)  # for reading directly from Google Sheets

# Authenticate (only required once)
#gs4_auth()  # Uncomment and run if this is your first time

# Import the tidy eigenvector file
pca <- read_tsv("~/Downloads/Vvulpestidy.evec")

# Rename columns to PC1, PC2, etc.
colnames(pca)[2:ncol(pca)] <- paste0("PC", 1:(ncol(pca)-1))

# OPTIONAL: If your PCA file doesn't have 'sample' column named explicitly
# colnames(pca)[1] <- "sample"

# Import metadata from Google Sheets (replace with your actual sheet URL or ID)
sheet_url <- "https://docs.google.com/spreadsheets/d/1KBqsCfPiuyno90iRJ-O0wV0HIijFncstvcIpGnpGY0k/edit?gid=1838883534#gid=1838883534"
metadata <- read_sheet(sheet_url) %>%
  mutate(sample = as.character(sample))  # <-- Fix here


# Make sure column names match between files for joining
# Rename if needed (e.g., "SampleID" → "sample")
# metadata <- metadata %>% rename(sample = SampleID)

# Join metadata to PCA data
pca_meta <- left_join(pca, metadata, by = "sample")

# Plot the PCA, color by continent
library(ggplot2)
library(dplyr)

# Assuming pca_df is your PCA result with 'sample' as the sample names and PC1, PC2 as the principal components
pca_df <- pca_meta  # replace with the actual name of your dataframe

# Create the PCA plot
library(ggrepel)

# Create the PCA plot using geom_text_repel to avoid overlapping labels
pca_plot <- ggplot(pca_meta, aes(x = PC1, y = PC2, color = continent)) +
  geom_point(size = 3) +  # Adjust the size of the points
  geom_text_repel(aes(label = sample), size = 3, box.padding = 0.35, max.overlaps =23) +  # Repel overlapping labels
  theme_minimal(base_size = 15) +  # Increase base font size for the plot
  labs(title = "PCA Plot by Continent", x = "Principal Component 1", y = "Principal Component 2") +
  theme(
    plot.title = element_text(size = 20, hjust = 0.5),  # Increase title size
    axis.title = element_text(size = 15),  # Increase axis title size
    axis.text = element_text(size = 12),  # Increase axis label size
    legend.text = element_text(size = 12),  # Increase legend text size
    legend.title = element_text(size = 14)  # Increase legend title size
  )

pca_plot

ggsave("~/Downloads/PCA_plot_large.png", plot = pca_plot, width = 12, height = 10)  # Save the plot with larger dimensions
