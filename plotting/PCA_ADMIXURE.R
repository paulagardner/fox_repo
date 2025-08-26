# --------- LOAD LIBRARIES ----------
library(tidyverse) #version in redfox conda: 2.0.0
library(googlesheets4) # 1.1.1
library(ggrepel) # 0.9.6
library(RColorBrewer) #1.1.3
library(grid) #4.4.3
library(gtable) #0.3.6
library(dplyr) #2.1.4
library(patchwork) #1.3.0

#to run this via the cluster: conda activate redfox
#R
# run this line by line. eg, 


# ---------------------- IMPORT & PREP google sheets data ----------------------
###running with version: 
# Import metadata from Google Sheets
sheet_url <- "https://docs.google.com/spreadsheets/d/1KBqsCfPiuyno90iRJ-O0wV0HIijFncstvcIpGnpGY0k/edit?gid=1838883534#gid=1838883534"
gs4_auth(path = "/gpfs/bio/xrq24scu/fox_repo/plotting/googlesheetskey.json")
metadata <- read_sheet(sheet_url) %>%
  mutate(sample = as.character(sample))

# Specify continent order manually (adjust as needed)
continent_order <- c("Europe", "Africa", "WestAsia", "Asia", "Americas")
metadata$continent <- factor(metadata$continent, levels = continent_order)
metadata

#-----------------PCA--------------


# Load PCA eigenvector file
pca <- read_tsv("/gpfs/data/bergstrom/paula/fox_repo/variant_calling/Eigenstrat/Vvulpes_rmoutliers_tidy.evec")

# Rename columns to PC1, PC2, etc.
colnames(pca)[2:ncol(pca)] <- paste0("PC", 1:(ncol(pca) - 1))

# Join metadata to PCA data
pca_meta <- left_join(pca, metadata, by = "sample")

# Remove specific samples from PCA plot (based on ADMIXTURE filtering script). 
excluded_samples <- readLines("/gpfs/data/bergstrom/paula/fox_repo/variant_calling/admixture/low_missingness_rm_outlier_removed_samples.txt")
#excluded_samples <- c(excluded_samples, "YPI1082")
#excluded_samples

pca_meta_filtered <- pca_meta %>%
  filter(!sample %in% excluded_samples)
pca_meta_filtered

# Define color palette for continents
#continents <- unique(metadata$continent)
#continent_colors <- setNames(RColorBrewer::brewer.pal(n = length(continents), name = "Set2"), continents)
continent_colors <- setNames(RColorBrewer::brewer.pal(n = length(continent_order), name = "Set2"), continent_order)



# Read eigenvalues and calculate variance explained
eigenvals <- scan("/gpfs/data/bergstrom/paula/fox_repo/variant_calling/Eigenstrat/Vvulpes_rmoutliers.pca.eval")
variance_explained <- 100 * eigenvals / sum(eigenvals)

# Format PC axis labels with variance
x_lab <- paste0("PC1 (", round(variance_explained[1], 1), "% variance)")
y_lab <- paste0("PC2 (", round(variance_explained[2], 1), "% variance)")


# ---------------------- PCA PLOT ----------------------

# PCA plot 
pca_plot <- ggplot(pca_meta_filtered, aes(x = PC1, y = PC2, color = continent)) +
  geom_point(size = 3) +
  geom_text_repel(aes(label = paste0(sample, " (", region, ")")), size = 3.5, box.padding = 0.5, max.overlaps = 35) +
  theme_minimal(base_size = 15) +
  labs(title = "PCA Plot by Continent", x = x_lab, y = y_lab) +
  theme(
    plot.title = element_text(size = 20, hjust = 0.5),
    axis.title = element_text(size = 15),
    axis.text = element_text(size = 12),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14)
  ) +
  scale_color_manual(values = continent_colors) +
  coord_fixed()   # <-- this makes 1 unit in x = 1 unit in y, but doesn't seem
  #to be working for the ggsave portion. Related to your concern the axes are not scaled
  #in a way that makes variance components make sense

pca_plot

# Show and save
#print(pca_plot)
ggsave("/gpfs/bio/xrq24scu/fox_repo/plotting/PCA_plot_large.png", plot = pca_plot, width = 10, height = 10)
#########################
#poster version- no sample names, larger labels
pca_plot_poster <- ggplot(pca_meta_filtered, aes(x = PC1, y = PC2, color = continent)) +
#  plot
  geom_point(size = 3) +
  geom_text_repel(aes(label = region), size = 5, box.padding = 0.5, max.overlaps = 20) +
  theme_minimal(base_size = 15) +
  labs(title = "PCA Plot by Continent", x = x_lab, y = y_lab) +
  theme(
    plot.title = element_text(size = 20, hjust = 0.5),
    axis.title = element_text(size = 15),
    axis.text = element_text(size = 12),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14)
  ) +
  scale_color_manual(values = continent_colors)

pca_plot_poster
ggsave("/gpfs/bio/xrq24scu/fox_repo/plotting/poster_PCA.png", plot = pca_plot_poster, width = 12, height = 10)


##########PCA WITH YPI1082 STILL INCLUDED
# code to generate plot without removing YPI1082
# ---------------------- IMPORT & PREP ----------------------

# Load PCA eigenvector file
pca2 <- read_tsv("/gpfs/data/bergstrom/paula/fox_repo/variant_calling/Eigenstrat/Vvulpestidy.evec")

# Rename columns to PC1, PC2, etc.
colnames(pca2)[2:ncol(pca2)] <- paste0("PC", 1:(ncol(pca2) - 1))

# Import metadata from Google Sheets
sheet_url <- "https://docs.google.com/spreadsheets/d/1KBqsCfPiuyno90iRJ-O0wV0HIijFncstvcIpGnpGY0k/edit?gid=1838883534#gid=1838883534"
gs4_auth(path = "/gpfs/bio/xrq24scu/fox_repo/plotting/googlesheetskey.json")
metadata <- read_sheet(sheet_url) %>%
  mutate(sample = as.character(sample))

# Join metadata to pca2 data
pca2_meta <- left_join(pca2, metadata, by = "sample")

# Read eigenvalues and calculate variance explained
eigenvals2 <- scan("/gpfs/data/bergstrom/paula/fox_repo/variant_calling/Eigenstrat/Vvulpes.pca.eval")
variance_explained2 <- 100 * eigenvals / sum(eigenvals)

# Format PC axis labels with variance
x_lab <- paste0("PC1 (", round(variance_explained2[1], 1), "% variance)")
y_lab <- paste0("PC2 (", round(variance_explained2[2], 1), "% variance)")


# ---------------------- PCA PLOT ----------------------

# PCA plot
pca2_plot <- ggplot(pca2_meta, aes(x = PC1, y = PC2, color = continent)) +
  geom_point(size = 3) +
  geom_text_repel(aes(label = paste0(sample, " (", region, ")")), size = 4, box.padding = 0.5, max.overlaps = 35) +
  theme_minimal(base_size = 15) +
  labs(title = "PCA Plot by Continent", x = x_lab, y = y_lab) +
  theme(
    plot.title = element_text(size = 20, hjust = 0.5),
    axis.title = element_text(size = 15),
    axis.text = element_text(size = 12),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14)
  ) +
  scale_color_manual(values = continent_colors)

# Show and save
#print(pca_plot)
ggsave("/gpfs/bio/xrq24scu/fox_repo/plotting/PCA_plot_large_withYPI1082.png", plot = pca2_plot, width = 12, height = 10)



###################################### CLUMPAK/ADMIXTURE Plotting #####################################
# This section generates CLUMPAK input files and produces ADMIXTURE barplots for different K values.
# It expects:
#   - Q files named CLUMPAK.<K>.Q in q_file_base
#   - A .fam file with sample order matching the Q files
#   - An order.tsv file with sample, region, and continent columns for plotting order (made manually to suit an 
      # order you thought was sensible)
#   - A metadata object with sample, region, and continent columns
#   - A list of excluded samples (to be filtered out of plots)
        ###NOTE HERE: If samples are removed because they are

# ----------- Set file paths -----------
q_file_base <- "/gpfs/data/bergstrom/paula/fox_repo/variant_calling/admixture/"
output_dir <- "/gpfs/data/bergstrom/paula/fox_repo/plotting/"

# ----------- Read sample order from .fam file -----------
# The .fam file's second column contains sample IDs, which matches the row order in Q files.
final_samples <- read.table(file.path(q_file_base, "low_missingness_rm_outlier_snpmanualfilter.fam"), header = FALSE)[,2]
final_samples <- tibble(sample = final_samples)

# ----------- Match metadata order to Q file sample order -----------
# Ensures that metadata rows are in the same order as the Q file rows.
qfileorder_df <- metadata %>%
  filter(sample %in% final_samples$sample) %>%
  slice(match(final_samples$sample, sample))

# ----------- Write numeric continent codes for CLUMPAK -----------
# CLUMPAK expects population codes as integers.
#popID <- as.numeric(factor(qfileorder_df$continent))
write.table(popID, file = file.path(q_file_base, "popID.txt"), col.names = FALSE, row.names = FALSE)

# ----------- Read plotting order and filter excluded samples -----------
# order.tsv should have columns: sample, region, continent
#order_df <- read.table(file.path(output_dir, "order.tsv"), sep = "\t", header = FALSE,
#                       col.names = c("sample", "region", "continent"), stringsAsFactors = FALSE)
#excluded_samples <- readLines(file.path(q_file_base, "low_missingness_rm_outlier_removed_samples.txt"))
#order_df_filtered <- order_df %>% filter(!sample %in% excluded_samples)
#order_df_filtered <- order_df_filtered %>%
#  mutate(population = as.integer(factor(region, levels = region_order)))
#write.table(order_df_filtered$populationID, file = file.path(q_file_base, "popID.txt"), col.names = FALSE, row.names = FALSE)

order_df <- read.table(file.path(output_dir, "order.tsv"), sep = "\t", header = FALSE,
                       col.names = c("sample", "region", "continent"), stringsAsFactors = FALSE)

excluded_samples <- readLines(file.path(q_file_base, "low_missingness_rm_outlier_removed_samples.txt"))

order_df_filtered <- order_df %>% 
  filter(!sample %in% excluded_samples) %>%
  mutate(population_ID = as.integer(factor(continent)))

# Write the population codes
write.table(order_df_filtered$population_ID,
            file = file.path(q_file_base, "popID.txt"),
            col.names = FALSE, row.names = FALSE, quote = FALSE)


###########              THIS NEEDS AN OVERHAUL  #########################
#probably, you should be filtering out samples with too much missing data before the q files stage,
#maybe??? you produce the removal text file at the qfile stage but double check if you shouldn't do 
#that earlier. the points with a lot off missing data behave well in the PCA, but I am just a bit
#uncomfortable with the organization of this pipeline. Basically overhaul the filtering/
#subsetting process so you can do it again on sex chromosomes, etc



# ----------- Set up color palettes -----------
max_k <- 6  # Only plot up to K=6
color_pal <- brewer.pal(min(12, max_k), "Paired")  # Cluster colors
continent_order <- unique(order_df_filtered$continent)
continent_colors <- setNames(brewer.pal(length(continent_order), "Set2"), continent_order)  # Continent colors

# ----------- Find Q files and extract K values -----------
q_files <- list.files(q_file_base, pattern = "^CLUMPAK\\.\\d+\\.Q$", full.names = TRUE)
k_values <- as.integer(str_extract(basename(q_files), "(?<=\\.)\\d+(?=\\.Q)"))

admix_plot_list <- list()  # Store plots for each K

# ----------- Loop over Q files and generate plots -----------
for (i in seq_along(q_files)) {
  k <- k_values[i]
  if (k > max_k) next  # Skip K values above max_k

  # Read Q file and add sample names as a column
  q_data <- read.table(q_files[i], header = FALSE)
  q_data <- bind_cols(final_samples, q_data)

  # Merge Q data with plotting order (order_df_filtered)
  # Ensures samples are plotted in the desired order
  q_merged <- left_join(order_df_filtered, q_data, by = "sample")
  stopifnot(nrow(q_merged) == nrow(order_df_filtered))  # Sanity check

  # Convert to long format for ggplot (one row per sample per cluster)
  kdf2 <- q_merged %>%
    pivot_longer(cols = starts_with("V"), names_to = "popGroup", values_to = "prop") %>%
    mutate(
      popGroup = as.factor(as.integer(str_remove(popGroup, "V"))),  # Cluster number as factor
      sample = factor(sample, levels = order_df_filtered$sample)    # Preserve plotting order
    )

  # Prepare label and divider positions for plot
  label_df <- order_df_filtered %>%
    mutate(sample = factor(sample, levels = order_df_filtered$sample),
           sample_index = as.numeric(sample))

  # Divider between continents
  dividers_df <- label_df %>%
    mutate(next_continent = lead(continent)) %>%
    filter(continent != next_continent) %>%
    transmute(sample = sample_index + 0.5, prop = 1)

  # ----------- Generate ADMIXTURE barplot -----------
  fill_colors <- setNames(color_pal[seq_len(k)], as.character(seq_len(k)))
  admix_final <- ggplot(kdf2, aes(x = sample, y = prop, fill = popGroup)) +
    geom_col(width = 1, color = "black") +
    # White divider between continents
    geom_col(data = dividers_df, aes(x = sample, y = prop), fill = "white", width = 0.3, inherit.aes = FALSE) +
    # Region labels below bars, colored by continent
    geom_text(data = label_df, aes(x = sample, y = -0.01, label = region, color = continent),
              angle = 45, hjust = 1, size = 6, inherit.aes = FALSE) +
    scale_fill_manual(values = fill_colors) +
    scale_color_manual(values = continent_colors) +
    coord_cartesian(clip = "off", ylim = c(-0.05, 1)) +
    theme_minimal() +
    labs(title = paste("ADMIXTURE Plot K =", k),
         x = "Samples", y = "Ancestry Proportion", fill = "Cluster") +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.line.x = element_blank(),
      axis.text.y = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      plot.margin = margin(30, 20, 20, 20),
      plot.title = element_text(size = 20, hjust = 0.5, margin = margin(b = 80)),
      axis.title.x = element_text(size = 14, margin = margin(t = 80)),
      axis.title.y = element_text(size = 14, margin = margin(r = 80)),
      legend.position = "right"
    )

  # Save individual plot
  ggsave(file.path(output_dir, paste0("admixplot_YPI1082removed", k, ".png")),
         plot = admix_final, width = 20, height = 7.5, bg = "transparent")
  admix_plot_list[[paste0("K", k)]] <- admix_final
}

# ----------- Combine all K plots vertically into a single figure -----------
# Remove geom_text (region labels) from all but the last plot for clarity
remove_geom_text <- function(plot) {
  text_layers <- which(sapply(plot$layers, function(l) inherits(l$geom, "GeomText")))
  if (length(text_layers) > 0) plot$layers <- plot$layers[-text_layers]
  plot
}
ks <- as.numeric(sub("K", "", names(admix_plot_list)))
admix_plot_list <- admix_plot_list[order(ks)]
n_plots <- length(admix_plot_list)
stripped_plots <- lapply(admix_plot_list[1:(n_plots - 1)], function(p) {
  remove_geom_text(p) +
    theme(legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank(),
          axis.text.x = element_blank(), plot.title = element_blank(), plot.margin = margin(5, 20, 5, 20))
})
last_plot <- admix_plot_list[[n_plots]] +
  theme(plot.title = element_blank(), plot.margin = margin(5, 20, 5, 20))
combined_plot <- wrap_plots(c(stripped_plots, list(last_plot)), ncol = 1) +
  plot_annotation(title = "ADMIXTURE Plots",
                  theme = theme(plot.title = element_text(size = 18, hjust = 0.5, margin = margin(b = 20))))
# Save combined plot
ggsave(file.path(output_dir, "admixplots_all_clean.png"),
       plot = combined_plot, bg = 'transparent', width = 16, height = 9, limitsize = TRUE)


# ==================
#   HETEROZYGOSITY 
# ==================

heterozygosity_url <- "https://docs.google.com/spreadsheets/d/1KBqsCfPiuyno90iRJ-O0wV0HIijFncstvcIpGnpGY0k/edit?gid=1838883534#gid=1838883534"
gs4_auth(path = "/gpfs/bio/xrq24scu/fox_repo/plotting/googlesheetskey.json")
heterozygosity_data <- read_sheet(heterozygosity_url, sheet = "HeterozygosityData") %>%
  mutate(sample = as.character(sample))

#get heterozygosity data merged w/ order_df_filtered (order.tsv with poor quality
#samples removed)
coverage_heterozygosity <- left_join(metadata, heterozygosity_data, by = "sample")

#Discard all data w/ a mean coverage below 15x, as indicated by which samples had 
#elevated heterozygosity (in your google plot)
filtered_cov_het <- coverage_heterozygosity[coverage_heterozygosity$mean_coverage > 10,]
filtered_cov_het

#now get this into a plotting order you like
merged.df<- left_join(order_df_filtered,filtered_cov_het, by = "sample")
merged.df

# Create a combined label for x-axis: "sample (region)"
merged.df$sample_region <- paste0(merged.df$sample, " (", merged.df$region.y, ")")
merged.df$sample_region <- factor(merged.df$sample_region, levels = merged.df$sample_region)

heterozygosity_plot <- ggplot(merged.df, aes(x = sample_region, y = Heterozygosity)) +
  geom_point(size = 3, color = "blue") +
  #plot 2* the error bars 
  geom_errorbar(aes(ymin = Heterozygosity - 2 * Jackknife_SE, ymax = Heterozygosity + 2 * Jackknife_SE),
                width = 0.2, color = "black") + 
  theme_minimal() +
  labs(title = "Heterozygosity by Sample", x = "Sample (Region)", y = "Heterozygosity") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("/gpfs/bio/xrq24scu/fox_repo/plotting/heterozygosity_plot.png", plot = heterozygosity_plot, width = 12, height = 6)



###### do it without statham samples












# =====================================================================================
# MSMC
# =====================================================================================
# make population codes for MSMC to run these as population groupings
# Generate population IDs from region, using the numbers we generated for CLUMPAK
write.table(
  order_df_filtered %>% select(sample, population),
  file = file.path(MSMC_dir, "pop_groups.txt"),
  col.names = FALSE, row.names = FALSE, quote = FALSE
)
#I THINK YOU MAY NEED TO RE-RUN THIS GROUPING NORTH AFRICA W/MIDDLE EAST


################now read in MSMC results

#######update this so that you know which pops are which, and color-code them accordingly.
#smartest thing would be to name the created MSMC files according to 

MSMC_dir <- "/gpfs/data/bergstrom/paula/fox_repo/MSMC/"
mu  <- 4.5e-9   # mutation rate from gray wolves, as used in the montane foxes 2024 paper
gen <- 2         # generation time in years

library(ggplot2)

# Read all msmc.final.txt files into one data frame
msmc_files <- list.files(MSMC_dir, pattern = "^pop[0-9]+\\.msmc\\.final\\.txt$", full.names = TRUE)


msmc <- do.call(rbind, lapply(msmc_files, function(f) {
  df <- read.table(f, header = TRUE)
  # get just the number after "pop"
  df$population_ID <- as.integer(sub(".*pop([0-9]+)\\.msmc\\.final\\.txt$", "\\1", basename(f)))
  df
}))

# make a lookup table: one continent per population_ID
pop_lookup <- unique(order_df_filtered[, c("population_ID", "continent")])

# merge onto msmc
msmc_annotated <- merge(msmc, pop_lookup, by = "population_ID", all.x = TRUE)


# Plot
MSMC_PLOT <- ggplot(msmc_annotated, aes(
  x = (left_time_boundary / mu) * gen,
  y = (1 / lambda_00) / (2 * mu),
  color = continent
)) +
  geom_step() +
  scale_x_log10(limits = c(15000, 1000000), breaks = scales::trans_breaks("log10", function(x) 10^x),
  labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  scale_y_log10(limits = c(10000, 500000), breaks = scales::trans_breaks("log10", function(x) 10^x),
  labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  theme_classic() +
  annotation_logticks() + # Default: log ticks on bottom and left
  labs (
    x = "Time (years)",
    y = "Effective population size (Ne)"
  )

ggsave(file.path(MSMC_dir, "msmc_plot.png"), MSMC_PLOT, width = 10, height = 6, dpi = 300)







############ testing out LARGE PLOT 
#  Required libraries
# Install if needed
library(ggplot2)
library(grid)
library(jpeg)
library(httr)

url <- "https://www.equal-earth.com/physical/Equal-Earth-Physical-Map-Raster-150E.jpg"
destfile <- "/gpfs/bio/xrq24scu/fox_repo/plotting/EqualEarth_150E.jpg"

if (!file.exists(destfile)) {
  GET(url, user_agent("Mozilla/5.0"), write_disk(destfile, overwrite = TRUE))
}

img <- readJPEG(destfile)
bg <- rasterGrob(img, width = unit(1,"npc"), height = unit(1,"npc"), interpolate = TRUE)

df <- data.frame(
  lon = c(139.6917, -77.0369, 151.2093),
  lat = c(35.6895, 38.9072, -33.8688)
)

eqearth_project <- function(lon, lat, lon0 = 150) {
  lon <- (lon - lon0) * pi/180
  lat <- lat * pi/180
  
  A1 <- 1.340264; A2 <- -0.081106; A3 <- 0.000893; A4 <- 0.003796
  theta <- asin(sqrt(3)/2 * sin(lat))
  
  x <- 2 * sqrt(3) / (3 * sqrt(pi)) * lon * cos(theta)
  y <- A1*theta + A2*sin(2*theta) + A3*sin(4*theta) + A4*sin(6*theta)
  
  data.frame(x = x, y = y)
}

proj <- eqearth_project(df$lon, df$lat, lon0 = 150)

# Rescale x/y to match annotation_custom extents
xlim <- c(-1, 1)
ylim <- c(-0.5, 0.5)

df$x <- (proj$x - min(proj$x)) / (max(proj$x) - min(proj$x)) * diff(xlim) + xlim[1]
df$y <- (proj$y - min(proj$y)) / (max(proj$y) - min(proj$y)) * diff(ylim) + ylim[1]


plot <- ggplot() +
  annotation_custom(bg, -1, 1, -0.5, 0.5) +
  geom_point(data = df, aes(x = x, y = y), color = "red", size = 3) +
  coord_fixed(xlim = c(-1, 1), ylim = c(-0.5, 0.5), expand = FALSE) +
  theme_void()
plot

# Save as PNG
ggsave(
  "/gpfs/bio/xrq24scu/fox_repo/plotting/EqualEarth_map.png",
  width = 12, height = 6, dpi = 300
)

# after: https://liamdbailey.com/posts/2024_08_31-azimuthalproj/
########################
#wget -O GGV_north.png -A "Mozilla/5.0" "https://upload.wikimedia.org/wikipedia/commons/7/7c/Gott-Goldberg-Vanderbei_Projection.png"

# ==============================
# Lambert Azimuthal Equal-Area Map (Pacific-centered)
# ==============================

# Load required packages
#had to conda install into the env:
# conda install -c conda-forge r-sf r-terra r-units r-s2 r-rnaturalearth r-rnaturalearthdata

# Example: 10m Physical raster (land/ocean)
library(sf)
library(ggplot2)
library(s2)

sf_use_s2(TRUE)  # enable spherical geometry

# Load countries
countries <- st_as_sf(s2_data_countries())

# Define hemisphere center (North Pole example)
lat0 <- 90
lon0 <- 0
## Create a proj string to describe the projection
## The projection alias is 'laea'
## The variables lat_0 and lon_0 represent our central point
crs <- sprintf("+proj=laea +lat_0=%f +lon_0=%f", lat0, lon0)

## Project our objects to this new projection
## Convert them into sf objects for easily plotting in ggplot
countries_laea <- st_transform(st_as_sfc(countries),
                               crs)

###clip to one hemisphere
buffer <- s2_buffer_cells(
  #central point
  as_s2_geography(sprintf("POINT(%f %f)", lon0, lat0)), 
  ## 9800km is the radius of the buffer to show one whole hemisphere of the earth
  distance = 9800000,
   ## How smooth is the buffer? Higher will be smoother (default: 1000)
  max_cells = 5000
)

countries_clip <- s2_intersection(buffer, countries)
## Reproject to laea
countries_clip_laea <- st_transform(st_as_sfc(countries_clip),
                                    crs)


## Create a plot
plot <- ggplot() +
  geom_sf(data = countries_clip_laea,
          fill = "grey50", colour = "grey10") +
  theme_void()

ggsave("/gpfs/bio/xrq24scu/fox_repo/plotting/LAEA_countries_hemisphere.png", plot, width = 8, height = 8, dpi = 300)


####################AFTER: https://milospopovic.net/display-raster-on-globe-in-r/

# libraries we need
libs <- c(
    "tidyverse", "sf", "giscoR",
    "mapview", "terra", "terrainr",
    "magick"
)
library(giscoR)
library(mapview)
library(terra)
library(terrainr)

# install missing libraries
installed_libs <- libs %in% rownames(installed.packages())
if (any(installed_libs == F)) {
    install.packages(libs[!installed_libs])
}
# load libraries
invisible(lapply(libs, library, character.only = T))



##################################################################
####Define where you want the map to be centered on for the azimuthal/orthographic 
#projection
center_lat <- 90.0
center_long <- 0.0

##equal area projection
equal_lat <- 0.0
equal_long <- 150

#First, we define two projections. One of them is the standard longitude/latitude projection,
# and the other is the Azimuthal orthographic projection. We use the former to project our map
# to WGS84 and then apply the latter to turn the flat Earth into a globe.
# define projections for terra
longlat_crs <- "+proj=longlat +datum=WGS84 +no_defs"
ortho_crs <- sprintf('+proj=ortho +lat_0=%f +lon_0=%f +x_0=0 +y_0=0 +R=6371000 +units=m +no_defs +type=crs', center_lat, center_long)
#define a third projection for the equal earth 
eqearth_crs <- sprintf('+proj=eqearth +lat_0=%f +lon_0=%f +x_0=0 +y_0=0 +R=6371000 +units=m +no_defs +type=crs', equal_lat, equal_long)


#In the next step, we retrieve the world polygon. We will use this polygon to crop the nightlight data.
get_flat_world_sf <- function() {
    world <- giscoR::gisco_get_countries(
        year = "2016",
        epsg = "4326",
        resolution = "10"
    ) %>%
        sf::st_transform(longlat_crs)

    world_vect <- terra::vect(world)

    return(world_vect)
}

world_vect <- get_flat_world_sf()


#Here we opt for the latest 2017 data in the lowest possible resolution 3600x1800 in order to make our life easier. 
#Thanks to terra package, we can import the file directly into R. Next, we crop the raster file so that it fits the world map coordinates.
# 2. NASA DATA
#-------------
get_nasa_data <- function() {
    ras <- terra::rast("/vsicurl/https://eoimages.gsfc.nasa.gov/images/imagerecords/144000/144898/BlackMarble_2016_01deg_geo.tif")
    rascrop <- terra::crop(x = ras, y = world_vect, snap = "in")
    ras_latlong <- terra::project(rascrop, longlat_crs)
    ras_ortho <- terra::project(ras_latlong, ortho_crs)
    return(ras_ortho)
}

ras_ortho <- get_nasa_data()

png("/gpfs/bio/xrq24scu/fox_repo/plotting/tutorial.png", width = 8*300, height = 8*300, res = 300)
plot(ras_ortho)
dev.off()




# Norwich coordinates (lon, lat)
norwich <- st_sf(
  name = "Norwich",
  geometry = st_sfc(st_point(c(1.2966, 52.6309)), crs = 4326)
)
# Project Norwich to the same CRS as your raster
norwich_ortho <- st_transform(norwich, crs = ortho_crs)
# Open PNG device
png("/gpfs/bio/xrq24scu/fox_repo/plotting/tutorial_norwich.png",
    width = 8*300, height = 8*300, res = 300)

# Plot the NASA raster in orthographic projection
plot(ras_ortho)

# Overlay Norwich as a red point
points(
  x = st_coordinates(norwich_ortho)[,1],
  y = st_coordinates(norwich_ortho)[,2],
  col = "red", pch = 19, cex = 1.5
)

# Optional: add a label
text(
  x = st_coordinates(norwich_ortho)[,1],
  y = st_coordinates(norwich_ortho)[,2],
  labels = "Norwich",
  pos = 4, col = "white", cex = 1
)

dev.off()

#Next, we project the satellite imagery into the longlat projection and retrieve a flat world map