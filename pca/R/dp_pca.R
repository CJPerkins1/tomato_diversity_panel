#pca on concordant SNPs of the diversity panel

library(SNPRelate)
setwd("/Users/cperkins/Desktop/palanivelu_lab/diversity_panel/dp_vcfs/plink")
# Convert PLINK to GDS format
snpgdsBED2GDS(
  bed.fn = "/Users/cperkins/Desktop/palanivelu_lab/diversity_panel/dp_vcfs/plink/data/concordant_popgen.bed",
  bim.fn = "/Users/cperkins/Desktop/palanivelu_lab/diversity_panel/dp_vcfs/plink/data/concordant_popgen.bim",
  fam.fn = "/Users/cperkins/Desktop/palanivelu_lab/diversity_panel/dp_vcfs/plink/data/concordant_popgen.fam",
  out.gdsfn = "concordant_popgen.gds"
)

# Open GDS file
genofile <- snpgdsOpen("concordant_popgen.gds")

# Run PCA
pca <- snpgdsPCA(genofile, num.thread = 8, autosome.only = FALSE)

# Close file
snpgdsClose(genofile)

# Extract results
pc_scores <- data.frame(
  sample.id = pca$sample.id,
  PC1 = pca$eigenvect[,1],
  PC2 = pca$eigenvect[,2],
  PC3 = pca$eigenvect[,3],
  PC4 = pca$eigenvect[,4],
  PC5 = pca$eigenvect[,5]
)

# Variance explained
variance <- pca$varprop * 100

# Basic plot
plot(pc_scores$PC1, pc_scores$PC2,
     xlab = paste0("PC1 (", round(variance[1], 2), "%)"),
     ylab = paste0("PC2 (", round(variance[2], 2), "%)"),
     main = "PCA of Population Structure",
     pch = 19, col = "steelblue")

# Scree plot
barplot(variance[1:10], 
        names.arg = 1:10,
        xlab = "Principal Component",
        ylab = "Variance Explained (%)",
        main = "Scree Plot")

dp_sheet <- read.csv("/Users/cperkins/Desktop/palanivelu_lab/diversity_panel/dp_vcfs/plink/data/dp_sheet.csv")

library(dplyr)
library(ggplot2)

# Extract the SRX ID from the sample.id column in pc_scores
pc_scores <- pc_scores %>%
  mutate(Experiment = sub("_.*", "", sample.id))  # Extract everything before the first underscore

# Merge with your metadata
pca_data <- pc_scores %>%
  left_join(dp_sheet, by = c("Experiment" = "Experiment"))

# Check the merge
head(pca_data)
table(pca_data$read_type)  # Should show paired/single counts

# Calculate variance explained
variance <- pca$varprop * 100  # If using SNPRelate
# OR if reading from PLINK output:
# eigenvalues <- scan("concordant_popgen_pca.eigenval")
# variance <- (eigenvalues / sum(eigenvalues)) * 100

# Remove Heinz reference
pca_data_no_ref <- pca_data %>%
  filter(name_CW != "CW0000")

# Plot PC1 vs PC2 colored by read type
ggplot(pca_data_no_ref, aes(x = PC1, y = PC2, color = read_type)) +
  geom_point(size = 3, alpha = 1) +
  scale_color_manual(values = c("paired" = "#FDA502", "single" = "#800680"),
                     name = "Read Type") +
  labs(
    x = paste0("PC1 (", round(variance[1], 2), "%)"),
    y = paste0("PC2 (", round(variance[2], 2), "%)"),
    title = "PCA of Population Structure by Read Type"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "right"
  )

#Adding cerasaform
library(googlesheets4)
sup_table <- read_sheet("1XUUZ9Uj6hJweuWM0BhO6vs88dXrs-o0V80yWVsH-NSA", 2)


# Extract genotype info from sup_table
genotype_info <- sup_table %>%
  select(`CW ID`, Genotype) %>%
  mutate(
    variety = case_when(
      grepl("SLL", Genotype) ~ "lycopersicum var. lycopersicum",
      grepl("SLC", Genotype) ~ "lycopersicum var. cerasiforme",
      TRUE ~ NA_character_
    )
  ) %>%
  rename(name_CW = `CW ID`)

# Merge with pca_data
pca_data <- pca_data %>%
  left_join(genotype_info %>% select(name_CW, variety), by = "name_CW")

# Create a combined species/variety column for plotting
pca_data <- pca_data %>%
  mutate(
    species_label = case_when(
      species == "lycopersicum" & !is.na(variety) ~ variety,
      species == "lycopersicum" & is.na(variety) ~ "lycopersicum var. lycopersicum",  # Default remaining to var. l
      !is.na(species) ~ species,
      TRUE ~ NA_character_
    )
  )

pca_data <- pca_data %>%
  select(-variety.x) %>%
  rename(variety = variety.y)

# Check the distribution
table(pca_data$species_label, useNA = "ifany")

# Set factor levels for legend order
pca_data$species_label <- factor(pca_data$species_label, 
                                 levels = c("lycopersicum var. lycopersicum",
                                            "lycopersicum var. cerasiforme",
                                            "pimpinellifolium", 
                                            "galapagense", 
                                            "cheesmaniae"))

# Create 5 colors from viridis (2 purples for lycopersicum varieties, then progression)
viridis_colors <- c(
  viridis(6, option = "D")[1],  # darkest purple for var. lycopersicum
  viridis(6, option = "D")[3],  # lighter purple for var. cerasiforme
  viridis(6, option = "D")[4],  # blue-green for pimpinellifolium
  viridis(6, option = "D")[5],  # green-yellow for galapagense
  viridis(6, option = "D")[6]   # yellow for cheesmaniae
)

# Plot
ggplot(pca_data_no_ref, aes(x = PC1, y = PC2, color = species_label)) +
  geom_point(size = 3, alpha = 1) +
  scale_color_manual(
    name = "Species",
    values = viridis_colors,
    na.value = "gray50"
  ) +
  labs(
    x = paste0("PC1 (", round(variance[1], 2), "%)"),
    y = paste0("PC2 (", round(variance[2], 2), "%)"),
    title = "PCA - Population Structure by Species"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "right"
  )

#Add in species
library(googlesheets4)

# Loading the data --------------------------------------------------------
accessions <- read.csv("/Users/cperkins/Desktop/palanivelu_lab/diversity_panel/data/accessions.csv")

# Safe merge - only keep name_CW and species, match to PCA samples
pca_data <- pca_data %>%
  left_join(accessions %>% select(name_CW, species), 
            by = "name_CW")

# pca_data <- pca_data %>%
#   select(-species.x) %>%
#   rename(species = species.y)

# Check the merge
cat("Total PCA samples:", nrow(pca_data), "\n")
cat("Samples with species info:", sum(!is.na(pca_data$species)), "\n")
cat("Samples missing species info:", sum(is.na(pca_data$species)), "\n\n")

# See what species you have
table(pca_data$species, useNA = "ifany")

# Check confounding with read type
cat("\nSpecies vs Read Type:\n")
table(pca_data$species, pca_data$read_type, useNA = "ifany")

# Plot PCA by species
library(viridis)

# Set species as ordered factor
pca_data$species <- factor(pca_data$species, 
                           levels = c("lycopersicum", "pimpinellifolium", 
                                      "galapagense", "cheesmaniae"))

# Get viridis colors (purple to yellow)
viridis_colors <- viridis(4, option = "D")  # D is the default viridis palette

# Plot with custom order and colors
ggplot(pca_data, aes(x = PC1, y = PC2, color = species)) +
  geom_point(size = 3, alpha = 1) +
  scale_color_manual(
    name = "Species",
    values = viridis_colors,
    na.value = "gray50",
    breaks = c("lycopersicum", "pimpinellifolium", "galapagense", "cheesmaniae")
  ) +
  labs(
    x = paste0("PC1 (", round(variance[1], 2), "%)"),
    y = paste0("PC2 (", round(variance[2], 2), "%)"),
    title = "PCA - Population Structure by Species"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "right"
  )


#now in 3D
library(plotly)

# Create 3D plot with custom camera angle
plot_ly(pca_data, 
        x = ~PC1, y = ~PC2, z = ~PC3,
        color = ~species_label,
        colors = viridis_colors,
        type = "scatter3d", 
        mode = "markers",
        marker = list(size = 5),
        text = ~paste("Sample:", sample.id, "<br>Species:", species_label),
        hoverinfo = "text") %>%
  layout(
    scene = list(
      xaxis = list(title = paste0("PC1 (", round(variance[1], 2), "%)")),
      yaxis = list(title = paste0("PC2 (", round(variance[2], 2), "%)")),
      zaxis = list(title = paste0("PC3 (", round(variance[3], 2), "%)")),
      camera = list(
        eye = list(x = 1, y = 1, z = 1.3),  # Camera position
        center = list(x = 0, y = 0, z = 0),     # What to look at
        up = list(x = 0, y = 0, z = 1)          # Up direction
      )
    ),
    title = "3D PCA - Population Structure by Species",
    legend = list(title = list(text = "Species"))
  )

#extracting PCs for corrected LD decay calculationsL

# Create sample lists for each species group
library(dplyr)

# 1. Lycopersicum (both varieties combined)
lycopersicum_samples <- pca_data %>%
  filter(species == "lycopersicum") %>%
  select(sample.id) %>%
  write.table("lycopersicum_samples.txt", 
              quote = FALSE, row.names = FALSE, col.names = FALSE)

# 2. Pimpinellifolium only
pimpi_samples <- pca_data %>%
  filter(species == "pimpinellifolium") %>%
  select(sample.id) %>%
  write.table("pimpinellifolium_samples.txt", 
              quote = FALSE, row.names = FALSE, col.names = FALSE)

# 3. All cultivated (lycopersicum + cerasiforme)
cultivated_samples <- pca_data %>%
  filter(species_label %in% c("lycopersicum var. lycopersicum", 
                              "lycopersicum var. cerasiforme")) %>%
  select(sample.id) %>%
  write.table("cultivated_samples.txt", 
              quote = FALSE, row.names = FALSE, col.names = FALSE)




#Stratified analysis: -------------------------------------------------

library(dplyr)
library(ggplot2)
library(data.table)
library(viridis)

setwd("/Users/cperkins/Desktop/palanivelu_lab/diversity_panel/LD/results/stratified")

# ============================================================================
# PART 1: Calculate LD Decay Statistics
# ============================================================================

# Function to calculate LD decay for one group/chromosome
calculate_ld_decay <- function(ld_file, group, chr) {
  cat("Processing:", ld_file, "\n")
  
  # Read LD file
  ld <- fread(ld_file)
  
  # Calculate distance between SNPs
  ld <- ld %>%
    mutate(
      distance_bp = abs(BP_B - BP_A),
      distance_kb = distance_bp / 1000
    )
  
  # Bin distances (every 10kb up to 1000kb)
  ld <- ld %>%
    mutate(
      distance_bin = cut(distance_kb, 
                         breaks = seq(0, 1000, by = 10),
                         labels = seq(5, 995, by = 10),  # Midpoint of each bin
                         include.lowest = TRUE)
    ) %>%
    filter(!is.na(distance_bin))
  
  # Calculate mean r2 per bin
  decay <- ld %>%
    group_by(distance_bin) %>%
    summarise(
      mean_r2 = mean(R2, na.rm = TRUE),
      median_r2 = median(R2, na.rm = TRUE),
      n_pairs = n(),
      .groups = "drop"
    ) %>%
    mutate(
      distance_kb = as.numeric(as.character(distance_bin)),
      group = group,
      chromosome = chr
    )
  
  return(decay)
}

# Function to find decay distance (where r2 drops below threshold)
find_decay_distance <- function(decay_data, threshold = 0.2) {
  # Find first distance where mean_r2 < threshold
  decay_point <- decay_data %>%
    filter(mean_r2 < threshold) %>%
    slice(1) %>%
    pull(distance_kb)
  
  if (length(decay_point) == 0) {
    return(NA)  # No decay within 1000kb
  } else {
    return(decay_point)
  }
}

# ============================================================================
# PART 2: Process All Files
# ============================================================================

# List all LD files
groups <- c("lycopersicum", "pimpinellifolium", "cultivated")
chromosomes <- paste0("SL4.0ch", sprintf("%02d", 1:12))

# Process all files
all_decay <- list()

for (group in groups) {
  for (chr in chromosomes) {
    file_name <- paste0("ld_", group, "_", chr, ".ld.gz")
    
    if (file.exists(file_name)) {
      decay <- calculate_ld_decay(file_name, group, chr)
      all_decay[[paste(group, chr, sep = "_")]] <- decay
    } else {
      cat("Warning: File not found:", file_name, "\n")
    }
  }
}

# Combine all results
decay_combined <- bind_rows(all_decay)

# ============================================================================
# PART 3: Calculate Summary Statistics
# ============================================================================

# Calculate decay distance for each group/chromosome
decay_summary <- decay_combined %>%
  group_by(group, chromosome) %>%
  summarise(
    decay_distance_kb = find_decay_distance(cur_data(), threshold = 0.2),
    mean_r2_100_500kb = mean(mean_r2[distance_kb >= 100 & distance_kb <= 500], na.rm = TRUE),
    mean_r2_0_100kb = mean(mean_r2[distance_kb <= 100], na.rm = TRUE),
    .groups = "drop"
  )

# Overall summary by group
group_summary <- decay_summary %>%
  group_by(group) %>%
  summarise(
    mean_decay_distance = mean(decay_distance_kb, na.rm = TRUE),
    median_decay_distance = median(decay_distance_kb, na.rm = TRUE),
    mean_r2_100_500kb = mean(mean_r2_100_500kb, na.rm = TRUE),
    n_chromosomes = n(),
    .groups = "drop"
  )

print("=== DECAY SUMMARY BY GROUP ===")
print(group_summary)

print("\n=== DECAY BY CHROMOSOME ===")
print(n=36, decay_summary %>% arrange(group, chromosome))

# ============================================================================
# PART 4: Visualization
# ============================================================================

# Plot 1: LD decay curves by group (genome-wide average)
decay_avg <- decay_combined %>%
  group_by(group, distance_kb) %>%
  summarise(
    mean_r2 = mean(mean_r2, na.rm = TRUE),
    se_r2 = sd(mean_r2, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  )

p1 <- ggplot(decay_avg, aes(x = distance_kb, y = mean_r2, color = group)) +
  geom_line(size = 1.5) +
  geom_ribbon(aes(ymin = mean_r2 - se_r2, ymax = mean_r2 + se_r2, fill = group), 
              alpha = 0.2, color = NA) +
  geom_hline(yintercept = 0.2, linetype = "dashed", color = "gray50") +
  annotate("text", x = 900, y = 0.22, label = "r² = 0.2 threshold", 
           color = "gray50", size = 3) +
  scale_color_viridis_d(option = "D", end = 0.8) +
  scale_fill_viridis_d(option = "D", end = 0.8) +
  labs(
    title = "LD Decay Comparison Across Species Groups",
    subtitle = "Genome-wide average with standard error",
    x = "Distance (kb)",
    y = "Mean r²",
    color = "Group",
    fill = "Group"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    legend.position = "right"
  )

print(p1)
ggsave("ld_decay_comparison.png", p1, width = 10, height = 6, dpi = 300)

# Plot 2: Faceted by chromosome
p2 <- ggplot(decay_combined, aes(x = distance_kb, y = mean_r2, color = group)) +
  geom_line(size = 0.8) +
  geom_hline(yintercept = 0.2, linetype = "dashed", color = "gray50", size = 0.5) +
  facet_wrap(~chromosome, ncol = 4) +
  scale_color_viridis_d(option = "D", end = 0.8) +
  labs(
    title = "LD Decay by Chromosome",
    x = "Distance (kb)",
    y = "Mean r²",
    color = "Group"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.position = "bottom",
    strip.text = element_text(face = "bold")
  )

print(p2)
ggsave("ld_decay_by_chromosome.png", p2, width = 14, height = 10, dpi = 300)

# Plot 3: Decay distance comparison
decay_summary_filtered <- decay_summary %>%
  filter(!is.na(decay_distance_kb))

p3 <- ggplot(decay_summary_filtered, aes(x = chromosome, y = decay_distance_kb, fill = group)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_viridis_d(option = "D", end = 0.8) +
  labs(
    title = "LD Decay Distance by Chromosome and Group",
    subtitle = "Distance where r² drops below 0.2",
    x = "Chromosome",
    y = "Decay Distance (kb)",
    fill = "Group"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right"
  )

print(p3)
ggsave("decay_distance_comparison.png", p3, width = 12, height = 6, dpi = 300)

# Plot 4: Mean r² in 100-500kb window
p4 <- ggplot(decay_summary, aes(x = chromosome, y = mean_r2_100_500kb, fill = group)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_viridis_d(option = "D", end = 0.8) +
  labs(
    title = "Mean r² in 100-500kb Window",
    subtitle = "Measure of long-range LD",
    x = "Chromosome",
    y = "Mean r² (100-500kb)",
    fill = "Group"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right"
  )

print(p4)
ggsave("mean_r2_100_500kb.png", p4, width = 12, height = 6, dpi = 300)

# ============================================================================
# PART 5: Statistical Comparisons
# ============================================================================

# Compare groups using t-test or Wilcoxon
cat("\n=== STATISTICAL COMPARISONS ===\n")

# Lycopersicum vs Pimpinellifolium
lyc_decay <- decay_summary %>% filter(group == "lycopersicum") %>% pull(decay_distance_kb)
pim_decay <- decay_summary %>% filter(group == "pimpinellifolium") %>% pull(decay_distance_kb)

if (sum(!is.na(lyc_decay)) > 0 & sum(!is.na(pim_decay)) > 0) {
  test_result <- wilcox.test(lyc_decay, pim_decay, na.action = na.omit)
  cat("\nLycopersicum vs Pimpinellifolium decay distance:\n")
  print(test_result)
}

# Export results
write.csv(decay_summary, "ld_decay_summary_by_chromosome.csv", row.names = FALSE)
write.csv(group_summary, "ld_decay_summary_by_group.csv", row.names = FALSE)
write.csv(decay_combined, "ld_decay_full_data.csv", row.names = FALSE)

cat("\n=== Analysis Complete! ===\n")
cat("Files saved:\n")
cat("  - ld_decay_comparison.png\n")
cat("  - ld_decay_by_chromosome.png\n")
cat("  - decay_distance_comparison.png\n")
cat("  - mean_r2_100_500kb.png\n")
cat("  - ld_decay_summary_by_chromosome.csv\n")
cat("  - ld_decay_summary_by_group.csv\n")



#Whole genome heatmap--------------------------------------------------------------------------

library(dplyr)
library(ggplot2)
library(data.table)
library(viridis)
library(tidyr)

# Load decay summary we already calculated
decay_summary <- read.csv("ld_decay_summary_by_chromosome.csv")

# ============================================================================
# Heatmap 1: Mean r² by chromosome and distance bin
# ============================================================================

# Read a subset of files to create distance-binned heatmap
# We'll use cultivated group as example, but can do all three

create_heatmap_data <- function(group_name) {
  chromosomes <- paste0("SL4.0ch", sprintf("%02d", 1:12))
  
  all_bins <- list()
  
  for (chr in chromosomes) {
    file_name <- paste0("ld_", group_name, "_", chr, ".ld.gz")
    
    if (file.exists(file_name)) {
      cat("Reading:", file_name, "\n")
      
      # Read and bin
      ld <- fread(file_name) %>%
        mutate(
          distance_kb = abs(BP_B - BP_A) / 1000,
          distance_bin = cut(distance_kb, 
                             breaks = c(0, 10, 25, 50, 100, 250, 500, 750, 1000),
                             labels = c("0-10", "10-25", "25-50", "50-100", 
                                        "100-250", "250-500", "500-750", "750-1000"))
        ) %>%
        filter(!is.na(distance_bin))
      
      # Calculate mean r2 per bin
      binned <- ld %>%
        group_by(distance_bin) %>%
        summarise(
          mean_r2 = mean(R2, na.rm = TRUE),
          .groups = "drop"
        ) %>%
        mutate(chromosome = chr)
      
      all_bins[[chr]] <- binned
    }
  }
  
  bind_rows(all_bins) %>%
    mutate(group = group_name)
}

# Create data for all three groups
heatmap_cultivated <- create_heatmap_data("cultivated")
heatmap_lycopersicum <- create_heatmap_data("lycopersicum")
heatmap_pimpinellifolium <- create_heatmap_data("pimpinellifolium")

# Combine
heatmap_data <- bind_rows(heatmap_cultivated, heatmap_lycopersicum, heatmap_pimpinellifolium)

# Reorder chromosomes properly
heatmap_data$chromosome <- factor(heatmap_data$chromosome, 
                                  levels = paste0("SL4.0ch", sprintf("%02d", 1:12)))

# Reorder distance bins
heatmap_data$distance_bin <- factor(heatmap_data$distance_bin,
                                    levels = c("0-10", "10-25", "25-50", "50-100", 
                                               "100-250", "250-500", "500-750", "750-1000"))

# Plot 1: Faceted heatmap by group
p_heatmap1 <- ggplot(heatmap_data, aes(x = distance_bin, y = chromosome, fill = mean_r2)) +
  geom_tile(color = "white", size = 0.5) +
  facet_wrap(~group, ncol = 1) +
  scale_fill_viridis_c(option = "plasma", 
                       limits = c(0, 0.6),
                       breaks = seq(0, 0.6, 0.1),
                       name = "Mean r²") +
  labs(
    title = "Genome-wide LD Heatmap by Species Group",
    subtitle = "Mean r² across distance bins",
    x = "Distance between SNPs (kb)",
    y = "Chromosome"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 16),
    plot.subtitle = element_text(hjust = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(face = "bold", size = 12),
    panel.grid = element_blank()
  )

print(p_heatmap1)
ggsave("ld_heatmap_by_group.png", p_heatmap1, width = 10, height = 10, dpi = 300)

# ============================================================================
# Heatmap 2: Side-by-side comparison (chromosomes vs groups)
# ============================================================================

# Pivot to wide format for easier comparison
heatmap_wide <- heatmap_data %>%
  select(chromosome, distance_bin, mean_r2, group) %>%
  pivot_wider(names_from = group, values_from = mean_r2, names_prefix = "r2_")

# Create difference heatmap: cultivated - pimpinellifolium
heatmap_diff <- heatmap_wide %>%
  mutate(ld_difference = r2_cultivated - r2_pimpinellifolium)

p_heatmap_diff <- ggplot(heatmap_diff, aes(x = distance_bin, y = chromosome, fill = ld_difference)) +
  geom_tile(color = "white", size = 0.5) +
  scale_fill_gradient2(
    low = "blue", mid = "white", high = "red",
    midpoint = 0,
    limits = c(-0.1, 0.4),
    name = "LD Difference\n(Cultivated - Wild)"
  ) +
  labs(
    title = "Domestication Impact on LD Across the Genome",
    subtitle = "Difference in mean r² between cultivated and wild tomatoes",
    x = "Distance between SNPs (kb)",
    y = "Chromosome"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
    plot.subtitle = element_text(hjust = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank()
  )

print(p_heatmap_diff)
ggsave("ld_heatmap_domestication_difference.png", p_heatmap_diff, width = 10, height = 8, dpi = 300)

# ============================================================================
# Heatmap 3: Simple chromosome x group mean r² (100-500kb window)
# ============================================================================

# Use the decay_summary data
p_heatmap_simple <- ggplot(decay_summary, 
                           aes(x = group, y = chromosome, fill = mean_r2_100_500kb)) +
  geom_tile(color = "white", size = 1) +
  geom_text(aes(label = round(mean_r2_100_500kb, 2)), color = "white", size = 4, fontface = "bold") +
  scale_fill_viridis_c(option = "plasma", 
                       limits = c(0.1, 0.5),
                       name = "Mean r²\n(100-500kb)") +
  labs(
    title = "Long-range LD Across the Genome",
    subtitle = "Mean r² in 100-500kb window",
    x = "Species Group",
    y = "Chromosome"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 16),
    plot.subtitle = element_text(hjust = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank()
  )

print(p_heatmap_simple)
ggsave("ld_heatmap_simple_100_500kb.png", p_heatmap_simple, width = 8, height = 8, dpi = 300)

# ============================================================================
# Heatmap 4: Decay distance heatmap
# ============================================================================

p_heatmap_decay <- ggplot(decay_summary, 
                          aes(x = group, y = chromosome, fill = decay_distance_kb)) +
  geom_tile(color = "white", size = 1) +
  geom_text(aes(label = ifelse(is.na(decay_distance_kb), ">1000", 
                               as.character(round(decay_distance_kb)))), 
            color = "white", size = 3.5, fontface = "bold") +
  scale_fill_viridis_c(option = "magma", 
                       na.value = "#440154",  # Darkest color for NA (>1000kb)
                       limits = c(0, 1000),
                       name = "Decay Distance\n(kb)") +
  labs(
    title = "LD Decay Distance Across the Genome",
    subtitle = "Distance where r² drops below 0.2 (>1000 = no decay)",
    x = "Species Group",
    y = "Chromosome"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 16),
    plot.subtitle = element_text(hjust = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank()
  )

print(p_heatmap_decay)
ggsave("ld_heatmap_decay_distance.png", p_heatmap_decay, width = 8, height = 8, dpi = 300)

# ============================================================================
# Summary stats for paper/presentation
# ============================================================================

cat("\n=== HEATMAP SUMMARY STATISTICS ===\n\n")

cat("Chr 11 r² values (100-500kb):\n")
chr11_stats <- decay_summary %>%
  filter(chromosome == "SL4.0ch11") %>%
  select(group, mean_r2_100_500kb)
print(chr11_stats)

cat("\nChr 9 r² values (100-500kb):\n")
chr9_stats <- decay_summary %>%
  filter(chromosome == "SL4.0ch09") %>%
  select(group, mean_r2_100_500kb)
print(chr9_stats)

cat("\nChromosomes with no decay in cultivated:\n")
no_decay <- decay_summary %>%
  filter(group == "cultivated", is.na(decay_distance_kb)) %>%
  pull(chromosome)
print(no_decay)

cat("\nMean difference in r² between cultivated and wild:\n")
mean_diff <- mean(heatmap_diff$ld_difference, na.rm = TRUE)
cat(sprintf("Average LD increase from domestication: %.3f\n", mean_diff))

cat("\n=== Files saved ===\n")
cat("  - ld_heatmap_by_group.png (main heatmap)\n")
cat("  - ld_heatmap_domestication_difference.png (cultivated - wild)\n")
cat("  - ld_heatmap_simple_100_500kb.png (simple grid with values)\n")
cat("  - ld_heatmap_decay_distance.png (decay distances)\n")

# 
# 
# 
# 
# library(SNPRelate)
# 
# # Reopen the GDS file
# genofile <- snpgdsOpen("concordant_popgen.gds")
# 
# # Run PCA with loadings - eigen.cnt needs to be numeric
# # Setting a higher number calculates loadings for more PCs
# pca_with_loadings <- snpgdsPCA(genofile, 
#                                num.thread = 8, 
#                                autosome.only = FALSE,
#                                eigen.cnt = 32)  # Number of eigenvectors to calculate
# 
# # Check what's available
# names(pca_with_loadings)
# str(pca_with_loadings)
# 
# # The loadings might be accessible through a different method
# # Let's see what we got
# print(head(pca_with_loadings$eigenvect))
# 
# snpgdsClose(genofile)
# 
# # Find top contributing SNPs for PC1
# pc1_loadings <- data.frame(
#   snp.id = pca$snp.id,
#   loading = loadings[, 1],
#   abs_loading = abs(loadings[, 1])
# ) %>%
#   arrange(desc(abs_loading))
# 
# head(pc1_loadings, 20)  # Top 20 SNPs contributing to PC1
# 
# # Top SNPs for PC2
# pc2_loadings <- data.frame(
#   snp.id = pca$snp.id,
#   loading = loadings[, 2],
#   abs_loading = abs(loadings[, 2])
# ) %>%
#   arrange(desc(abs_loading))
# 
# head(pc2_loadings, 20)
# 
# # Visualize loadings distribution
# par(mfrow = c(2, 2))
# hist(loadings[, 1], breaks = 50, main = "PC1 Loadings", xlab = "Loading")
# hist(loadings[, 2], breaks = 50, main = "PC2 Loadings", xlab = "Loading")
# hist(loadings[, 3], breaks = 50, main = "PC3 Loadings", xlab = "Loading")
# hist(loadings[, 4], breaks = 50, main = "PC4 Loadings", xlab = "Loading")
# 
# 
# 
# 
# 
# 
# # Multiple PC comparisons
# library(GGally)
# ggpairs(pca_data, 
#         columns = c("PC1", "PC2", "PC3", "PC4", "PC5"),
#         aes(color = read_type, alpha = 0.6),
#         upper = list(continuous = "points"),
#         lower = list(continuous = wrap("points", size = 1)),
#         diag = list(continuous = wrap("densityDiag", alpha = 0.5))) +
#   scale_color_manual(values = c("paired" = "#E69F00", "single" = "#56B4E9")) +
#   scale_fill_manual(values = c("paired" = "#E69F00", "single" = "#56B4E9")) +
#   theme_minimal()
# 
# # PC3 vs PC4
# ggplot(pca_data, aes(x = PC3, y = PC4, color = read_type)) +
#   geom_point(size = 3, alpha = 0.7) +
#   scale_color_manual(values = c("paired" = "#E69F00", "single" = "#56B4E9"),
#                      name = "Read Type") +
#   labs(
#     x = paste0("PC3 (", round(variance[3], 2), "%)"),
#     y = paste0("PC4 (", round(variance[4], 2), "%)"),
#     title = "PC3 vs PC4 by Read Type"
#   ) +
#   theme_minimal(base_size = 14)
# 
# # Summary statistics
# pca_data %>%
#   group_by(read_type) %>%
#   summarise(
#     n = n(),
#     mean_PC1 = mean(PC1),
#     mean_PC2 = mean(PC2),
#     sd_PC1 = sd(PC1),
#     sd_PC2 = sd(PC2)
#   )
# 
# # The issue: Bray-Curtis distance doesn't work with PCA scores (which can be negative)
# # Use Euclidean distance instead
# 
# # Correct PERMANOVA test
# adonis_result <- adonis2(pca_matrix ~ read_type, 
#                          data = pca_data, 
#                          method = "euclidean",  # Specify Euclidean distance
#                          permutations = 999)
# print(adonis_result)
# 
# # Alternative: Use betadisper to test for dispersion differences
# dist_matrix <- vegdist(pca_matrix, method = "euclidean")
# dispersion <- betadisper(dist_matrix, pca_data$read_type)
# anova(dispersion)
# permutest(dispersion, pairwise = TRUE)
# 
# # Visual check of dispersion
# plot(dispersion, hull = FALSE, ellipse = TRUE)
# 
# # Another approach: Simple MANOVA
# manova_result <- manova(cbind(PC1, PC2, PC3, PC4, PC5) ~ read_type, 
#                         data = pca_data)
# summary(manova_result)
# summary.aov(manova_result)  # Univariate tests for each PC
# 
# # T-tests for individual PCs
# t.test(PC1 ~ read_type, data = pca_data)
# t.test(PC2 ~ read_type, data = pca_data)
# 
# # Visualize the distributions
# library(ggplot2)
# 
# # Boxplots
# library(tidyr)
# pca_long <- pca_data %>%
#   select(read_type, PC1, PC2, PC3, PC4, PC5) %>%
#   pivot_longer(cols = starts_with("PC"), 
#                names_to = "PC", 
#                values_to = "Score")
# 
# ggplot(pca_long, aes(x = PC, y = Score, fill = read_type)) +
#   geom_boxplot(alpha = 0.7) +
#   scale_fill_manual(values = c("paired" = "#E69F00", "single" = "#56B4E9")) +
#   labs(title = "PC Scores Distribution by Read Type",
#        y = "PC Score",
#        fill = "Read Type") +
#   theme_minimal() +
#   theme(axis.text.x = element_text(angle = 0))
# 
# # Density plots for PC1 and PC2
# ggplot(pca_data, aes(x = PC1, fill = read_type)) +
#   geom_density(alpha = 0.5) +
#   scale_fill_manual(values = c("paired" = "#E69F00", "single" = "#56B4E9")) +
#   labs(title = "PC1 Distribution by Read Type") +
#   theme_minimal()
# 
# ggplot(pca_data, aes(x = PC2, fill = read_type)) +
#   geom_density(alpha = 0.5) +
#   scale_fill_manual(values = c("paired" = "#E69F00", "single" = "#56B4E9")) +
#   labs(title = "PC2 Distribution by Read Type") +
#   theme_minimal()
# 
# 
# 
# 
# 
