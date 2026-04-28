#creating inputs for SNP GWAS

library(SNPRelate)

setwd("/Users/cperkins/Desktop/palanivelu_lab/diversity_panel/dp_vcfs/plink")

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



#input list for gwas



# ============================================================================
# 1. Create covariate file (PCs)
# ============================================================================

# GEMMA expects: FID IID PC1 PC2 PC3 ... (tab-separated, no header)
covariates <- pca_data %>%
  arrange(sample.id) %>%  # Match order to .fam file
  mutate(FID = sample.id,  # Create FID column
         IID = sample.id) %>%  # Create IID column (same as FID with --double-id)
  select(FID, IID, PC1, PC2, PC3, PC4, PC5)

# Check it looks right
head(covariates)
# Write covariate file
write.table(covariates, 
            "gemma_covariates.txt",
            quote = FALSE, 
            row.names = FALSE, 
            col.names = FALSE,
            sep = "\t")

cat("Created gemma_covariates.txt with", nrow(covariates), "samples\n")

# ============================================================================
# 2. Create phenotype file 
# ============================================================================


#pheno files

pheno_dir <- "/Users/cperkins/Desktop/palanivelu_lab/kmers_gwas/nextflow/gwas/input"

library(dplyr)
library(data.table)

setwd("/Users/cperkins/Desktop/palanivelu_lab/diversity_panel/dp_vcfs/plink")

# ============================================================================
# Load all phenotype files (excluding "short" versions)
# ============================================================================

pheno_dir <- "/Users/cperkins/Desktop/palanivelu_lab/kmers_gwas/nextflow/gwas/input"

# List all .pheno files, excluding ones with "short" in the name
pheno_files <- list.files(pheno_dir, pattern = "\\.pheno$", full.names = TRUE)
pheno_files <- pheno_files[!grepl("short", basename(pheno_files))]

cat("Found", length(pheno_files), "phenotype files (after excluding 'short'):\n")
print(basename(pheno_files))

# ============================================================================
# Load and merge all phenotype files
# ============================================================================

all_phenotypes <- list()

for (pheno_file in pheno_files) {
  # Extract phenotype name from filename
  pheno_name <- gsub("\\.pheno$", "", basename(pheno_file))
  
  # Read phenotype file - skip first row (it's a header)
  pheno_data <- read.table(pheno_file, header = FALSE, stringsAsFactors = FALSE, skip = 1)
  
  # Column 1 = sample ID (name_CW format like CW0045)
  # Column 2 = phenotype value
  colnames(pheno_data) <- c("sample_id", pheno_name)
  
  # Convert phenotype to numeric
  pheno_data[[pheno_name]] <- as.numeric(pheno_data[[pheno_name]])
  
  all_phenotypes[[pheno_name]] <- pheno_data
  
  cat("Loaded:", pheno_name, "-", nrow(pheno_data), "samples\n")
}

all_phenotypes[["locule"]] <- NULL

# Merge all phenotypes into one dataframe
pheno_merged <- all_phenotypes[[1]]

if (length(all_phenotypes) > 1) {
  for (i in 2:length(all_phenotypes)) {
    pheno_merged <- full_join(pheno_merged, all_phenotypes[[i]], by = "sample_id")
  }
}

cat("\nMerged phenotype data:\n")
cat("  Samples:", nrow(pheno_merged), "\n")
cat("  Phenotypes:", ncol(pheno_merged) - 1, "\n")
cat("\nPhenotype names:\n")
print(colnames(pheno_merged)[-1])

# Keep only samples that are in the PCA/genotype data
pheno_merged <- pheno_merged %>%
  filter(sample_id %in% pca_data$name_CW)

cat("\nAfter filtering to genotyped samples:\n")
cat("  Samples:", nrow(pheno_merged), "\n")

# ============================================================================
# Match phenotype sample IDs to your PCA sample IDs
# ============================================================================

cat("\n========================================\n")
cat("Checking sample ID overlap...\n")
cat("========================================\n")

# Phenotype IDs are in CW format (e.g., CW0045)
# PCA data has name_CW column with same format
cat("\nFirst few phenotype sample IDs:\n")
print(head(pheno_merged$sample_id))

cat("\nFirst few PCA name_CW IDs:\n")
print(head(pca_data$name_CW))

# Check for matches with name_CW
n_overlap <- sum(pheno_merged$sample_id %in% pca_data$name_CW)
cat("\nSamples in both phenotype and PCA data:", n_overlap, "\n")
cat("Samples only in phenotype data:", sum(!(pheno_merged$sample_id %in% pca_data$name_CW)), "\n")
cat("Samples only in PCA data:", sum(!(pca_data$name_CW %in% pheno_merged$sample_id)), "\n")

# ============================================================================
# Create GEMMA phenotype file
# ============================================================================

# Merge phenotypes with PCA data using name_CW
# This ensures we only keep samples that have BOTH genotype and phenotype data
gemma_pheno <- pca_data %>%
  select(sample.id, name_CW) %>%
  left_join(pheno_merged, by = c("name_CW" = "sample_id")) %>%
  arrange(sample.id) %>%
  mutate(FID = sample.id,
         IID = sample.id) %>%
  select(FID, IID, everything(), -sample.id, -name_CW)

# Check for missing data
cat("\nMissing data summary:\n")
missing_summary <- colSums(is.na(gemma_pheno))
print(missing_summary)

# Replace NA with -999 for GEMMA
gemma_pheno[is.na(gemma_pheno)] <- -999

# Preview
cat("\nFinal GEMMA phenotype file preview:\n")
print(head(gemma_pheno))
cat("\nDimensions:", nrow(gemma_pheno), "samples x", ncol(gemma_pheno)-2, "phenotypes\n")

# Write phenotype file
write.table(gemma_pheno, 
            "gemma_phenotypes.txt",
            quote = FALSE, 
            row.names = FALSE, 
            col.names = FALSE,
            sep = "\t")

cat("\n========================================\n")
cat("Files created successfully!\n")
cat("========================================\n")
cat("1. gemma_covariates.txt (", nrow(covariates), " samples x 5 PCs)\n", sep = "")
cat("2. gemma_phenotypes.txt (", nrow(gemma_pheno), " samples x ", 
    ncol(gemma_pheno)-2, " phenotypes)\n", sep = "")
cat("\nPhenotypes included:\n")
print(colnames(gemma_pheno)[-(1:2)])

# ============================================================================
# Create a phenotype info file for reference
# ============================================================================

pheno_info <- data.frame(
  phenotype_number = 1:(ncol(gemma_pheno)-2),
  phenotype_name = colnames(gemma_pheno)[-(1:2)],
  n_samples = colSums(gemma_pheno[,-(1:2)] != -999),
  mean_value = apply(gemma_pheno[,-(1:2)], 2, function(x) {
    vals <- x[x != -999]
    if(length(vals) > 0) mean(vals, na.rm = TRUE) else NA
  }),
  sd_value = apply(gemma_pheno[,-(1:2)], 2, function(x) {
    vals <- x[x != -999]
    if(length(vals) > 0) sd(vals, na.rm = TRUE) else NA
  }),
  min_value = apply(gemma_pheno[,-(1:2)], 2, function(x) {
    vals <- x[x != -999]
    if(length(vals) > 0) min(vals, na.rm = TRUE) else NA
  }),
  max_value = apply(gemma_pheno[,-(1:2)], 2, function(x) {
    vals <- x[x != -999]
    if(length(vals) > 0) max(vals, na.rm = TRUE) else NA
  })
)

write.csv(pheno_info, "phenotype_info.csv", row.names = FALSE)

cat("\nPhenotype summary:\n")
print(pheno_info)

# ============================================================================
# Verify files are ready for upload
# ============================================================================

cat("\n========================================\n")
cat("Ready to upload to HPC!\n")
cat("========================================\n")
cat("\nFiles created in:", getwd(), "\n")
cat("  - gemma_covariates.txt\n")
cat("  - gemma_phenotypes.txt\n")
cat("  - phenotype_info.csv (for reference)\n")
cat("\nUpload commands:\n")
cat("cd", getwd(), "\n")
cat("scp gemma_covariates.txt cjperkins1@hpc.arizona.edu:/xdisk/yadegari/cjperkins1/GWAS/covariates/\n")
cat("scp gemma_phenotypes.txt cjperkins1@hpc.arizona.edu:/xdisk/yadegari/cjperkins1/GWAS/covariates/\n")
cat("\nThen on HPC, run:\n")
cat("cd /xdisk/yadegari/cjperkins1/GWAS/slurm\n")
cat("sbatch gemma_gwas.sh\n")

library(dplyr)

setwd("/Users/cperkins/Desktop/palanivelu_lab/diversity_panel/dp_vcfs/plink")

# Load your PCA data
# ... (your existing code to load pca_data)

# REMOVE THE REFERENCE BEFORE CREATING FILES
pca_data_no_ref <- pca_data %>%
  filter(name_CW != "CW0000")

# Now create covariates using the filtered data
covariates <- pca_data_no_ref %>%
  arrange(sample.id) %>%
  mutate(FID = sample.id, IID = sample.id) %>%
  select(FID, IID, PC1, PC2, PC3, PC4, PC5)

write.table(covariates, 
            "gemma_covariates.txt",
            quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")

# Create phenotypes using the filtered data
gemma_pheno <- pca_data_no_ref %>%
  select(sample.id, name_CW) %>%
  left_join(pheno_merged, by = c("name_CW" = "sample_id")) %>%
  arrange(sample.id) %>%
  mutate(FID = sample.id, IID = sample.id) %>%
  select(FID, IID, everything(), -sample.id, -name_CW)

# Replace NA with -999
gemma_pheno[is.na(gemma_pheno)] <- -999

write.table(gemma_pheno, 
            "gemma_phenotypes.txt",
            quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")


# Find the reference sample's PCs
ref_pcs <- pca_data %>%
  filter(grepl("SRX10558138", sample.id)) %>%
  select(sample.id, PC1, PC2, PC3, PC4, PC5)

print(ref_pcs)

# Create the covariate line
ref_line <- ref_pcs %>%
  mutate(FID = sample.id, IID = sample.id) %>%
  select(FID, IID, PC1, PC2, PC3, PC4, PC5)

# Show the line to add
cat(paste(ref_line$FID, ref_line$IID, ref_line$PC1, ref_line$PC2, ref_line$PC3, ref_line$PC4, ref_line$PC5, sep = "\t"), "\n")

# You'll need to create this based on your actual phenotype data
# GEMMA expects: FID IID pheno1 pheno2 ... (tab-separated, no header)
# Use -999 or NA for missing values

# Example if you have phenotype data:
# phenotypes <- pca_data %>%
#   left_join(your_phenotype_data, by = "sample.id") %>%
#   arrange(sample.id) %>%
#   select(sample.id, sample.id, trait1, trait2, trait3) %>%
#   rename(FID = sample.id...1, IID = sample.id...2)
# 
# write.table(phenotypes,
#             "gemma_phenotypes.txt",
#             quote = FALSE,
#             row.names = FALSE,
#             col.names = FALSE,
#             sep = "\t",
#             na = "-999")

# For now, create a dummy phenotype file as template:
# (Replace this with your actual phenotype data!)
dummy_pheno <- covariates %>%
  select(FID, IID) %>%
  mutate(phenotype1 = rnorm(n()),  # Replace with real phenotype!
         phenotype2 = rnorm(n()))   # Can have multiple phenotypes

write.table(dummy_pheno,
            "gemma_phenotypes_TEMPLATE.txt",
            quote = FALSE,
            row.names = FALSE,
            col.names = FALSE,
            sep = "\t")

cat("\nIMPORTANT: Replace gemma_phenotypes_TEMPLATE.txt with your actual phenotype data!\n")
cat("Format: FID IID pheno1 pheno2 ... (tab-separated, no header, -999 for missing)\n")