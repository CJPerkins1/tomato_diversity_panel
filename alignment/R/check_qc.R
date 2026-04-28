#Analyzing Bwamem2 aligned read metrics


library(dplyr)
library(stringr)
library(ggplot2)

setwd("~/Desktop/palanivelu_lab/diversity_panel/dp_alignment")


# ============================================
# 1. PREPARE DATA
# ============================================

short_reads <- read.delim("bam_qc/Illumina_all_samples_bam_qc.tsv", sep = '\t')

long_reads <- read.delim("bam_qc/ONT_all_samples_bam_qc.tsv", sep = '\t')

metadata <- read.csv("/Users/cperkins/Desktop/palanivelu_lab/diversity_panel/data/supplementary_table.csv")

# Method 1: Using str_extract to pull out SRR ID
long_reads <- long_reads %>%
  mutate(Run = str_extract(sample, "SRR[0-9]+"))

# Method 2: Alternative using gsub if str_extract doesn't work
# long_reads <- long_reads %>%
#   mutate(Run = gsub(".*_(SRR[0-9]+)_.*", "\\1", sample))

# Now join with your metadata
long_reads_annotated <- long_reads %>%
  left_join(metadata, by = "Run")

# Check the join worked
head(long_reads_annotated)
nrow(long_reads)  # should equal
nrow(long_reads_annotated)  # should equal

# Check for any samples that didn't match
long_reads_annotated %>%
  filter(is.na(Accession)) %>%
  select(sample, Run)


# Extract Run IDs from short reads
short_reads <- short_reads %>%
  mutate(Run = str_extract(sample, "SRR[0-9]+"))

# Join with metadata
short_reads_annotated <- short_reads %>%
  left_join(metadata, by = "Run") %>%
  filter(!is.na(coverage_mean),
         !is.na(mapped_pct),
         !is.na(mean_mapq),
         !is.na(mq0_pct),
         !is.na(pct_zero_coverage))

# Add dataset labels
long_reads_annotated$tech <- "ONT"
short_reads_annotated$tech <- "Illumina"

# Combine for comparison
combined_qc <- bind_rows(
  long_reads_annotated %>% 
    select(sample, Run, tech, coverage_mean, mapped_pct, 
           mean_mapq, median_mapq, mq0_pct, pct_zero_coverage),
  short_reads_annotated %>% 
    select(sample, Run, tech, coverage_mean, mapped_pct, 
           mean_mapq, median_mapq, mq0_pct, pct_zero_coverage)
)

# ============================================
# 2. DETERMINE --max-depth VALUE
# ============================================

cat("\n========== MAX DEPTH RECOMMENDATION ==========\n")

# Calculate per-tech coverage stats
depth_stats <- combined_qc %>%
  group_by(tech) %>%
  summarise(
    n = n(),
    mean_cov = mean(coverage_mean),
    median_cov = median(coverage_mean),
    sd_cov = sd(coverage_mean),
    max_cov = max(coverage_mean),
    q95_cov = quantile(coverage_mean, 0.95),
    .groups = 'drop'
  )

print(depth_stats)

# Calculate recommended max-depth (2.5x the 95th percentile)
max_depth_ont <- ceiling(depth_stats %>% filter(tech == "ONT") %>% pull(q95_cov) * 2.5)
max_depth_illumina <- ceiling(depth_stats %>% filter(tech == "Illumina") %>% pull(q95_cov) * 2.5)

cat("\nRecommended --max-depth values:\n")
cat(sprintf("  ONT: %d (2.5x the 95th percentile coverage)\n", max_depth_ont))
cat(sprintf("  Illumina: %d (2.5x the 95th percentile coverage)\n", max_depth_illumina))
cat(sprintf("  Conservative for both: %d\n", max(max_depth_ont, max_depth_illumina)))

# Visualize coverage distributions
ggplot(combined_qc, aes(x = coverage_mean, fill = tech)) +
  geom_histogram(bins = 30, alpha = 0.6, position = "identity") +
  geom_vline(data = depth_stats, aes(xintercept = q95_cov, color = tech), 
             linetype = "dashed", size = 1) +
  labs(title = "Coverage Distribution by Technology",
       subtitle = "Dashed lines = 95th percentile",
       x = "Mean Coverage",
       y = "Count") +
  theme_minimal()

# ============================================
# 3. DETERMINE --min-MQ VALUE
# ============================================

cat("\n========== MIN MAPQ RECOMMENDATION ==========\n")

mapq_stats <- combined_qc %>%
  group_by(tech) %>%
  summarise(
    n = n(),
    mean_mapq = mean(mean_mapq),
    median_mapq = median(median_mapq),
    mean_mq0_pct = mean(mq0_pct),
    median_mq0_pct = median(mq0_pct),
    .groups = 'drop'
  )

print(mapq_stats)

# Calculate what percentage of reads would be kept at different thresholds
# Approximate based on MQ distribution (assuming reads spread across MQ bins)
ont_mq0 <- combined_qc %>% filter(tech == "ONT") %>% pull(mq0_pct) %>% mean()
illumina_mq0 <- combined_qc %>% filter(tech == "Illumina") %>% pull(mq0_pct) %>% mean()

cat("\nMean % of mapped reads with MQ=0:\n")
cat(sprintf("  ONT: %.1f%%\n", ont_mq0))
cat(sprintf("  Illumina: %.1f%%\n", illumina_mq0))

cat("\nRecommended --min-MQ values:\n")
if (ont_mq0 > 30) {
  cat("  ONT: 10-15 (high MQ0 rate, be lenient)\n")
} else if (ont_mq0 > 20) {
  cat("  ONT: 20 (moderate MQ0 rate, standard threshold)\n")
} else {
  cat("  ONT: 25 (low MQ0 rate, can be stricter)\n")
}

if (illumina_mq0 > 20) {
  cat("  Illumina: 20 (moderate MQ0 rate)\n")
} else {
  cat("  Illumina: 25 (low MQ0 rate, use standard threshold)\n")
}

# Visualize MAPQ distributions
ggplot(combined_qc, aes(x = tech, y = mean_mapq, fill = tech)) +
  geom_boxplot() +
  geom_hline(yintercept = c(20, 25, 30), linetype = "dashed", alpha = 0.5) +
  labs(title = "Mean MAPQ Distribution",
       subtitle = "Dashed lines at MQ 20, 25, 30",
       y = "Mean MAPQ") +
  theme_minimal()

# ============================================
# 4. SAMPLE-LEVEL QC FOR BAM LIST
# ============================================
# ============================================
# UPDATED QC THRESHOLDS (Based on Your Data)
# ============================================

qc_thresholds <- list(
  # Sample-level filters (BEFORE variant calling)
  ont = list(
    mapping_min = 70,      # ONT mapping varies widely, 70% is reasonable floor
    coverage_min = 10,     # Mean is 44x, so 10x minimum is conservative
    zero_cov_max = 10,     # Allow up to 10% gaps
    mq0_max = 40           # Mean is 30%, cap at 40% (exclude worst ~15%)
  ),
  
  illumina = list(
    mapping_min = 70,      # Illumina generally maps well
    coverage_min = 8,      # Mean is 17x, so 8x is reasonable minimum
    zero_cov_max = 5,      # Stricter - expect better genome coverage
    mq0_max = 25           # Mean is 17%, cap at 25% (exclude worst ~10%)
  ),
  
  # Variant calling parameters
  variant_calling = list(
    ont_max_depth = 175,   # 2.5x the 95th percentile (69x)
    ont_min_mq = 20,       # Standard threshold despite 30% MQ0
    ont_min_bq = 20,       # Standard base quality
    
    illumina_max_depth = 100,  # 2.5x the 95th percentile (33x)
    illumina_min_mq = 25,      # Can be stricter (only 17% MQ0)
    illumina_min_bq = 20       # Standard base quality
  )
)

# Apply QC filters
combined_qc <- combined_qc %>%
  mutate(
    # ONT-specific flags
    flag_low_mapping = ifelse(tech == "ONT", 
                              mapped_pct < qc_thresholds$ont$mapping_min,
                              mapped_pct < qc_thresholds$illumina$mapping_min),
    
    flag_low_coverage = ifelse(tech == "ONT",
                               coverage_mean < qc_thresholds$ont$coverage_min,
                               coverage_mean < qc_thresholds$illumina$coverage_min),
    
    flag_high_zero_cov = ifelse(tech == "ONT",
                                pct_zero_coverage > qc_thresholds$ont$zero_cov_max,
                                pct_zero_coverage > qc_thresholds$illumina$zero_cov_max),
    
    flag_high_mq0 = ifelse(tech == "ONT",
                           mq0_pct > qc_thresholds$ont$mq0_max,
                           mq0_pct > qc_thresholds$illumina$mq0_max),
    
    total_flags = flag_low_mapping + flag_low_coverage + 
      flag_high_zero_cov + flag_high_mq0,
    
    qc_status = case_when(
      total_flags == 0 ~ "PASS",
      total_flags == 1 ~ "WARNING", 
      total_flags >= 2 ~ "FAIL"
    )
  )

# Summary
cat("\n========== SAMPLE QC RESULTS ==========\n")
qc_summary <- combined_qc %>%
  group_by(tech, qc_status) %>%
  summarise(n = n(), .groups = 'drop') %>%
  tidyr::pivot_wider(names_from = qc_status, values_from = n, values_fill = 0)

print(qc_summary)

# Samples passing QC (for variant calling)
ont_pass <- combined_qc %>%
  filter(tech == "ONT", qc_status %in% c("PASS", "WARNING")) %>%
  arrange(sample)

illumina_pass <- combined_qc %>%
  filter(tech == "Illumina", qc_status %in% c("PASS", "WARNING")) %>%
  arrange(sample)

cat(sprintf("\nSamples for variant calling:\n"))
cat(sprintf("  ONT: %d samples (including %d with warnings)\n", 
            nrow(ont_pass),
            sum(combined_qc$tech == "ONT" & combined_qc$qc_status == "WARNING")))
cat(sprintf("  Illumina: %d samples (including %d with warnings)\n", 
            nrow(illumina_pass),
            sum(combined_qc$tech == "Illumina" & combined_qc$qc_status == "WARNING")))

# Save passing sample lists
write.table(ont_pass$sample, 
            "ont_samples_pass_qc.txt", 
            row.names = FALSE, col.names = FALSE, quote = FALSE)

write.table(illumina_pass$sample, 
            "illumina_samples_pass_qc.txt", 
            row.names = FALSE, col.names = FALSE, quote = FALSE)

cat("\n========== FINAL PARAMETERS ==========\n")
cat("\nONT variant calling:\n")
cat(sprintf("  --max-depth %d\n", qc_thresholds$variant_calling$ont_max_depth))
cat(sprintf("  --min-MQ %d\n", qc_thresholds$variant_calling$ont_min_mq))
cat(sprintf("  --min-BQ %d\n", qc_thresholds$variant_calling$ont_min_bq))
cat(sprintf("  Samples: %d\n", nrow(ont_pass)))

cat("\nIllumina variant calling:\n")
cat(sprintf("  --max-depth %d\n", qc_thresholds$variant_calling$illumina_max_depth))
cat(sprintf("  --min-MQ %d\n", qc_thresholds$variant_calling$illumina_min_mq))
cat(sprintf("  --min-BQ %d\n", qc_thresholds$variant_calling$illumina_min_bq))
cat(sprintf("  Samples: %d\n", nrow(illumina_pass)))

# Check which Illumina samples failed
illumina_failed <- combined_qc %>%
  filter(tech == "Illumina", qc_status == "FAIL") %>%
  select(sample, Run, mapped_pct, coverage_mean, pct_zero_coverage, mq0_pct,
         flag_low_mapping, flag_low_coverage, flag_high_zero_cov, flag_high_mq0)

cat("Illumina samples that failed QC:\n")
print(illumina_failed)

# Summary of why they failed
cat("\nFailure reasons:\n")
cat(sprintf("  Low mapping (<%d%%): %d samples\n", 
            qc_thresholds$illumina$mapping_min,
            sum(illumina_failed$flag_low_mapping)))
cat(sprintf("  Low coverage (<%.0fx): %d samples\n", 
            qc_thresholds$illumina$coverage_min,
            sum(illumina_failed$flag_low_coverage)))
cat(sprintf("  High zero coverage (>%d%%): %d samples\n", 
            qc_thresholds$illumina$zero_cov_max,
            sum(illumina_failed$flag_high_zero_cov)))
cat(sprintf("  High MQ0 (>%d%%): %d samples\n", 
            qc_thresholds$illumina$mq0_max,
            sum(illumina_failed$flag_high_mq0)))

# Also check WARNING samples
illumina_warning <- combined_qc %>%
  filter(tech == "Illumina", qc_status == "WARNING") %>%
  select(sample, Run, mapped_pct, coverage_mean, pct_zero_coverage, mq0_pct,
         flag_low_mapping, flag_low_coverage, flag_high_zero_cov, flag_high_mq0)

cat(sprintf("\nIllumina samples with warnings: %d\n", nrow(illumina_warning)))
if (nrow(illumina_warning) > 0) {
  print(illumina_warning)
}

# Recreate combined_qc with population info
combined_qc <- bind_rows(
  long_reads_annotated %>% 
    select(sample, Run, tech, Accession, assigned_population_filled1,
           coverage_mean, mapped_pct, mean_mapq, median_mapq, 
           mq0_pct, pct_zero_coverage),
  short_reads_annotated %>% 
    select(sample, Run, tech, Accession, assigned_population_filled1,
           coverage_mean, mapped_pct, mean_mapq, median_mapq, 
           mq0_pct, pct_zero_coverage)
) %>%
  filter(complete.cases(coverage_mean, mapped_pct, mean_mapq, 
                        mq0_pct, pct_zero_coverage))

# Re-apply QC flags
combined_qc <- combined_qc %>%
  mutate(
    flag_low_mapping = ifelse(tech == "ONT", 
                              mapped_pct < 70,
                              mapped_pct < 70),
    flag_low_coverage = ifelse(tech == "ONT",
                               coverage_mean < 10,
                               coverage_mean < 8),
    flag_high_zero_cov = ifelse(tech == "ONT",
                                pct_zero_coverage > 10,
                                pct_zero_coverage > 5),
    flag_high_mq0 = ifelse(tech == "ONT",
                           mq0_pct > 40,
                           mq0_pct > 25),
    total_flags = flag_low_mapping + flag_low_coverage + 
      flag_high_zero_cov + flag_high_mq0,
    qc_status = case_when(
      total_flags == 0 ~ "PASS",
      total_flags == 1 ~ "WARNING", 
      total_flags >= 2 ~ "FAIL"
    )
  )

# Now check populations
cat("========== FAILED ILLUMINA SAMPLES BY POPULATION ==========\n")

illumina_failed_pops <- combined_qc %>%
  filter(tech == "Illumina", qc_status == "FAIL") %>%
  select(sample, Run, Accession, assigned_population_filled1, 
         coverage_mean, pct_zero_coverage, mq0_pct)

print(illumina_failed_pops)

# Summary by population
cat("\nFailed samples by population:\n")
failed_pop_summary <- illumina_failed_pops %>%
  group_by(assigned_population_filled1) %>%
  summarise(n = n(), .groups = 'drop') %>%
  arrange(desc(n))
print(failed_pop_summary)

# Overall distribution
cat("\n========== OVERALL ILLUMINA POPULATION DISTRIBUTION ==========\n")
overall_pop_summary <- combined_qc %>%
  filter(tech == "Illumina") %>%
  group_by(assigned_population_filled1, qc_status) %>%
  summarise(n = n(), .groups = 'drop') %>%
  tidyr::pivot_wider(names_from = qc_status, values_from = n, values_fill = 0) %>%
  mutate(total = PASS + WARNING + FAIL,
         fail_pct = round(100 * FAIL / total, 1))

print(overall_pop_summary)

# Compare coverage and mapping by population
cat("========== COVERAGE & MAPPING BY POPULATION ==========\n")

pop_comparison <- combined_qc %>%
  filter(tech == "Illumina") %>%
  group_by(assigned_population_filled1) %>%
  summarise(
    n = n(),
    mean_coverage = mean(coverage_mean),
    mean_mapping = mean(mapped_pct),
    mean_zero_cov = mean(pct_zero_coverage),
    mean_mq0 = mean(mq0_pct),
    .groups = 'drop'
  ) %>%
  arrange(desc(mean_mq0))

print(pop_comparison)

# Specifically compare SLL vs SP
cat("\n========== SLL vs SP COMPARISON ==========\n")

sll_vs_sp <- combined_qc %>%
  filter(tech == "Illumina") %>%
  mutate(group = case_when(
    grepl("^SLL", assigned_population_filled1) ~ "SLL (cultivated)",
    grepl("^SP", assigned_population_filled1) ~ "SP (wild)",
    grepl("^SLC", assigned_population_filled1) ~ "SLC (intermediate)",
    TRUE ~ "Other"
  )) %>%
  group_by(group) %>%
  summarise(
    n = n(),
    coverage = mean(coverage_mean),
    mapping = mean(mapped_pct),
    zero_cov = mean(pct_zero_coverage),
    mq0 = mean(mq0_pct),
    .groups = 'drop'
  )

print(sll_vs_sp)

# Apply species-aware QC
combined_qc <- combined_qc %>%
  mutate(
    species_group = case_when(
      grepl("^SP", assigned_population_filled1) ~ "wild",
      grepl("^SLC", assigned_population_filled1) ~ "intermediate",
      grepl("^SLL", assigned_population_filled1) ~ "cultivated",
      TRUE ~ "other"
    ),
    
    # Adjusted thresholds accounting for biological divergence
    flag_high_zero_cov_adj = case_when(
      species_group == "wild" & tech == "Illumina" ~ pct_zero_coverage > 8,
      species_group == "intermediate" & tech == "Illumina" ~ pct_zero_coverage > 6,
      tech == "Illumina" ~ pct_zero_coverage > 5,
      TRUE ~ pct_zero_coverage > 10  # ONT threshold
    ),
    
    flag_high_mq0_adj = case_when(
      species_group == "wild" & tech == "Illumina" ~ mq0_pct > 30,
      tech == "Illumina" ~ mq0_pct > 25,
      TRUE ~ mq0_pct > 40  # ONT threshold
    ),
    
    # Recalculate with adjusted flags
    total_flags_adj = flag_low_mapping + flag_low_coverage + 
      flag_high_zero_cov_adj + flag_high_mq0_adj,
    
    qc_status_adj = case_when(
      total_flags_adj == 0 ~ "PASS",
      total_flags_adj == 1 ~ "WARNING",
      total_flags_adj >= 2 ~ "FAIL"
    )
  )

# Compare original vs adjusted QC
cat("========== QC STATUS COMPARISON ==========\n")
qc_comparison <- combined_qc %>%
  filter(tech == "Illumina") %>%
  group_by(species_group, qc_status, qc_status_adj) %>%
  summarise(n = n(), .groups = 'drop') %>%
  arrange(species_group, qc_status)

print(qc_comparison)

# Final counts by species
cat("\n========== FINAL SAMPLE COUNTS (ADJUSTED) ==========\n")
final_counts <- combined_qc %>%
  filter(tech == "Illumina", qc_status_adj %in% c("PASS", "WARNING")) %>%
  group_by(species_group) %>%
  summarise(n = n(), .groups = 'drop')

print(final_counts)

cat(sprintf("\nTotal Illumina samples for calling: %d\n", 
            sum(combined_qc$tech == "Illumina" & 
                  combined_qc$qc_status_adj %in% c("PASS", "WARNING"))))

# Save final sample lists
ont_for_calling <- combined_qc %>%
  filter(tech == "ONT", qc_status_adj %in% c("PASS", "WARNING"))

illumina_for_calling <- combined_qc %>%
  filter(tech == "Illumina", qc_status_adj %in% c("PASS", "WARNING"))

write.table(ont_for_calling$sample, 
            "ont_samples_for_calling.txt", 
            row.names = FALSE, col.names = FALSE, quote = FALSE)

write.table(illumina_for_calling$sample, 
            "illumina_samples_for_calling.txt", 
            row.names = FALSE, col.names = FALSE, quote = FALSE)

cat("\nSample lists saved:\n")
cat("  ont_samples_for_calling.txt\n")
cat("  illumina_samples_for_calling.txt\n")

# Identify samples to EXCLUDE (FAIL status only, with adjusted thresholds)
ont_exclude <- combined_qc %>%
  filter(tech == "ONT", qc_status_adj == "FAIL") %>%
  pull(sample)

illumina_exclude <- combined_qc %>%
  filter(tech == "Illumina", qc_status_adj == "FAIL") %>%
  pull(sample)

cat("========== SAMPLES TO EXCLUDE ==========\n")
cat(sprintf("ONT samples to exclude: %d\n", length(ont_exclude)))
cat(sprintf("Illumina samples to exclude: %d\n", length(illumina_exclude)))

# Save exclude lists
write.table(ont_exclude, 
            "ont_samples_exclude.txt", 
            row.names = FALSE, col.names = FALSE, quote = FALSE)

write.table(illumina_exclude, 
            "illumina_samples_exclude.txt", 
            row.names = FALSE, col.names = FALSE, quote = FALSE)

# Show which samples are being excluded
if (length(ont_exclude) > 0) {
  cat("\nONT samples excluded:\n")
  ont_exclude_details <- combined_qc %>%
    filter(tech == "ONT", qc_status_adj == "FAIL") %>%
    select(sample, coverage_mean, mapped_pct, pct_zero_coverage, mq0_pct)
  print(ont_exclude_details)
}

if (length(illumina_exclude) > 0) {
  cat("\nIllumina samples excluded:\n")
  illumina_exclude_details <- combined_qc %>%
    filter(tech == "Illumina", qc_status_adj == "FAIL") %>%
    select(sample, assigned_population_filled1, coverage_mean, 
           mapped_pct, pct_zero_coverage, mq0_pct)
  print(illumina_exclude_details)
}

cat("\n========== FINAL COUNTS ==========\n")
cat(sprintf("ONT samples for calling: %d (total %d - excluded %d)\n", 
            sum(combined_qc$tech == "ONT") - length(ont_exclude),
            sum(combined_qc$tech == "ONT"),
            length(ont_exclude)))
cat(sprintf("Illumina samples for calling: %d (total %d - excluded %d)\n", 
            sum(combined_qc$tech == "Illumina") - length(illumina_exclude),
            sum(combined_qc$tech == "Illumina"),
            length(illumina_exclude)))

cat("\nExclude lists saved:\n")
cat("  ont_samples_exclude.txt\n")
cat("  illumina_samples_exclude.txt\n")

# ONT mapping rate statistics
ont_mapping_stats <- combined_qc %>%
  filter(tech == "ONT") %>%
  summarise(
    n = n(),
    mean_mapping = mean(mapped_pct),
    median_mapping = median(mapped_pct),
    sd_mapping = sd(mapped_pct),
    min_mapping = min(mapped_pct),
    max_mapping = max(mapped_pct),
    q25 = quantile(mapped_pct, 0.25),
    q75 = quantile(mapped_pct, 0.75)
  )

print(ont_mapping_stats)

# Also show distribution
cat("\nONT Mapping Rate Distribution:\n")
ont_mapping_hist <- combined_qc %>%
  filter(tech == "ONT") %>%
  mutate(mapping_bin = cut(mapped_pct, 
                           breaks = c(0, 70, 80, 90, 95, 100),
                           labels = c("<70%", "70-80%", "80-90%", "90-95%", "95-100%"))) %>%
  group_by(mapping_bin) %>%
  summarise(n = n(), .groups = 'drop')

print(ont_mapping_hist)

# Compare total reads between technologies
cat("========== TOTAL READS COMPARISON ==========\n")

# ONT total reads stats
ont_reads_stats <- long_reads_annotated %>%
  summarise(
    n_samples = n(),
    mean_total_reads = mean(total_reads),
    median_total_reads = median(total_reads),
    sd_total_reads = sd(total_reads),
    min_total_reads = min(total_reads),
    max_total_reads = max(total_reads)
  )

cat("\nONT Total Reads:\n")
print(ont_reads_stats)

# Illumina total reads stats
illumina_reads_stats <- short_reads_annotated %>%
  summarise(
    n_samples = n(),
    mean_total_reads = mean(total_reads),
    median_total_reads = median(total_reads),
    sd_total_reads = sd(total_reads),
    min_total_reads = min(total_reads),
    max_total_reads = max(total_reads)
  )

cat("\nIllumina Total Reads:\n")
print(illumina_reads_stats)

# Compare side by side
comparison <- data.frame(
  metric = c("Mean total reads", "Median total reads", "Mean coverage", 
             "Reads per 1x coverage"),
  ONT = c(
    mean(long_reads_annotated$total_reads),
    median(long_reads_annotated$total_reads),
    mean(long_reads_annotated$coverage_mean),
    mean(long_reads_annotated$total_reads) / mean(long_reads_annotated$coverage_mean)
  ),
  Illumina = c(
    mean(short_reads_annotated$total_reads),
    median(short_reads_annotated$total_reads),
    mean(short_reads_annotated$coverage_mean),
    mean(short_reads_annotated$total_reads) / mean(short_reads_annotated$coverage_mean)
  )
)

cat("\n========== SIDE-BY-SIDE COMPARISON ==========\n")
print(comparison)

cat("\nRatio (ONT / Illumina):\n")
cat(sprintf("  Total reads ratio: %.2fx\n", 
            mean(long_reads_annotated$total_reads) / mean(short_reads_annotated$total_reads)))
cat(sprintf("  Coverage ratio: %.2fx\n",
            mean(long_reads_annotated$coverage_mean) / mean(short_reads_annotated$coverage_mean)))

# Also get mapping rate for ONT while we're at it
ont_mapping_stats <- combined_qc %>%
  filter(tech == "ONT") %>%
  summarise(
    mean_mapping = mean(mapped_pct),
    median_mapping = median(mapped_pct),
    sd_mapping = sd(mapped_pct),
    min_mapping = min(mapped_pct),
    max_mapping = max(mapped_pct)
  )

cat("\n========== ONT MAPPING RATE (EXACT) ==========\n")
print(ont_mapping_stats)
