library(dplyr)
library(ggplot2)
library(data.table)

setwd("/Users/cperkins/Desktop/palanivelu_lab/diversity_panel/snp_gwas")

# ============================================================================
# Load GWAS results
# ============================================================================

# Load results for one phenotype (adjust number based on which you want)
gwas <- fread("gwas_pheno4_lmm.assoc.txt")

# Check the structure
head(gwas)
colnames(gwas)

# GEMMA output columns:
# chr, rs, ps, n_miss, allele1, allele0, af, beta, se, logl_H1, l_remle, p_wald, p_lrt, p_score

# ============================================================================
# Load phenotype info to know what each phenotype is
# ============================================================================

pheno_info <- read.csv("/Users/cperkins/Desktop/palanivelu_lab/diversity_panel/dp_vcfs/plink/phenotype_info.csv")
print(pheno_info)

# ============================================================================
# Manhattan Plot
# ============================================================================

# Prepare data
gwas_plot <- gwas %>%
  filter(!is.na(p_lrt)) %>%  # Use p_lrt from LMM
  mutate(
    chr_num = as.numeric(gsub("SL4.0ch", "", chr)),
    logP = -log10(p_lrt)
  ) %>%
  arrange(chr_num, ps)

# Calculate chromosome positions for x-axis
chr_lengths <- gwas_plot %>%
  group_by(chr_num) %>%
  summarise(max_pos = max(ps), .groups = "drop") %>%
  mutate(
    chr_start = cumsum(lag(max_pos, default = 0)),
    chr_mid = chr_start + max_pos/2
  )

gwas_plot <- gwas_plot %>%
  left_join(chr_lengths %>% select(chr_num, chr_start), by = "chr_num") %>%
  mutate(genome_pos = chr_start + ps)

# Significance thresholds based on your LD analysis
suggestive <- -log10(1e-5)  # Suggestive
significant <- -log10(1e-6)  # Genome-wide significant (based on your LD)
bonferroni <- -log10(0.05 / nrow(gwas_plot))  # Bonferroni correction

# Get phenotype name
pheno_num <- 4  # Change this based on which phenotype
pheno_name <- pheno_info$phenotype_name[pheno_num]

# Manhattan plot
p_manhattan <- ggplot(gwas_plot, aes(x = genome_pos, y = logP)) +
  geom_point(aes(color = as.factor(chr_num)), alpha = 0.8, size = 1) +
  geom_hline(yintercept = significant, color = "red", linetype = "dashed", size = 1) +
  geom_hline(yintercept = suggestive, color = "blue", linetype = "dashed", size = 0.8) +
  scale_color_manual(values = rep(c("#440154", "#2C728E"), 6), guide = "none") +
  scale_x_continuous(
    breaks = chr_lengths$chr_mid,
    labels = 1:12
  ) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, max(gwas_plot$logP) * 1.05)) +
  labs(
    title = paste0("GWAS Manhattan Plot - ", pheno_name),
    subtitle = paste0("n = 208 individuals | ", 
                      nrow(gwas_plot), " SNPs"),
    x = "Chromosome",
    y = expression(-log[10](p))
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank()
  )

print(p_manhattan)
ggsave(paste0("manhattan_", pheno_name, ".png"), p_manhattan, 
       width = 14, height = 6, dpi = 300)

# ============================================================================
# QQ Plot
# ============================================================================

observed <- sort(-log10(gwas$p_lrt[!is.na(gwas$p_lrt)]))
expected <- -log10(ppoints(length(observed)))

qq_data <- data.frame(
  expected = expected,
  observed = observed
)

# Calculate lambda (genomic inflation factor)
chisq <- qchisq(1 - gwas$p_lrt, 1)
lambda <- median(chisq, na.rm = TRUE) / qchisq(0.5, 1)

p_qq <- ggplot(qq_data, aes(x = expected, y = observed)) +
  geom_point(alpha = 0.5, size = 1.5) +
  geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") +
  annotate("text", x = max(expected) * 0.2, y = max(observed) * 0.9,
           label = paste0("λ = ", round(lambda, 3)),
           size = 5, fontface = "bold") +
  labs(
    title = paste0("QQ Plot - ", pheno_name),
    subtitle = "Expected vs. Observed -log10(p)",
    x = expression(Expected~~-log[10](p)),
    y = expression(Observed~~-log[10](p))
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5)
  )

print(p_qq)
ggsave(paste0("qq_", pheno_name, ".png"), p_qq, width = 8, height = 8, dpi = 300)

# ============================================================================
# Extract Top Hits
# ============================================================================

# Significant SNPs
sig_snps <- gwas %>%
  filter(p_lrt < 1e-6) %>%
  arrange(p_lrt) %>%
  select(chr, rs, ps, allele1, allele0, af, beta, se, p_lrt)

cat("\nNumber of genome-wide significant SNPs (p < 1e-6):", nrow(sig_snps), "\n")
if(nrow(sig_snps) > 0) {
  print(head(sig_snps, 20))
  write.csv(sig_snps, paste0("significant_snps_", pheno_name, ".csv"), 
            row.names = FALSE)
}

# Top 100 SNPs
top_snps <- gwas %>%
  arrange(p_lrt) %>%
  head(100) %>%
  select(chr, rs, ps, allele1, allele0, af, beta, se, p_lrt)

write.csv(top_snps, paste0("top100_snps_", pheno_name, ".csv"), row.names = FALSE)

# ============================================================================
# Summary by chromosome
# ============================================================================

chr_summary <- gwas %>%
  mutate(chr_num = as.numeric(gsub("SL4.0ch", "", chr))) %>%
  filter(p_lrt < 1e-5) %>%
  group_by(chr_num) %>%
  summarise(
    n_suggestive = n(),
    n_significant = sum(p_lrt < 1e-6),
    min_p = min(p_lrt),
    top_snp_pos = ps[which.min(p_lrt)],
    .groups = "drop"
  ) %>%
  arrange(desc(n_significant))

cat("\nSummary of hits by chromosome:\n")
print(chr_summary)

write.csv(chr_summary, paste0("chromosome_summary_", pheno_name, ".csv"), 
          row.names = FALSE)

# ============================================================================
# Lambda interpretation
# ============================================================================

cat("\n=== Genomic Inflation Factor (λ) ===\n")
cat("λ =", round(lambda, 3), "\n")
if(lambda < 0.95) {
  cat("Interpretation: Possible over-correction (λ < 0.95)\n")
} else if(lambda <= 1.05) {
  cat("Interpretation: Good! Well-controlled for population structure\n")
} else if(lambda <= 1.10) {
  cat("Interpretation: Slight inflation, acceptable\n")
} else {
  cat("Interpretation: High inflation (λ > 1.10), residual population structure\n")
}

cat("\n=== Analysis Complete! ===\n")
cat("Files saved:\n")
cat("  - manhattan_", pheno_name, ".png\n", sep = "")
cat("  - qq_", pheno_name, ".png\n", sep = "")
if(nrow(sig_snps) > 0) {
  cat("  - significant_snps_", pheno_name, ".csv\n", sep = "")
}
cat("  - top100_snps_", pheno_name, ".csv\n", sep = "")
cat("  - chromosome_summary_", pheno_name, ".csv\n", sep = "")
