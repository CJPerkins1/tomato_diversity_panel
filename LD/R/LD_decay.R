# LD Decay Analysis for Multiple Chromosomes
# Processes PLINK --r2 output for chromosomes 1-12
# Provides combined analysis for GWAS linkage block estimation

library(ggplot2)
library(dplyr)
library(tidyr)

# ===== CORE FUNCTIONS =====

# Function to read and process PLINK LD output
read_ld_data <- function(file_path) {
  cat("Reading file:", file_path, "\n")
  
  # Handle both gzipped and regular files
  if (grepl("\\.gz$", file_path)) {
    ld_data <- read.table(gzfile(file_path), header = TRUE, stringsAsFactors = FALSE)
  } else {
    ld_data <- read.table(file_path, header = TRUE, stringsAsFactors = FALSE)
  }
  
  # Calculate distance between SNPs
  ld_data$DIST <- abs(ld_data$BP_B - ld_data$BP_A)
  
  cat("Total SNP pairs:", nrow(ld_data), "\n")
  
  return(ld_data)
}

# Function to bin LD values by distance
bin_ld_decay <- function(ld_data, bin_size = 5000, max_dist = 1000000) {
  ld_data %>%
    filter(DIST <= max_dist, DIST > 0) %>%
    mutate(dist_bin = cut(DIST, 
                          breaks = seq(0, max_dist, by = bin_size),
                          labels = seq(bin_size/2, max_dist - bin_size/2, by = bin_size))) %>%
    group_by(dist_bin) %>%
    summarize(
      mean_r2 = mean(R2, na.rm = TRUE),
      median_r2 = median(R2, na.rm = TRUE),
      n = n(),
      .groups = 'drop'
    ) %>%
    mutate(distance = as.numeric(as.character(dist_bin)))
}

# Function to plot LD decay
plot_ld_decay <- function(binned_data, output_file, title = "LD Decay") {
  cat("Creating plot:", output_file, "\n")
  
  p <- ggplot(binned_data, aes(x = distance/1000, y = mean_r2)) +
    geom_line(color = "blue", linewidth = 1) +
    geom_point(size = 2, alpha = 0.6) +
    labs(
      title = title,
      x = "Distance (kb)",
      y = expression(paste("Mean ", r^2))
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
      axis.title = element_text(size = 14),
      axis.text = element_text(size = 12)
    ) +
    geom_hline(yintercept = 0.2, linetype = "dashed", color = "red") +
    annotate("text", x = max(binned_data$distance/1000) * 0.8, y = 0.25, 
             label = "r² = 0.2 threshold", color = "red")
  
  ggsave(output_file, p, width = 10, height = 6, dpi = 300)
  
  return(p)
}

# Function to find decay distance
find_decay_distance <- function(binned_data, threshold = 0.2) {
  decay_point <- binned_data %>% 
    filter(mean_r2 < threshold) %>% 
    slice(1) %>% 
    pull(distance)
  
  if(length(decay_point) > 0) {
    return(decay_point)
  } else {
    return(NA)
  }
}

# ===== MAIN ANALYSIS FUNCTION =====

analyze_chromosome_ld <- function(chr, data_dir, output_dir, 
                                  bin_size = 10000, 
                                  max_dist = 5000000,
                                  file_pattern = "ld_SL4.0ch%02d.ld.gz") {
  
  cat("\n========================================\n")
  cat("ANALYZING CHROMOSOME", chr, "\n")
  cat("========================================\n")
  
  # Construct file path
  file_name <- sprintf(file_pattern, chr)
  file_path <- file.path(data_dir, file_name)
  
  # Check if file exists
  if (!file.exists(file_path)) {
    cat("WARNING: File not found:", file_path, "\n")
    return(NULL)
  }
  
  # Read data
  ld_data <- read_ld_data(file_path)
  
  # Bin the data
  cat("Binning LD values (bin size:", bin_size, "bp, max distance:", max_dist/1e6, "Mb)\n")
  binned_data <- bin_ld_decay(ld_data, bin_size = bin_size, max_dist = max_dist)
  
  # Create output file name
  output_file <- file.path(output_dir, sprintf("ld_decay_chr%02d.png", chr))
  
  # Plot
  plot_ld_decay(binned_data, output_file, 
                title = sprintf("LD Decay - Chromosome %d (0-%.0f Mb)", chr, max_dist/1e6))
  
  # Find decay distance
  decay_dist <- find_decay_distance(binned_data, threshold = 0.2)
  
  if(!is.na(decay_dist)) {
    cat("LD decays to r² < 0.2 at:", round(decay_dist/1000, 1), "kb\n")
  } else {
    cat("LD does NOT decay below r² = 0.2 within", max_dist/1e6, "Mb\n")
  }
  
  # Summary statistics
  cat("\n--- Distance Statistics ---\n")
  cat("Mean r² at 100-500kb:", round(mean(ld_data$R2[ld_data$DIST >= 100000 & ld_data$DIST < 500000], na.rm = TRUE), 3), "\n")
  cat("Mean r² at 0.5-1 Mb:", round(mean(ld_data$R2[ld_data$DIST >= 500000 & ld_data$DIST < 1000000], na.rm = TRUE), 3), "\n")
  cat("Mean r² at 1-2 Mb:", round(mean(ld_data$R2[ld_data$DIST >= 1000000 & ld_data$DIST < 2000000], na.rm = TRUE), 3), "\n")
  
  # Return results
  return(list(
    chr = chr,
    binned_data = binned_data,
    decay_distance = decay_dist,
    raw_data = ld_data,
    n_pairs = nrow(ld_data)
  ))
}

# ===== COMBINED ANALYSIS FUNCTION =====

run_all_chromosomes <- function(data_dir, output_dir, 
                                chromosomes = 1:12,
                                bin_size = 10000,
                                max_dist = 5000000,
                                file_pattern = "ld_SL4.0ch%02d.ld.gz") {
  
  cat("\n" , rep("=", 60), "\n", sep = "")
  cat("GENOME-WIDE LD DECAY ANALYSIS\n")
  cat(rep("=", 60), "\n\n", sep = "")
  
  # Create output directory if it doesn't exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Initialize results list
  all_results <- list()
  
  # Process each chromosome
  for (chr in chromosomes) {
    result <- analyze_chromosome_ld(chr, data_dir, output_dir, 
                                    bin_size, max_dist, file_pattern)
    if (!is.null(result)) {
      all_results[[paste0("chr", chr)]] <- result
    }
  }
  
  # ===== COMBINED ANALYSIS =====
  cat("\n" , rep("=", 60), "\n", sep = "")
  cat("COMBINED GENOME-WIDE ANALYSIS\n")
  cat(rep("=", 60), "\n\n")
  
  # Extract decay distances
  decay_distances <- sapply(all_results, function(x) x$decay_distance)
  valid_decays <- decay_distances[!is.na(decay_distances)]
  
  cat("--- Decay Distance Summary (r² < 0.2) ---\n")
  if (length(valid_decays) > 0) {
    cat("Chromosomes with decay:\n")
    for (i in which(!is.na(decay_distances))) {
      cat(sprintf("  Chr %d: %.1f kb\n", all_results[[i]]$chr, decay_distances[i]/1000))
    }
    cat("\nSummary statistics:\n")
    cat("  Mean:", round(mean(valid_decays)/1000, 1), "kb\n")
    cat("  Median:", round(median(valid_decays)/1000, 1), "kb\n")
    cat("  Min:", round(min(valid_decays)/1000, 1), "kb\n")
    cat("  Max:", round(max(valid_decays)/1000, 1), "kb\n")
  } else {
    cat("No chromosomes showed decay below r² = 0.2\n")
  }
  
  # Combine all binned data for genome-wide plot
  cat("\n--- Creating genome-wide average plot ---\n")
  combined_binned <- do.call(rbind, lapply(names(all_results), function(chr_name) {
    all_results[[chr_name]]$binned_data %>%
      mutate(chromosome = all_results[[chr_name]]$chr)
  }))
  
  # Calculate genome-wide average
  genome_wide <- combined_binned %>%
    group_by(distance) %>%
    summarize(
      mean_r2 = mean(mean_r2, na.rm = TRUE),
      median_r2 = mean(median_r2, na.rm = TRUE),
      n_chr = n(),
      .groups = 'drop'
    )
  
  # Plot genome-wide average
  plot_ld_decay(genome_wide, 
                file.path(output_dir, "ld_decay_genome_wide.png"),
                title = "Genome-Wide Average LD Decay (All Chromosomes)")
  
  # Plot all chromosomes together
  cat("Creating multi-chromosome comparison plot...\n")
  p_multi <- ggplot(combined_binned, aes(x = distance/1000, y = mean_r2, 
                                         color = factor(chromosome))) +
    geom_line(linewidth = 0.8, alpha = 0.7) +
    labs(
      title = "LD Decay - All Chromosomes",
      x = "Distance (kb)",
      y = expression(paste("Mean ", r^2)),
      color = "Chromosome"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
      axis.title = element_text(size = 14),
      axis.text = element_text(size = 12),
      legend.position = "right"
    ) +
    geom_hline(yintercept = 0.2, linetype = "dashed", color = "red") +
    scale_color_viridis_d()
  
  ggsave(file.path(output_dir, "ld_decay_all_chromosomes.png"), 
         p_multi, width = 12, height = 6, dpi = 300)
  
  # Find genome-wide decay distance
  genome_decay <- find_decay_distance(genome_wide, threshold = 0.2)
  
  cat("\n" , rep("=", 60), "\n", sep = "")
  cat("GWAS RECOMMENDATIONS\n")
  cat(rep("=", 60), "\n\n")
  
  if (!is.na(genome_decay)) {
    cat("Genome-wide LD decay distance (r² < 0.2):", round(genome_decay/1000, 1), "kb\n\n")
    cat("RECOMMENDED HAPLOTYPE BLOCK SETTINGS:\n")
    cat("  Conservative (short blocks):", round(genome_decay/1000 * 0.5, 1), "kb\n")
    cat("  Standard (at decay point):", round(genome_decay/1000, 1), "kb\n")
    cat("  Liberal (extended blocks):", round(genome_decay/1000 * 1.5, 1), "kb\n")
  } else {
    # Use median of valid decays
    if (length(valid_decays) > 0) {
      recommended <- median(valid_decays)
      cat("Genome-wide average does not decay, but individual chromosomes do.\n")
      cat("Using median chromosome decay:", round(recommended/1000, 1), "kb\n\n")
      cat("RECOMMENDED HAPLOTYPE BLOCK SETTINGS:\n")
      cat("  Conservative (short blocks):", round(recommended/1000 * 0.5, 1), "kb\n")
      cat("  Standard (at decay point):", round(recommended/1000, 1), "kb\n")
      cat("  Liberal (extended blocks):", round(recommended/1000 * 1.5, 1), "kb\n")
    } else {
      cat("WARNING: High LD persists across all chromosomes!\n")
      cat("This may indicate:\n")
      cat("  - High population structure\n")
      cat("  - Small effective population size\n")
      cat("  - Recent bottleneck\n")
      cat("  - Self-pollinating species\n\n")
      cat("RECOMMENDED ACTIONS:\n")
      cat("  1. Check population structure (PCA, admixture)\n")
      cat("  2. Consider using 500 kb - 1 Mb blocks conservatively\n")
      cat("  3. May need to account for structure in GWAS model\n")
    }
  }
  
  # Create summary table
  summary_table <- data.frame(
    Chromosome = sapply(all_results, function(x) x$chr),
    SNP_Pairs = sapply(all_results, function(x) x$n_pairs),
    Decay_Distance_kb = round(decay_distances/1000, 1),
    Mean_r2_100_500kb = sapply(all_results, function(x) {
      round(mean(x$raw_data$R2[x$raw_data$DIST >= 100000 & 
                                 x$raw_data$DIST < 500000], na.rm = TRUE), 3)
    })
  )
  
  cat("\n--- Summary Table ---\n")
  print(summary_table)
  
  # Save summary table
  write.csv(summary_table, 
            file.path(output_dir, "ld_decay_summary.csv"), 
            row.names = FALSE)
  
  cat("\nSummary table saved to:", file.path(output_dir, "ld_decay_summary.csv"), "\n")
  
  # Return all results
  return(list(
    chromosome_results = all_results,
    genome_wide = genome_wide,
    summary_table = summary_table,
    decay_distances = decay_distances
  ))
}

# ===== EXAMPLE USAGE =====

# Set your paths
# data_dir <- "/Users/cperkins/Desktop/palanivelu_lab/diversity_panel/LD/data"
# output_dir <- "/Users/cperkins/Desktop/palanivelu_lab/diversity_panel/LD/results"

# Run analysis for all chromosomes
# results <- run_all_chromosomes(
#   data_dir = data_dir,
#   output_dir = output_dir,
#   chromosomes = 1:12,
#   bin_size = 10000,      # 10 kb bins
#   max_dist = 5000000,    # Analyze up to 5 Mb
#   file_pattern = "ld_SL4.0ch%02d.ld.gz"  # Adjust pattern to match your files
# )

# ===== QUICK START FUNCTION =====

quick_ld_analysis <- function(data_dir, output_dir = NULL) {
  # If no output directory specified, create one in data directory
  if (is.null(output_dir)) {
    output_dir <- file.path(dirname(data_dir), "LD_analysis_results")
  }
  
  # Set working directory
  setwd(dirname(data_dir))
  
  # Run analysis
  results <- run_all_chromosomes(
    data_dir = data_dir,
    output_dir = output_dir,
    chromosomes = 1:12,
    bin_size = 10000,
    max_dist = 5000000,
    file_pattern = "ld_SL4.0ch%02d.ld.gz"
  )
  
  cat("\n✓ Analysis complete! Results saved to:", output_dir, "\n")
  
  return(results)
}

# ===== TO RUN THE ANALYSIS =====
# Uncomment and modify these lines:

setwd("/Users/cperkins/Desktop/palanivelu_lab/diversity_panel/LD")
results <- quick_ld_analysis(
  data_dir = "/Users/cperkins/Desktop/palanivelu_lab/diversity_panel/LD/results"
)
