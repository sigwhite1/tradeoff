library(vcfR)
library(adegenet)
library(ggplot2)
library(tidyr)
library(dplyr)
library(tidyverse)
library(reshape2)
library(ape)
library(clusterProfiler)  # For gene set enrichment
library(biomaRt)
library(SNPRelate)
library(wordcloud)
library(tm)
library(phangorn)
library(patchwork)
library(gridExtra)
library(grid)

###################Figure 4 - heatmap###########################
setwd("C:/Users/lojoh/OneDrive/Documents/BootsLab/moths/data") ##Change this to your working directory

vars <- read.delim("inbred_intersect_cds.out", header=FALSE, stringsAsFactors = FALSE)
vars <- vars %>%
  # Extract the gene name (delimited by `gene=...;`)
  mutate(gene_name = str_extract(vars[[30]], "(?<=gene=)[^;]+")) %>%
  # Extract the product description (delimited by `product=...`)
  mutate(product_description = str_extract(vars[[30]], "(?<=product=)[^;]+")) 
# vars[, 30] <- gsub(".*gene=([A-Za-z0-9_.]+).*", "\\1", vars[, 30])
vars <- vars[, -c(3,7,8,22, 27, 28, 29, 30)]
colnames(vars) <- c("chrom","pos","ref", "alt", "qual", "format", "IL-1", "IL-24", "IL-29", 
                    "IL-36", "IL-2", "IL-4", "IL-7", "IL-9", "IL-10","IL-14","IL-17","IL-19",
                    "source", "type", "start", "end", "gene", "annot")

geno_mat <- vars[,7:18]

clean_genotype <- function(geno) {
  gsub(":.*", "", geno)  
}

# Apply cleaning function to genotype matrix
geno_mat_clean <- as.matrix(apply(geno_mat, 2, clean_genotype))


genotype_numeric <- geno_mat_clean
genotype_numeric[geno_mat_clean == "./."] <- 0
genotype_numeric[geno_mat_clean == "0/1"] <- 1
genotype_numeric[geno_mat_clean == "1/1"] <- 2

genotype_numeric <- apply(genotype_numeric, 2, as.numeric)

genotype_numeric <- cbind(genotype_numeric, rep(0, nrow(genotype_numeric)))

colnames(genotype_numeric)[ncol(genotype_numeric)] <- "ref"

pairwise_distances <- dist.gene(t(genotype_numeric), 
                                method = "pairwise", pairwise.deletion = TRUE)

# Convert to a table format
pairwise_table <- as.matrix(pairwise_distances)



# Reshape the distance matrix for ggplot
distance_melted <- melt(as.matrix(pairwise_table))

# Rename columns for clarity
colnames(distance_melted) <- c("Inbred_1", "Inbred_2", "Distance")

# Plot the heatmap
ggplot(distance_melted, aes(x = Inbred_1, y = Inbred_2, fill = Distance)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "red") +  # Adjust colors as needed
  theme_minimal() +
  labs(
    title = "Pairwise Genetic Distance Heatmap",
    x = "Inbred Line",
    y = "Inbred Line",
    fill = "Distance"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),  # Rotate x-axis labels
    axis.text.y = element_text(size = 10)
  )
#####################################################

#################figure 6 - manhattan plots###############
setwd("C:/Users/lojoh/OneDrive/Documents/BootsLab/moths/data/snprelate") ##Change this to your working directory

vcf_files <- list.files("windows/", pattern = "\\.vcf$", full.names = TRUE, recursive = TRUE)

# Exclude NW contigs
vcf_files <- vcf_files[!grepl("^window_NW_", basename(vcf_files))]

# Output cache
dir.create("wza_cache", showWarnings = FALSE)

# LD50 phenotype table
ld50 <- data.frame(
  sample_id = c("IL-1", "IL-24", "IL-29", "IL-36", "IL-2", "IL-4", "IL-7", "IL-9", "IL-10", "IL-14", "IL-17", "IL-19"),
  ld50 = c(0.04955, 0.08315, 0.0536, 0.063, 0.0773, 0.04625, 0.08253, 0.05304, 0.06401, 0.06311, 0.08253, 0.05316),
  stringsAsFactors = FALSE
)

# Sample mapping
vcf_samples <- c(
  "/global/scratch/users/lojohnathan6/projects/moth_seq_23/data/merged_bam/I1.bam",
  "/global/scratch/users/lojohnathan6/projects/moth_seq_23/data/merged_bam/I10.bam",
  "/global/scratch/users/lojohnathan6/projects/moth_seq_23/data/merged_bam/I11.bam",
  "/global/scratch/users/lojohnathan6/projects/moth_seq_23/data/merged_bam/I12.bam",
  "/global/scratch/users/lojohnathan6/projects/moth_seq_23/data/merged_bam/I2.bam",
  "/global/scratch/users/lojohnathan6/projects/moth_seq_23/data/merged_bam/I3.bam",
  "/global/scratch/users/lojohnathan6/projects/moth_seq_23/data/merged_bam/I4.bam",
  "/global/scratch/users/lojohnathan6/projects/moth_seq_23/data/merged_bam/I5.bam",
  "/global/scratch/users/lojohnathan6/projects/moth_seq_23/data/merged_bam/I6.bam",
  "/global/scratch/users/lojohnathan6/projects/moth_seq_23/data/merged_bam/I7.bam",
  "/global/scratch/users/lojohnathan6/projects/moth_seq_23/data/merged_bam/I8.bam",
  "/global/scratch/users/lojohnathan6/projects/moth_seq_23/data/merged_bam/I9.bam"
)
ld50_ids <- c("IL-1", "IL-24", "IL-29", "IL-36", "IL-2", "IL-4", "IL-7", "IL-9", "IL-10", "IL-14", "IL-17", "IL-19")
sample_map <- setNames(ld50_ids, vcf_samples)

all_snp_data <- list()

for (vcf_path in vcf_files) {
  window_id <- basename(vcf_path)
  out_txt <- file.path("wza_cache", paste0(window_id, ".txt"))
  
  if (file.exists(out_txt)) {
    cat("Skipping cached:", window_id, "\n")
    next
  }
  
  cat("Processing:", window_id, "\n")
  vcf <- read.vcfR(vcf_path, verbose = FALSE)
  
  # Extract raw GT strings
  gt_raw <- extract.gt(vcf, element = "GT", as.numeric = FALSE)
  
  # Parse GT strings
  gt <- apply(gt_raw, c(1, 2), function(x) {
    if (is.na(x) || x == "" || x == "." || x == ".:." || x == "./.") return(0L)
    gt_part <- strsplit(x, ":")[[1]][1]
    if (gt_part %in% c("0/0", "0|0")) return(0L)
    if (gt_part %in% c("0/1", "1/0", "0|1", "1|0")) return(1L)
    if (gt_part %in% c("1/1", "1|1")) return(2L)
    return(NA_integer_)
  })
  
  # Rename columns to LD50-compatible names
  samples <- colnames(gt)
  colnames(gt) <- sample_map[samples]
  
  # Subset LD50 data to match genotype columns
  ld50_sub <- ld50[match(colnames(gt), ld50$sample_id), ]
  
  # Drop SNPs with any NA values
  gt <- gt[complete.cases(gt), , drop = FALSE]
  if (nrow(gt) == 0) next
  
  # MAF without NA checks
  maf_vec <- rowMeans(gt) / 2
  maf_vec <- pmin(maf_vec, 1 - maf_vec)
  
  # Filter SNPs by MAF >= 0.05
  keep_idx <- which(maf_vec >= 0.05)
  if (length(keep_idx) == 0) next
  
  gt <- gt[keep_idx, , drop = FALSE]
  maf_vec <- maf_vec[keep_idx]
  snp_ids <- rownames(gt)[keep_idx]
  
  # Correlation without any more checks
  cor_p <- t(apply(gt, 1, function(geno) {
    test <- cor.test(ld50_sub$ld50, geno, method = "pearson")
    c(cor = test$estimate, p = test$p.value)
  }))
  
  cor_p <- as.data.frame(cor_p)
  cor_p$snp_id <- snp_ids
  cor_p$maf <- maf_vec
  cor_p$window <- window_id
  
  write.table(cor_p, file = out_txt, sep = "\t", row.names = FALSE, quote = FALSE)
  all_snp_data[[vcf_path]] <- cor_p
}

# Load all cached SNP-level data
cached_files <- list.files("wza_cache", pattern = "*.txt$", full.names = TRUE)
snp_data <- bind_rows(lapply(cached_files, read.delim))

# WZA calculations
snp_data$empirical_p <- rank(snp_data$p, ties.method = "average") / nrow(snp_data)
snp_data$z <- qnorm(1 - snp_data$empirical_p)
snp_data$He <- snp_data$maf * (1 - snp_data$maf)

wza_results <- snp_data %>%
  group_by(window) %>%
  summarise(
    n_SNPs = n(),
    He_sq_sum = sum(He^2),
    ZW = sum(He * z) / sqrt(He_sq_sum)
  ) %>%
  filter(is.finite(ZW))  # filters out -Inf, Inf, NaN


# Normalize ZW by SNP count
mean_model <- lm(ZW ~ poly(n_SNPs, 2), data = wza_results)
sd_model   <- lm(abs(ZW - predict(mean_model)) ~ poly(n_SNPs, 2), data = wza_results)

wza_results$meanZW  <- predict(mean_model)
wza_results$sdZW    <- predict(sd_model)
wza_results$ZW_norm <- (wza_results$ZW - wza_results$meanZW) / wza_results$sdZW
wza_results$p_value <- 2 * (1 - pnorm(abs(wza_results$ZW_norm)))
wza_results$log10p  <- -log10(wza_results$p_value)

# Final WZA results
print(head(wza_results))
qq_data <- wza_results %>%
  filter(is.finite(p_value), !is.na(p_value)) %>%
  arrange(p_value) %>%
  mutate(
    observed = -log10(p_value),
    expected = -log10(ppoints(n()))
  )

# Plot
ggplot(qq_data, aes(x = expected, y = observed)) +
  geom_point(size = 1.5, alpha = 0.7, color = "blue") +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
  theme_minimal() +
  labs(
    title = "WZA QQ",
    x = "Expected -log10(P)",
    y = "Observed -log10(P)"
  )

manhattan_data <- wza_results %>%
  mutate(
    chrom = sub("window_((NC|NW|CM)_[0-9]+\\.[0-9]+)_.*", "\\1", window),
    position = as.numeric(sub(".*_([0-9]+)_([0-9]+)\\.vcf$", "\\1", window)),
    end_position = as.numeric(sub(".*_([0-9]+)_([0-9]+)\\.vcf$", "\\2", window)),
    chrom_index = as.numeric(factor(chrom, levels = unique(chrom))),
    color = as.factor(chrom_index %% 2),
    logP = -log10(p_value)
  )

# Optional: define a significance threshold (e.g., FDR-adjusted or Bonferroni)
fdr_threshold_logP <- -log10(0.05 / nrow(manhattan_data))

# Manhattan plot
# ggplot(manhattan_data, aes(x = position, y = logP, color = color)) +
#   geom_point(alpha = 0.75, size = 1.5) +
#   geom_hline(yintercept = fdr_threshold_logP, linetype = "dashed", color = "red", size = 1) +
#   facet_grid(~ chrom, scales = "free_x", space = "free_x") +
#   scale_color_manual(values = c("gray30", "gray60"), guide = "none") +
#   labs(
#     title = "WZA scores by chromosome",
#     x = "Chromosome",
#     y = "-log10(WZA P-value)"
#   ) +
#   theme_minimal() +
#   theme(
#     legend.position = "none",
#     strip.text = element_text(size = 8, angle = 45, hjust = 1, vjust = 1),
#     axis.text.x = element_blank(),
#     axis.ticks.x = element_blank()
#   )

max_window_pos <- manhattan_data %>%
  group_by(chrom) %>%
  summarise(max_pos = max(end_position, na.rm = TRUE)) %>%
  ungroup()

global_max_pos <- max(max_window_pos$max_pos)

# Add fake invisible points at global_max_pos per chrom to enforce equal width
padding_rows <- manhattan_data %>%
  group_by(chrom) %>%
  summarise() %>%
  mutate(
    position = global_max_pos,
    logP = NA,
    color = "gray",
    window = NA,
    end_position = global_max_pos,
    chrom_index = as.numeric(factor(chrom, levels = unique(manhattan_data$chrom)))
  )

# Bind with original data
manhattan_plot_data <- bind_rows(manhattan_data, padding_rows)

manhattan_plot_data <- manhattan_plot_data %>%
  mutate(
    color = as.factor(ifelse(is.na(color), "padding", as.character(color)))
  )


# Plot
p_wza <- ggplot(manhattan_plot_data, aes(x = position, y = logP, color = color)) +
  geom_point(alpha = 0.75, size = 1.5, na.rm = TRUE) +
  geom_hline(yintercept = fdr_threshold_logP, linetype = "dashed", color = "red", size = 1) +
  facet_grid(~ chrom, scales = "fixed", space = "fixed") +
  scale_color_manual(
    values = c("0" = "gray30", "1" = "gray60", "padding" = NA),
    guide = "none"
  )+
  labs(
    title = "WZA scores by chromosome",
    x = "Chromosome",
    y = "-log10(WZA P-value)"
  ) +
  theme_minimal() +
  theme(
    legend.position = "none",
    strip.text = element_text(size = 8, angle = 45, hjust = 1, vjust = 1),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )

#################snprelate PCA###############################

# List VCF files
vcf_files <- list.files("windows/", pattern = "*.vcf$", full.names = TRUE, recursive = TRUE)
vcf_files <- vcf_files[!grepl("_NW_", vcf_files)]

# Initialize an empty data frame for storing PCA results
ordered <- NULL
output_file <- "snprelate_1mb.csv"

# Function to write results to CSV
write_results_to_csv <- function(data, file, append = FALSE) {
  write.table(data, file, append = append, sep = ",", col.names = !append, row.names = FALSE)
}

# Loop through VCF files, convert to GDS, and perform PCA
for (i in 1:length(vcf_files)) {
  vcf.fn <- vcf_files[i]
  print(paste(i, "Processing file:", vcf.fn))
  
  # Extract file information
  only_file_name <- basename(vcf.fn)  # Use basename for file name extraction
  print(paste("File name:", only_file_name))
  
  dataset <- only_file_name
  print(paste("Dataset:", dataset))
  
  # Convert VCF to GDS
  gds.fn <- sub("\\.vcf$", ".gds", vcf.fn)  # Replace .vcf with .gds
  if(!file.exists(gds.fn)){
    snpgdsVCF2GDS(vcf.fn, gds.fn, method = "biallelic.only")
    print(paste("Converted VCF to GDS:", gds.fn))
  }
  
  # Open GDS file
  genofile <- tryCatch({
    snpgdsOpen(gds.fn)
  }, error = function(e) {
    message(paste("Error opening GDS file:", gds.fn, ":", e$message))
    return(NULL)
  })
  
  if (is.null(genofile)) {
    next
  }
  
  print("GDS file opened successfully")
  
  # Perform PCA
  pca_window <- tryCatch({
    snpgdsPCA(genofile, autosome.only = FALSE)
  }, error = function(e) {
    message(paste("Error performing PCA for", gds.fn, ":", e$message))
    snpgdsClose(genofile)
    return(NULL)
  })
  
  if (is.null(pca_window)) {
    next
  }
  
  print("PCA completed successfully")
  
  # Extract sample IDs and PCA results
  sample_ids <- pca_window$sample.id
  pca_values <- as.double(pca_window$eigenvect[, 1])
  
  # Initialize the data frame with sample IDs if it's the first successful iteration
  if (is.null(ordered)) {
    ordered <- data.frame(sample_id = sample_ids)
  }
  
  # # Ensure the sample IDs match
  # if (!all(ordered$sample_id == sample_ids)) {
  #   warning(paste("Sample IDs do not match for", gds.fn, "- Skipping this file"))
  #   snpgdsClose(genofile)
  #   next
  # }
  
  # Add PCA results to the data frame
  column_name <- paste("slidingwindow_", dataset, sep = "")
  ordered[[column_name]] <- pca_values
  
  # Write the updated data frame to CSV
  write_results_to_csv(ordered, output_file, append = FALSE)
  
  cat("PCA values appended to data frame \n", file="snprelate.log", append=T)
  
  # Close the GDS file
  snpgdsClose(genofile)
}

print("Processing completed.")

sample_ids_converted <- c("IL-1", "IL-24", "IL-29", "IL-36", "IL-2", "IL-4", "IL-7", "IL-9", "IL-10","IL-14","IL-17","IL-19")

ordered$sample_id <- sample_ids_converted

ld50_table <- cbind(sample_ids_converted, c(.04955, 0.08315, .0536, .063, .0773, .04625, .08253, .05304, .06401, .06311, .08253, .05316))
colnames(ld50_table) <- c("sample_id", "ld50")

# Merge by Sample ID
merged_data <- merge(ordered, ld50_table, by = "sample_id")

# Move LD50 column to the first column for clarity
merged_data <- merged_data %>%
  dplyr::select(sample_id, ld50, everything())

# Remove Sample ID column after merging
rownames(merged_data) <- merged_data$sample_id
merged_data <- merged_data[,-1]
merged_data$ld50 <- as.numeric(merged_data$ld50)

cor_values <- numeric(ncol(merged_data)-1)
p_values <- numeric(ncol(merged_data)-1)

# Loop through each PCA score column and compute correlation with LD50
for (i in 1:(ncol(merged_data)-1)) {
  test_result <- cor.test(merged_data$ld50, merged_data[, i+1], method = "pearson")
  cor_values[i] <- test_result$estimate  # Pearson correlation coefficient
  p_values[i] <- test_result$p.value     # P-value
}

# Create a results dataframe
cor_results <- data.frame(window = colnames(merged_data)[-1], 
                          correlation = cor_values, 
                          p_value = p_values)

# Print first few results
head(cor_results)

cor_results <- cor_results[-ncol(merged_data),]

cor_results$fdr_adj <- p.adjust(cor_results$p_value, method = "fdr")

cor_results <- na.omit(cor_results)
alpha <- 0.05  # FDR level
sorted_pvalues <- sort(cor_results$p_value)  # Sort p-values in ascending order
threshold_index <- max(which(sorted_pvalues <= (alpha * (1:nrow(cor_results)) / nrow(cor_results))))  # Find largest significant index

fdr_threshold_logP <- -log10(alpha/nrow(cor_results))

# Create a dataframe for plotting
observed_pvals_log <- -log10(sort(cor_results$p_value))
expected_pvals <- -log10(ppoints(length(observed_pvals_log)))

# Create a dataframe for plotting
qq_data <- data.frame(Expected = expected_pvals, Observed = observed_pvals_log)

# Generate QQ plot using ggplot2
ggplot(qq_data, aes(x = Expected, y = Observed)) +
  geom_point(alpha = 0.7, size = 2, color = "blue") +  # Scatter points
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +  # Reference line
  theme_minimal() +
  labs(title = "1 Mb sliding window QQ",
       x = "Expected -log10(P)",
       y = "Observed -log10(P)") +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10))

#######chromplots#########
cor_results <- cor_results %>%
  mutate(chrom = sub("slidingwindow_window_((NC|NW)_[0-9]+\\.[0-9]+)_.*", "\\1", window))


# Extract Start Position
cor_results <- cor_results %>%
  mutate(position = as.numeric(sub(".*_([0-9]+)_([0-9]+)\\.vcf$", "\\1", window)))

# Extract End Position (Optional)
cor_results <- cor_results %>%
  mutate(end_position = as.numeric(sub(".*_([0-9]+)_([0-9]+)\\.vcf$", "\\2", window)))

cor_results <- cor_results %>%
  mutate(chrom_index = as.numeric(factor(chrom, levels = unique(chrom)))) %>%  # Assign numeric index to each chromosome
  mutate(color = as.factor(chrom_index %% 2)) 

cor_results <- cor_results %>%
  mutate(logP = -log10(p_value))


ggplot(cor_results, aes(x = position, y = logP, color = logP)) +
  geom_point(alpha = 0.75, size = 1.5) +  # Points based on Position
  geom_hline(yintercept = fdr_threshold_logP, linetype = "dashed", color = "red", size = 1) +  # Add FDR threshold line
  facet_grid(~ chrom, scales = "free_x", space = "free_x") +  # Group by chromosome, keep positions spread out
  scale_color_gradient(low = "blue", high = "red") +  # Color by P-value (low = blue, high = red)
  labs(title = "Manhattan Plot of PCA-LD50 Correlation, all SNPs",
       x = "Chromosome",
       y = "-log10(P-value)",
       color = "Significance") +  # Legend label
  theme_minimal() +
  theme(
    legend.position = "right",  # Keep legend for significance
    strip.text = element_text(size = 8, angle = 45, hjust = 1, vjust = 1),  # Rotate chromosome names
    axis.text.x = element_blank(),  # Remove position labels
    axis.ticks.x = element_blank()  # Remove x-axis ticks
  )

pca_manhattan_data <- cor_results %>%
  mutate(
    chrom = sub("slidingwindow_window_((NC|NW)_[0-9]+\\.[0-9]+)_.*", "\\1", window),
    position = as.numeric(sub(".*_([0-9]+)_([0-9]+)\\.vcf$", "\\1", window)),
    end_position = as.numeric(sub(".*_([0-9]+)_([0-9]+)\\.vcf$", "\\2", window)),
    chrom_index = as.numeric(factor(chrom, levels = unique(chrom))),
    color = as.factor(chrom_index %% 2),
    logP = -log10(p_value)
  )

# Compute max end_position per chromosome
max_window_pos <- pca_manhattan_data %>%
  group_by(chrom) %>%
  summarise(max_pos = max(end_position, na.rm = TRUE), .groups = "drop")

global_max_pos <- max(max_window_pos$max_pos)

# Generate padding rows
padding_rows <- pca_manhattan_data %>%
  distinct(chrom, chrom_index) %>%
  mutate(
    position = global_max_pos,
    end_position = global_max_pos,
    logP = NA,
    color = "padding",
    window = NA
  )
pca_plot_data <- bind_rows(pca_manhattan_data, padding_rows)

pca_plot_data <- pca_plot_data %>%
  mutate(
    color = factor(color, levels = c("0", "1", "padding"))
  )
p_pca <- ggplot(pca_plot_data, aes(x = position, y = logP, color = color)) +
  geom_point(alpha = 0.75, size = 1.5, na.rm = TRUE) +
  geom_hline(yintercept = fdr_threshold_logP, linetype = "dashed", color = "red", size = 1) +
  facet_grid(~ chrom, scales = "fixed", space = "fixed") +
  scale_color_manual(
    values = c("0" = "gray30", "1" = "gray60", "padding" = NA),
    guide = "none"
  ) +
  labs(
    title = "Manhattan Plot of PCA-LD50 Correlation, all SNPs",
    x = "Chromosome",
    y = expression(-log[10](italic(P)))
  ) +
  theme_minimal() +
  theme(
    legend.position = "none",
    strip.text = element_text(size = 8, angle = 45, hjust = 1, vjust = 1),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )


#####################combined wza/pca plot##################
combined_plot <- p_pca / p_wza + plot_layout(heights = c(1, 1))
combined_plot

##################table of windows###########################
# Assume wza_results and cor_results (PCA) are loaded and processed

# Normalize window names to align between tables
cor_results$clean_window <- gsub("slidingwindow_window_", "", cor_results$window)
wza_results$clean_window <- gsub("^window_", "", wza_results$window)

# Apply FDR correction
wza_results$wza_fdr <- p.adjust(wza_results$p_value, method = "fdr")
cor_results$pca_fdr <- p.adjust(cor_results$p_value, method = "fdr")

# Rank each by significance
wza_results$wza_rank <- rank(wza_results$p_value, ties.method = "min")
cor_results$pca_rank <- rank(cor_results$p_value, ties.method = "min")

# Merge on cleaned window name
combined_df <- merge(
  wza_results[, c("clean_window", "wza_fdr", "wza_rank")],
  cor_results[, c("clean_window", "pca_fdr", "pca_rank")],
  by = "clean_window"
)

# Get top 10 windows by WZA
top10_combined <- combined_df[order(combined_df$wza_rank), ][1:10, ]

# View or write table
print(top10_combined)
# write.csv(top10_combined, "top10_combined_table.csv", row.names = FALSE)

# Format for display
top10_combined_display <- top10_combined
colnames(top10_combined_display) <- c("Window", "WZA FDR", "WZA Rank", "PCA FDR", "PCA Rank")

# Round numeric values for aesthetics
top10_combined_display$`WZA FDR` <- signif(top10_combined_display$`WZA FDR`, 3)
top10_combined_display$`PCA FDR` <- signif(top10_combined_display$`PCA FDR`, 3)

# Create the table grob
table_grob <- tableGrob(top10_combined_display, rows = NULL)

# Draw the table with a title
grid.newpage()
grid.draw(table_grob)
grid.text("Top 10 Genomic Windows by WZA",
          y = unit(1, "npc") - unit(2, "mm"),
          gp = gpar(fontsize = 14, fontface = "bold"))
