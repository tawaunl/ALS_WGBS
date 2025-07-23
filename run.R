# run_analysis.R
# Main script to perform the cell-free DNA analysis.

# Load necessary libraries (already handled by setup.R, but good for standalone execution)
# You can uncomment these if you want to run this file directly without setup.R first
# library(Rsamtools)
# library(ggplot2)
# library(dplyr)
# library(tidyr)
# library(moments)
# library(GenomicAlignments)
# library(GenomicRanges)
# library(rtracklayer)
# library(GenomicFeatures)
# library(TxDb.Hsapiens.UCSC.hg38.knownGene)
# library(BSgenome.Hsapiens.UCSC.hg38)
# library(ggseqlogo)
# library(caret)
# library(pROC)

# Source the functions script
source("functions.R")

# Set up directories
# Make sure to adjust these paths based on where you place your data and figures
data.dir <- "../celfie" # Assuming 'celfie' folder is parallel to 'Analysis'
fig.dir <- "../Figures" # Assuming 'Figures' folder is parallel to 'Analysis'

# Create figures directory if it doesn't exist
if (!dir.exists(fig.dir)) {
  dir.create(fig.dir, recursive = TRUE)
}

message("Starting cfDNA analysis...")

# Read Metadata
meta <- read.csv(file.path(data.dir, "celfie_cfDNA_ss.csv"))
ctrl.samples <- meta$Run[which(meta$disease_status == "ctrl")]
als.sample <- meta$Run[which(meta$disease_status == "als")]

bam_files_ctrl <- c()
for (name in ctrl.samples) {
  bam_files_ctrl <- c(bam_files_ctrl, list.files(path = file.path(data.dir, "bam"),
                                                 pattern = paste0("^", name, ".*\\.bam$"), recursive = TRUE))
}
bam_files_als <- c()
for (name in als.sample) {
  bam_files_als <- c(bam_files_als, list.files(path = file.path(data.dir, "bam"),
                                               pattern = paste0("^", name, ".*\\.bam$"), recursive = TRUE))
}

# Combine for easier iteration and plotting
all_bam_files <- c(bam_files_ctrl, bam_files_als)
sample_groups <- c(rep("Control", length(bam_files_ctrl)),
                   rep("ALS", length(bam_files_als)))
sample_names <- c(ctrl.samples, als.sample)

message("--- Insert Size Analysis ---")
# Insert Analysis
max_insert_size <- 1000
min_insert_size <- 0

all_insert_sizes_list <- list()

for (i in seq_along(all_bam_files)) {
  bam_path <- file.path(data.dir, "bam", all_bam_files[i])
  sample_name <- sample_names[i]
  group <- sample_groups[i]
  
  sizes <- get_insert_sizes(bam_path, min_len = min_insert_size, max_len = max_insert_size)
  
  all_insert_sizes_list[[sample_name]] <- tibble(
    InsertSize = sizes,
    Sample = sample_name,
    Group = group
  )
}

all_insert_sizes_df <- bind_rows(all_insert_sizes_list)

# Generate Frequency Plot (Histogram) by sample
pdf(file.path(fig.dir, "InsertSizeDensity_bySample.pdf"), width = 10, height = 10)
ggplot(all_insert_sizes_df, aes(x = InsertSize, fill = Group)) +
  geom_histogram(binwidth = 1, alpha = 0.7, position = "identity") +
  facet_wrap(~ Sample, scales = "free_y", ncol = 2) +
  labs(title = "Insert Size Distribution per Sample",
       x = "Insert Size (bp)",
       y = "Frequency") +
  theme_minimal() +
  theme(legend.position = "bottom")
dev.off()
message("Generated InsertSizeDensity_bySample.pdf")

# Plotting by group
pdf(file.path(fig.dir, "InsertSizeDensity_byGroup.pdf"))
ggplot(all_insert_sizes_df, aes(x = InsertSize, color = Group, fill = Group)) +
  geom_density(alpha = 0.3) +
  labs(title = "Combined Insert Size Distribution by Group",
       x = "Insert Size (bp)",
       y = "Density") +
  theme_minimal() +
  theme(legend.position = "bottom") +
  xlim(min_insert_size, max_insert_size)
dev.off()
message("Generated InsertSizeDensity_byGroup.pdf")


# Analyze mono-nucleosomal and di-nucleosomal fragments
mono_min <- 150
mono_max <- 180
di_min <- 300
di_max <- 380

df_fragment_types <- all_insert_sizes_df %>%
  mutate(
    FragmentType = case_when(
      InsertSize >= mono_min & InsertSize <= mono_max ~ "Mono-nucleosomal",
      InsertSize >= di_min & InsertSize <= di_max ~ "Di-nucleosomal",
      TRUE ~ "Other_Fragment"
    )
  )

fragment_counts_proportions <- df_fragment_types %>%
  group_by(Sample, Group, FragmentType) %>%
  summarise(Count = n(), .groups = 'drop_last') %>%
  mutate(Proportion = Count / sum(Count)) %>%
  ungroup()

mono_di_proportions <- fragment_counts_proportions %>%
  filter(FragmentType %in% c("Mono-nucleosomal", "Di-nucleosomal")) %>%
  group_by(Sample, Group) %>%
  summarise(
    Mono_Count = sum(Count[FragmentType == "Mono-nucleosomal"]),
    Di_Count = sum(Count[FragmentType == "Di-nucleosomal"]),
    Total_MonoDi = Mono_Count + Di_Count,
    Mono_Proportion = Mono_Count / Total_MonoDi,
    Di_Proportion = Di_Count / Total_MonoDi,
    Mono_Di_Ratio = Mono_Count / Di_Count,
    .groups = 'drop'
  ) %>%
  mutate(Mono_Di_Ratio = ifelse(is.infinite(Mono_Di_Ratio), NA, Mono_Di_Ratio))

wilcox_mono_prop <- wilcox.test(Mono_Proportion ~ Group, data = mono_di_proportions)
message("\nWilcoxon Rank-Sum Test for Mono-nucleosomal Proportion difference:")
print(wilcox_mono_prop)

# Box plot for Mono_Proportion with outliers flagged
mono_di_proportions_flagged <- identify_outliers_iqr(
  data = mono_di_proportions,
  variable_col = Mono_Proportion,
  group_col = Group,
  k = 1.5
)
ggplot(mono_di_proportions_flagged, aes(x = Group, y = Mono_Proportion, fill = Group)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_point(aes(color = is_outlier),
             position = position_jitter(width = 0.1),
             size = 2,
             alpha = 0.8) +
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "black"),
                     labels = c("TRUE" = "Outlier", "FALSE" = "Not Outlier")) +
  labs(title = "Mono-nucleosomal Fragment Proportion with Outliers Highlighted",
       x = "Group",
       y = "Proportion of Fragments (140-180bp)",
       color = "Data Point Type") +
  theme_minimal() +
  theme(legend.position = "bottom",
        plot.title = element_text(hjust = 0.5))
ggsave(file.path(fig.dir, "MononucleosomalProportion_outliers_boxplot.pdf"), width = 10, height = 10)
message("Generated MononucleosomalProportion_outliers_boxplot.pdf")


# Remove outliers for further statistical tests and clean plots
mono_di_proportions_cleaned <- remove_outliers_iqr(
  data = mono_di_proportions,
  variable_col = Mono_Proportion,
  group_col = Group,
  k = 1.5
)

wilcox_mono_prop_cleaned <- wilcox.test(Mono_Proportion ~ Group, data = mono_di_proportions_cleaned)
message("\nWilcoxon Rank-Sum Test for Mono-nucleosomal Proportion difference (after outlier removal):")
print(wilcox_mono_prop_cleaned)

pdf(file.path(fig.dir, "MononucleosomalProportionCLEAN_boxplot.pdf"), width = 10, height = 10)
ggplot(mono_di_proportions_cleaned, aes(x = Group, y = Mono_Proportion, fill = Group)) +
  geom_boxplot(alpha = 0.7) +
  geom_point(position = position_jitter(width = 0.1), alpha = 0.8, color = "black") +
  labs(title = "Mono-nucleosomal Fragment Proportion by Group (Outliers Removed)",
       x = "Group",
       y = "Proportion of Fragments (150-180bp)") +
  theme_minimal() +
  theme(legend.position = "none")
dev.off()
message("Generated MononucleosomalProportionCLEAN_boxplot.pdf")

wilcox_mono_di_ratio <- wilcox.test(Mono_Di_Ratio ~ Group, data = mono_di_proportions_cleaned)
message("\nWilcoxon Rank-Sum Test for Mono/Di-nucleosomal Ratio difference:")
print(wilcox_mono_di_ratio)

pdf(file.path(fig.dir, "Mono_Di_Ratio_boxplot.pdf"), width = 10, height = 10)
ggplot(mono_di_proportions_cleaned, aes(x = Group, y = Mono_Di_Ratio, fill = Group)) +
  geom_boxplot(alpha = 0.7) +
  geom_point(position = position_jitter(width = 0.1), alpha = 0.8, color = "black") +
  labs(title = "Mono/Di-nucleosomal Fragment Ratio by Group",
       x = "Group",
       y = "Ratio of Fragments (Mono- / Di-nucleosomal)") +
  theme_minimal() +
  theme(legend.position = "none")
dev.off()
message("Generated Mono_Di_Ratio_boxplot.pdf")


message("--- Distribution Tests (Skewness, Kurtosis) ---")
# Distribution tests
sample_moments <- all_insert_sizes_df %>%
  group_by(Sample, Group) %>%
  summarise(
    Mean_InsertSize = mean(InsertSize, na.rm = TRUE),
    Median_InsertSize = median(InsertSize, na.rm = TRUE),
    SD_InsertSize = sd(InsertSize, na.rm = TRUE),
    Skewness_InsertSize = moments::skewness(InsertSize, na.rm = TRUE), # Explicitly call moments::skewness
    Kurtosis_InsertSize = moments::kurtosis(InsertSize, na.rm = TRUE), # Explicitly call moments::kurtosis
    .groups = 'drop'
  )

write.csv(sample_moments,
          file = file.path(fig.dir, "fragment_length_summary_statistics_chr21.csv"),
          row.names = FALSE)
message("Generated fragment_length_summary_statistics_chr21.csv")

dat <- pivot_longer(sample_moments, cols = c("Mean_InsertSize", "Median_InsertSize",
                                             "SD_InsertSize", "Skewness_InsertSize",
                                             "Kurtosis_InsertSize"))

wilcox_skewness <- wilcox.test(Skewness_InsertSize ~ Group, data = sample_moments)
message("\nWilcoxon Test for Skewness_InsertSize:")
print(wilcox_skewness)
wilcox_Kurtosis <- wilcox.test(Kurtosis_InsertSize ~ Group, data = sample_moments) # Corrected variable name
message("\nWilcoxon Test for Kurtosis_InsertSize:")
print(wilcox_Kurtosis)
wilcox_SD_InsertSize <- wilcox.test(SD_InsertSize ~ Group, data = sample_moments)
message("\nWilcoxon Test for SD_InsertSize:")
print(wilcox_SD_InsertSize)

ggplot(dat, aes(x = Group, y = value, fill = Group)) +
  geom_boxplot() + geom_point(position = position_jitter(width = 0.1)) +
  facet_wrap(~ name, ncol = 5, scales = "free_y") + theme_minimal() +
  theme(strip.text = element_text(size = 12, face = 'bold'),
        legend.position = "none")
ggsave(file.path(fig.dir, "FragmentLengthMoments_boxplots.pdf"), width = 12, height = 6)
message("Generated FragmentLengthMoments_boxplots.pdf")

message("--- Short Fragments Test ---")
max_insert_size <- 120
min_insert_size <- 0

all_insert_sizes_list <- list()

for (i in seq_along(all_bam_files)) {
  bam_path <- file.path(data.dir, "bam", all_bam_files[i])
  sample_name <- sample_names[i]
  group <- sample_groups[i]
  
  sizes <- get_insert_sizes(bam_path, min_len = min_insert_size, max_len = max_insert_size)
  
  all_insert_sizes_list[[sample_name]] <- tibble(
    InsertSize = sizes,
    Sample = sample_name,
    Group = group
  )
}
all_insert_sizes_df <- bind_rows(all_insert_sizes_list)

ggplot(all_insert_sizes_df, aes(x = InsertSize, color = Group, fill = Group)) +
  geom_density(alpha = 0.3) +
  labs(title = "Combined Insert Size Distribution of Short Fragments by Group",
       x = "Insert Size (bp)",
       y = "Density") +
  theme_minimal() +
  theme(legend.position = "bottom") +
  xlim(min_insert_size, max_insert_size)
ggsave(file.path(fig.dir, "ShortInsertSizeDensity_byGroup.pdf"), width = 8, height = 6)
message("Generated ShortInsertSizeDensity_byGroup.pdf")


short_fragment_range_min <- 85
short_fragment_range_max <- 105

short_fragment_props <- all_insert_sizes_df %>%
  mutate(IsShortFragment = InsertSize >= short_fragment_range_min & InsertSize <= short_fragment_range_max) %>%
  group_by(Sample, Group) %>%
  summarise(
    ShortFragmentCount = sum(IsShortFragment),
    TotalFragments = n(),
    ShortFragmentProportion = ShortFragmentCount / TotalFragments,
    .groups = 'drop'
  )

message("\nShort fragment proportions:")
print(short_fragment_props)
wilcox_short_fragment_prop <- wilcox.test(ShortFragmentProportion ~ Group, data = short_fragment_props)
message("\nWilcoxon Test for ShortFragmentProportion:")
print(wilcox_short_fragment_prop)
ggplot(short_fragment_props, aes(x = Group, y = ShortFragmentProportion, fill = Group)) +
  geom_boxplot() + geom_point(position = position_jitter(width = 0.1)) +
  ggtitle("Short Fragment Enrichment (85-105bp)") + theme_minimal() +
  theme(legend.position = "none")
ggsave(file.path(fig.dir, "ShortFragmentEnrichment_boxplot.pdf"), width = 6, height = 6)
message("Generated ShortFragmentEnrichment_boxplot.pdf")

message("--- Fragment Coverage and Feature Overlap Analysis ---")
all_fragments_list <- list()
for (i in seq_along(bam_files_ctrl)) {
  sample_name <- ctrl.samples[i]
  bam_path <- file.path(data.dir, "bam", bam_files_ctrl[i])
  fragments_gr <- get_fragments_granges(bam_path)
  mcols(fragments_gr)$Sample <- sample_name
  mcols(fragments_gr)$Group <- "Control"
  all_fragments_list[[sample_name]] <- fragments_gr
}
for (i in seq_along(bam_files_als)) {
  sample_name <- als.sample[i]
  bam_path <- file.path(data.dir, "bam", bam_files_als[i])
  fragments_gr <- get_fragments_granges(bam_path)
  mcols(fragments_gr)$Sample <- sample_name
  mcols(fragments_gr)$Group <- "ALS"
  all_fragments_list[[sample_name]] <- fragments_gr
}

cpg_islands_chr21 <- import(file.path(data.dir, "cpgIslandExt.bed"))
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
promoters_full_genome <- promoters(txdb, upstream = 2000, downstream = 500)
promoters_chr21 <- promoters_full_genome[seqnames(promoters_full_genome) == "chr21"]
promoters_chr21 <- keepSeqlevels(promoters_chr21, "chr21", pruning.mode = "coarse")

all_feature_insert_sizes_df <- list()

for (sample_name in names(all_fragments_list)) {
  fragments_gr <- all_fragments_list[[sample_name]]
  
  cpg_df <- get_feature_insert_sizes(fragments_gr, cpg_islands_chr21, "CpG_Island")
  if (!is.null(cpg_df)) {
    all_feature_insert_sizes_df[[paste0(sample_name, "_CpG")]] <- cpg_df
  }
  
  prom_df <- get_feature_insert_sizes(fragments_gr, promoters_chr21, "Promoter")
  if (!is.null(prom_df)) {
    all_feature_insert_sizes_df[[paste0(sample_name, "_Promoter")]] <- prom_df
  }
}

all_feature_insert_sizes_combined_df <- bind_rows(all_feature_insert_sizes_df)

ggplot(all_feature_insert_sizes_combined_df, aes(x = InsertSize, color = Group, fill = Group)) +
  geom_density(alpha = 0.3) +
  facet_wrap(~ FeatureType, scales = "free_y") +
  labs(title = "Insert Size Distribution by Group within Genomic Features",
       x = "Insert Size (bp)",
       y = "Density") +
  theme_minimal() +
  xlim(0, 400)
ggsave(file.path(fig.dir, "InsertSizeDistribution_Features_chr21.pdf"), width = 10, height = 8)
message("Generated InsertSizeDistribution_Features_chr21.pdf")


median_insert_sizes_per_feature <- all_feature_insert_sizes_combined_df %>%
  group_by(Sample, Group, FeatureType) %>%
  summarise(MedianInsertSize = median(InsertSize, na.rm = TRUE), .groups = 'drop')

message("\nMedian Insert Sizes per Sample within Genomic Features:")
print(median_insert_sizes_per_feature)

feature_types_to_test <- unique(median_insert_sizes_per_feature$FeatureType)
for (ft in feature_types_to_test) {
  message(paste0("\n--- Wilcoxon Test for Median Insert Size in ", ft, " ---"))
  data_subset <- filter(median_insert_sizes_per_feature, FeatureType == ft)
  if (nrow(data_subset) > 1 && length(unique(data_subset$Group)) > 1) {
    wilcox_result <- wilcox.test(MedianInsertSize ~ Group, data = data_subset)
    print(wilcox_result)
  } else {
    message("Not enough data to perform test for this feature type.")
  }
}

ggplot(median_insert_sizes_per_feature, aes(y = MedianInsertSize, x = Group, fill = Group)) +
  geom_boxplot(outlier.shape = NA) + geom_jitter(width = 0.4, na.rm = TRUE, outlier.shape = NA) +
  facet_wrap(~FeatureType) + theme_minimal() +
  labs(title = "Median Insert Size by Group within Genomic Features",
       x = "Group", y = "Median Insert Size (bp)") +
  theme(legend.position = "none")
ggsave(file.path(fig.dir, "MedianInsertSize_Features_boxplot.pdf"), width = 8, height = 6)
message("Generated MedianInsertSize_Features_boxplot.pdf")

message("--- Short/Long Fragment Ratio per Genomic Bin ---")
genome_info_chr21 <- seqinfo(BSgenome.Hsapiens.UCSC.hg38)["chr21"]
bin_size <- 1000000 # 1 Mb bins
genome_bins <- tileGenome(genome_info_chr21, tilewidth = bin_size, cut.last.tile.in.chrom = TRUE)

short_frag_min <- 100
short_frag_max <- 150
long_frag_min <- 151
long_frag_max <- 220

all_bin_ratios_list <- list()
for (sample_name in names(all_fragments_list)) {
  fragments_gr <- all_fragments_list[[sample_name]]
  group_name <- unique(fragments_gr$Group)
  sample_bin_ratios <- calculate_bin_ratios(fragments_gr, genome_bins, sample_name, group_name)
  all_bin_ratios_list[[sample_name]] <- sample_bin_ratios
}
all_bin_ratios_df <- bind_rows(all_bin_ratios_list)

ggplot(all_bin_ratios_df, aes(x = start / 1e6, y = Ratio, color = Group)) +
  geom_line(aes(group = Group)) +
  geom_smooth(method = "loess", se = FALSE) +
  labs(title = paste0("Median Short/Long Fragment Ratio across Chr21 (", bin_size / 1e6, "Mb bins)"),
       x = "Genomic Position (Mb)",
       y = "Short/Long Fragment Ratio") +
  theme_minimal()
ggsave(file.path(fig.dir, "ShortLongRatio_Chr21_Bins.pdf"), width = 10, height = 6)
message("Generated ShortLongRatio_Chr21_Bins.pdf")

wilcox_bin_ratio <- wilcox.test(Ratio ~ Group, data = all_bin_ratios_df)
message("\nWilcoxon Test for Short/Long Fragment Ratio per Genomic Bin:")
print(wilcox_bin_ratio)


message("--- Fragment End Motif Analysis ---")
motif_length <- 6
genome <- BSgenome.Hsapiens.UCSC.hg38
valid_chroms <- "chr21" # Assuming chr21 data for analysis

all_five_prime_motifs <- list()
all_three_prime_motifs <- list()

for (i in seq_along(bam_files_ctrl)) {
  bam_path <- file.path(data.dir, "bam", bam_files_ctrl[i])
  sample_name <- ctrl.samples[i]
  group_label <- "Control"
  
  motifs <- get_fragment_end_motifs(bam_path, ref_genome = genome, motif_len = motif_length, valid_chroms = valid_chroms)
  all_five_prime_motifs[[paste0(sample_name, "_", group_label)]] <- motifs$five_prime_motifs
  all_three_prime_motifs[[paste0(sample_name, "_", group_label)]] <- motifs$three_prime_motifs
}

for (i in seq_along(bam_files_als)) {
  bam_path <- file.path(data.dir, "bam", bam_files_als[i])
  sample_name <- als.sample[i]
  group_label <- "ALS"
  
  motifs <- get_fragment_end_motifs(bam_path, ref_genome = genome, motif_len = motif_length, valid_chroms = valid_chroms)
  all_five_prime_motifs[[paste0(sample_name, "_", group_label)]] <- motifs$five_prime_motifs
  all_three_prime_motifs[[paste0(sample_name, "_", group_label)]] <- motifs$three_prime_motifs
}

# Pool motifs by Group - 5' ends
filtered_5prime_control_list <- filter_empty_DNAStringSets(all_five_prime_motifs[grep("_Control", names(all_five_prime_motifs))])
pooled_5prime_control <- do.call(c, unname(filtered_5prime_control_list))

filtered_5prime_als_list <- filter_empty_DNAStringSets(all_five_prime_motifs[grep("_ALS", names(all_five_prime_motifs))])
pooled_5prime_als <- do.call(c, unname(filtered_5prime_als_list))

# Pool motifs by Group - 3' ends
filtered_3prime_control_list <- filter_empty_DNAStringSets(all_three_prime_motifs[grep("_Control", names(all_three_prime_motifs))])
pooled_3prime_control <- do.call(c, unname(filtered_3prime_control_list))

filtered_3prime_als_list <- filter_empty_DNAStringSets(all_three_prime_motifs[grep("_ALS", names(all_three_prime_motifs))])
pooled_3prime_als <- do.call(c, unname(filtered_3prime_als_list))


# Count k-mers
counts_5prime_control <- oligonucleotideFrequency(pooled_5prime_control, width = motif_length, simplify.as="matrix")
counts_5prime_als <- oligonucleotideFrequency(pooled_5prime_als, width = motif_length, simplify.as="matrix")

df_5prime_control <- data.frame(Motif = colnames(counts_5prime_control), Count_Control = colSums(counts_5prime_control))
df_5prime_als <- data.frame(Motif = colnames(counts_5prime_als), Count_ALS = colSums(counts_5prime_als))

merged_5prime <- full_join(df_5prime_control, df_5prime_als, by = "Motif") %>%
  mutate(Count_Control = replace_na(Count_Control, 0),
         Count_ALS = replace_na(Count_ALS, 0),
         Total_Control = sum(Count_Control),
         Total_ALS = sum(Count_ALS),
         Prop_Control = Count_Control / Total_Control,
         Prop_ALS = Count_ALS / Total_ALS) %>%
  rowwise() %>%
  mutate(Log2FC = log2((Prop_ALS + 1e-9) / (Prop_Control + 1e-9))) %>%
  ungroup() %>%
  arrange(desc(abs(Log2FC)))

message("\nTop differentially abundant 5' end motifs (by absolute Log2FC):")
print(head(merged_5prime, 10))


counts_3prime_control <- oligonucleotideFrequency(pooled_3prime_control, width = motif_length, simplify.as="matrix")
counts_3prime_als <- oligonucleotideFrequency(pooled_3prime_als, width = motif_length, simplify.as="matrix")

df_3prime_control <- data.frame(Motif = colnames(counts_3prime_control), Count_Control = colSums(counts_3prime_control))
df_3prime_als <- data.frame(Motif = colnames(counts_3prime_als), Count_ALS = colSums(counts_3prime_als))

merged_3prime <- full_join(df_3prime_control, df_3prime_als, by = "Motif") %>%
  mutate(Count_Control = replace_na(Count_Control, 0),
         Count_ALS = replace_na(Count_ALS, 0),
         Total_Control = sum(Count_Control),
         Total_ALS = sum(Count_ALS),
         Prop_Control = Count_Control / Total_Control,
         Prop_ALS = Count_ALS / Total_ALS) %>%
  rowwise() %>%
  mutate(Log2FC = log2((Prop_ALS + 1e-9) / (Prop_Control + 1e-9))) %>%
  ungroup() %>%
  arrange(desc(abs(Log2FC)))


# Fisher's Exact Test and FDR correction for 5' motifs
fisher_results_5prime <- merged_5prime %>%
  rowwise() %>%
  mutate(
    p_value = tryCatch({
      test_table <- matrix(c(Count_Control, Total_Control - Count_Control,
                             Count_ALS, Total_ALS - Count_ALS), nrow = 2, byrow = TRUE)
      fisher.test(test_table)$p.value
    }, error = function(e) NA)
  ) %>%
  ungroup() %>%
  mutate(p_adj = p.adjust(p_value, method = "BH")) %>%
  arrange(p_adj)

message("\nTop 5' end motifs by adjusted p-value:")
print(head(fisher_results_5prime, 10))
write.csv(fisher_results_5prime,
          file = file.path(fig.dir, "Motif_5primeDE.csv"), row.names = FALSE)
message("Generated Motif_5primeDE.csv")

# Fisher's Exact Test and FDR correction for 3' motifs
fisher_results_3prime <- merged_3prime %>%
  rowwise() %>%
  mutate(
    p_value = tryCatch({
      test_table <- matrix(c(Count_Control, Total_Control - Count_Control,
                             Count_ALS, Total_ALS - Count_ALS), nrow = 2, byrow = TRUE)
      fisher.test(test_table)$p.value
    }, error = function(e) NA)
  ) %>%
  ungroup() %>%
  mutate(p_adj = p.adjust(p_value, method = "BH")) %>%
  arrange(p_adj)

message("\nTop 3' end motifs by adjusted p-value:")
print(head(fisher_results_3prime, 10))
write.csv(fisher_results_3prime,
          file = file.path(fig.dir, "Motif_3primeDE.csv"), row.names = FALSE)
message("Generated Motif_3primeDE.csv")

# Plotting Sequence Logos
top_als_5prime_motifs <- fisher_results_5prime %>%
  filter(Log2FC > 0) %>%
  arrange(desc(Log2FC)) %>%
  head(5) %>% pull(Motif)

if (length(top_als_5prime_motifs) > 0) {
  p_5prime_als <- ggplot() + geom_logo(top_als_5prime_motifs) +
    labs(title = paste0("Top 5' End Motifs Enriched in ALS (", motif_length, "bp)")) + theme_logo()
  ggsave(file.path(fig.dir, paste0("top_5prime_als_motifs_", motif_length, "bp.png")), plot = p_5prime_als, width = 8, height = 4, dpi = 300)
  message(paste0("Generated top_5prime_als_motifs_", motif_length, "bp.png"))
} else {
  message("No 5' motifs enriched in ALS found for plotting.")
}

top_control_5prime_motifs <- fisher_results_5prime %>%
  filter(Log2FC < 0) %>%
  arrange(Log2FC) %>%
  head(5) %>% pull(Motif)

if (length(top_control_5prime_motifs) > 0) {
  p_5prime_control <- ggplot() + geom_logo(top_control_5prime_motifs) +
    labs(title = paste0("Top 5' End Motifs Enriched in Control (", motif_length, "bp)")) + theme_logo()
  ggsave(file.path(fig.dir, paste0("top_5prime_control_motifs_", motif_length, "bp.png")), plot = p_5prime_control, width = 8, height = 4, dpi = 300)
  message(paste0("Generated top_5prime_control_motifs_", motif_length, "bp.png"))
} else {
  message("No 5' motifs enriched in Control found for plotting.")
}

top_als_3prime_motifs <- fisher_results_3prime %>%
  filter(Log2FC > 0) %>%
  arrange(desc(Log2FC)) %>%
  head(5) %>% pull(Motif)

if (length(top_als_3prime_motifs) > 0) {
  p_3prime_als <- ggplot() + geom_logo(top_als_3prime_motifs) +
    labs(title = paste0("Top 3' End Motifs Enriched in ALS (", motif_length, "bp)")) + theme_logo()
  ggsave(file.path(fig.dir, paste0("top_3prime_als_motifs_", motif_length, "bp.png")), plot = p_3prime_als, width = 8, height = 4, dpi = 300)
  message(paste0("Generated top_3prime_als_motifs_", motif_length, "bp.png"))
} else {
  message("No 3' motifs enriched in ALS found for plotting.")
}

top_control_3prime_motifs <- fisher_results_3prime %>%
  filter(Log2FC < 0 & p_adj < 0.05) %>%
  arrange(Log2FC) %>%
  head(5) %>% pull(Motif)

if (length(top_control_3prime_motifs) > 0) {
  p_3prime_control <- ggplot() + geom_logo(top_control_3prime_motifs) +
    labs(title = paste0("Top 3' End Motifs Enriched in Control (", motif_length, "bp)")) + theme_logo()
  ggsave(file.path(fig.dir, paste0("top_3prime_control_motifs_", motif_length, "bp.png")), plot = p_3prime_control, width = 8, height = 4, dpi = 300)
  message(paste0("Generated top_3prime_control_motifs_", motif_length, "bp.png"))
} else {
  message("No 3' motifs enriched in Control found for plotting.")
}

message("--- Fragment Start/End Position Analysis around Features ---")
flanking_window <- 1000
bin_size_bp <- 5
tss_gr <- transcripts(txdb)
tss_gr <- resize(tss_gr, width = 1, fix = "start")
tss_gr <- tss_gr[seqnames(tss_gr) == valid_chroms]
seqlevelsStyle(tss_gr) <- "UCSC"
tss_gr <- keepSeqlevels(tss_gr, valid_chroms, pruning.mode = "coarse")


promoter_centers <- promoters_chr21
start(promoter_centers) <- start(promoters_chr21) + floor(width(promoters_chr21) / 2)
end(promoter_centers) <- start(promoter_centers)
seqinfo(promoter_centers) <- seqinfo(genome)[seqlevels(promoter_centers)]

cpg_centers <- cpg_islands_chr21
start(cpg_centers) <- start(cpg_islands_chr21) + floor(width(cpg_islands_chr21) / 2)
end(cpg_centers) <- start(cpg_centers)

all_end_pos_dfs <- list()

for (sample_name in names(all_fragments_list)) {
  fragments_gr_sample <- all_fragments_list[[sample_name]]
  
  prom_end_pos_df <- calculate_relative_distances(
    fragments_gr_sample, promoter_centers, "Promoter_TSS", window_size = flanking_window
  )
  if (!is.null(prom_end_pos_df)) {
    all_end_pos_dfs[[paste0(sample_name, "_Promoter")]] <- prom_end_pos_df
  }
  
  cpg_end_pos_df <- calculate_relative_distances(
    fragments_gr_sample, cpg_centers, "CpG_Island_Center", window_size = flanking_window
  )
  if (!is.null(cpg_end_pos_df)) {
    all_end_pos_dfs[[paste0(sample_name, "_CpG")]] <- cpg_end_pos_df
  }
  
}

final_end_pos_df <- bind_rows(all_end_pos_dfs)

# Plot Fragment End Distributions
if (nrow(final_end_pos_df) > 0) {
  ggplot(final_end_pos_df, aes(x = RelativePosition, color = Group, fill = Group)) +
    geom_histogram(aes(y = after_stat(density)), binwidth = bin_size_bp, alpha = 0.3, position = "identity") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
    facet_grid(EndType ~ FeatureType, scales = "free_y") +
    labs(title = paste0("Fragment End Distribution around Genomic Features (", flanking_window, "bp window)"),
         x = paste0("Distance from Feature Center (bp)"),
         y = "Density") +
    theme_minimal() +
    theme(legend.position = "bottom") +
    xlim(-flanking_window, flanking_window)
  ggsave(file.path(fig.dir, "fragment_end_distributions_chr21.png"), width = 10, height = 8, dpi = 300)
  message("Generated fragment_end_distributions_chr21.png")
} else {
  message("No fragment ends found overlapping specified features for this analysis to plot initial distributions.")
}

# Stats for Fragment End Distributions
breaks <- seq(-flanking_window, flanking_window, by = bin_size_bp)
labels <- seq(-flanking_window + bin_size_bp/2, flanking_window - bin_size_bp/2, by = bin_size_bp)

binned_counts <- final_end_pos_df %>%
  mutate(Bin = cut(RelativePosition,
                   breaks = breaks,
                   labels = labels,
                   include.lowest = TRUE,
                   right = TRUE)) %>%
  filter(!is.na(Bin)) %>%
  group_by(Group, FeatureType, EndType, Bin) %>%
  summarise(Count = n(), .groups = 'drop_last') %>%
  ungroup()

all_combinations <- expand_grid(
  Group = unique(binned_counts$Group),
  FeatureType = unique(binned_counts$FeatureType),
  EndType = unique(binned_counts$EndType),
  Bin = unique(binned_counts$Bin)
)

binned_counts_full <- left_join(all_combinations, binned_counts, by = c("Group", "FeatureType", "EndType", "Bin")) %>%
  mutate(Count = replace_na(Count, 0))

total_counts_per_group_end_feature <- binned_counts_full %>%
  group_by(Group, FeatureType, EndType) %>%
  summarise(TotalCount = sum(Count), .groups = 'drop')

binned_proportions <- left_join(binned_counts_full, total_counts_per_group_end_feature, by = c("Group", "FeatureType", "EndType")) %>%
  mutate(Proportion = Count / TotalCount)

fisher_results_bins <- data.frame(binned_proportions) %>%
  dplyr::select(FeatureType, EndType, Bin, Group, Count, TotalCount) %>%
  pivot_wider(names_from = Group, values_from = c(Count, TotalCount)) %>%
  rowwise() %>%
  mutate(
    Count_Control = replace_na(Count_Control, 0),
    Count_ALS = replace_na(Count_ALS, 0),
    TotalCount_Control = replace_na(TotalCount_Control, 0),
    TotalCount_ALS = replace_na(TotalCount_ALS, 0),
    p_value = tryCatch({
      if (TotalCount_Control == 0 && TotalCount_ALS == 0) {
        return(NA_real_)
      }
      test_table <- matrix(c(Count_Control, TotalCount_Control - Count_Control,
                             Count_ALS, TotalCount_ALS - Count_ALS),
                           nrow = 2, byrow = TRUE)
      fisher.test(test_table)$p.value
    }, error = function(e) NA_real_)
  ) %>%
  ungroup() %>%
  mutate(p_adj = p.adjust(p_value, method = "BH"))


message("\nTop bins with statistically significant differences (p_adj < 0.05):")
significant_bins <- fisher_results_bins %>% filter(p_value < 0.05)
print(significant_bins)

if (nrow(final_end_pos_df) > 0 && nrow(significant_bins) > 0) {
  highlight_rects <- significant_bins %>%
    mutate(
      xmin = as.numeric(as.character(Bin)) - bin_size_bp / 2,
      xmax = as.numeric(as.character(Bin)) + bin_size_bp / 2,
      ymin = -Inf,
      ymax = Inf
    ) %>%
    dplyr::select(FeatureType, EndType, xmin, xmax, ymin, ymax) %>%
    distinct()
  
  p <- ggplot(final_end_pos_df, aes(x = RelativePosition, color = Group, fill = Group)) +
    geom_rect(data = highlight_rects,
              aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
              fill = "red",
              alpha = 0.5,
              inherit.aes = FALSE) +
    geom_histogram(aes(y = after_stat(density)),
                   binwidth = bin_size_bp,
                   alpha = 0.3,
                   position = "identity") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
    facet_grid(EndType ~ FeatureType, scales = "free_y") +
    labs(title = paste0("Fragment End Distribution around Genomic Features (", flanking_window, "bp window)"),
         x = "Distance from Feature Center (bp)",
         y = "Density") +
    theme_minimal() +
    theme(legend.position = "bottom",
          plot.title = element_text(hjust = 0.5)) +
    xlim(-flanking_window, flanking_window)
  
  ggsave(file.path(fig.dir, "fragment_end_distributions_significant_bins_chr21.png"), plot = p, width = 12, height = 8, dpi = 300)
  message("Generated fragment_end_distributions_significant_bins_chr21.png")
} else if (nrow(final_end_pos_df) > 0) {
  message("No significant bins found to highlight for fragment end distributions. Plotting original density profiles.")
  p <- ggplot(final_end_pos_df, aes(x = RelativePosition, color = Group, fill = Group)) +
    geom_histogram(aes(y = after_stat(density)), binwidth = bin_size_bp, alpha = 0.3, position = "identity") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
    facet_grid(EndType ~ FeatureType, scales = "free_y") +
    labs(title = paste0("Fragment End Distribution around Genomic Features (", flanking_window, "bp window)"),
         x = "Distance from Feature Center (bp)",
         y = "Density") +
    theme_minimal() +
    theme(legend.position = "bottom",
          plot.title = element_text(hjust = 0.5)) +
    xlim(-flanking_window, flanking_window)
  ggsave(file.path(fig.dir, "fragment_end_distributions_chr21_no_highlights.png"), plot = p, width = 10, height = 8, dpi = 300)
  message("Generated fragment_end_distributions_chr21_no_highlights.png")
} else {
  message("No data available in final_end_pos_df to plot fragment end distributions.")
}


message("--- Classification Analysis ---")
feature_data <- sample_moments %>%
  dplyr::select(Sample, Group, Mean_InsertSize, Median_InsertSize,
                SD_InsertSize, Skewness_InsertSize, Kurtosis_InsertSize)

feature_data <- feature_data %>%
  left_join(mono_di_proportions %>% dplyr::select(Sample, Mono_Proportion, Mono_Di_Ratio), by = "Sample")

feature_data <- feature_data %>%
  left_join(short_fragment_props %>% dplyr::select(Sample, ShortFragmentCount, TotalFragments, ShortFragmentProportion), by = "Sample")

for (col_name in names(feature_data)) {
  if (is.numeric(feature_data[[col_name]]) && any(is.na(feature_data[[col_name]]))) {
    feature_data[[col_name]][is.na(feature_data[[col_name]])] <- mean(feature_data[[col_name]], na.rm = TRUE)
  }
}

feature_data$Group <- factor(feature_data$Group)
classification_df <- feature_data %>% dplyr::select(-Sample)

message("Prepared classification data (first 5 rows):")
print(head(classification_df, 5))
message(paste0("Total samples for classification: ", nrow(classification_df)))
message(paste0("Number of features: ", ncol(classification_df) - 1))

train_control <- trainControl(
  method = "repeatedcv",
  number = 10,
  repeats = 50,
  classProbs = TRUE,
  summaryFunction = twoClassSummary,
  savePredictions = "final"
)

message("--- Training Logistic Regression Model ---")
model_glm <- train(
  Group ~ .,
  data = classification_df,
  method = "glm",
  family = "binomial",
  trControl = train_control,
  metric = "ROC"
)
message("Logistic Regression Model Training Complete.")
print(model_glm)

message("--- Training Random Forest Model ---")
model_rf <- train(
  Group ~ .,
  data = classification_df,
  method = "rf",
  trControl = train_control,
  metric = "ROC",
  tuneLength = 4
)
message("Random Forest Model Training Complete.")
print(model_rf)

# Evaluate and Present Metrics
message("\nLogistic Regression Performance:")
print(model_glm$results)
message(paste0("Best Logistic Regression ROC: ", round(max(model_glm$results$ROC), 3)))
message(paste0("Best Logistic Regression Sensitivity: ", round(model_glm$results$Sens[which.max(model_glm$results$ROC)], 3)))
message(paste0("Best Logistic Regression Specificity: ", round(model_glm$results$Spec[which.max(model_glm$results$ROC)], 3)))

message("\nRandom Forest Performance:")
print(model_rf$results)
message(paste0("Best Random Forest ROC: ", round(max(model_rf$results$ROC), 3)))
message(paste0("Best Random Forest Sensitivity: ", round(model_rf$results$Sens[which.max(model_rf$results$ROC)], 3)))
message(paste0("Best Random Forest Specificity: ", round(model_rf$results$Spec[which.max(model_rf$results$ROC)], 3)))


if ("mtry" %in% names(model_glm$bestTune)) { # This is not typically present for GLM but keeping for robustness
  glm_pred <- model_glm$pred %>% filter(mtry == model_glm$bestTune$mtry)
} else {
  glm_pred <- model_glm$pred # No mtry for glm, use all predictions
}

message("\nLogistic Regression Confusion Matrix:")
cm_glm <- confusionMatrix(data = glm_pred$pred, reference = glm_pred$obs, positive = "ALS")
print(cm_glm)
message(paste0("Logistic Regression Precision: ", round(cm_glm$byClass["Pos Pred Value"], 3)))
message(paste0("Logistic Regression F1-Score: ", round(cm_glm$byClass["F1"], 3)))

rf_pred <- model_rf$pred %>% filter(mtry == model_rf$bestTune$mtry)

message("\nRandom Forest Confusion Matrix:")
cm_rf <- confusionMatrix(data = rf_pred$pred, reference = rf_pred$obs, positive = "ALS")
print(cm_rf)
message(paste0("Random Forest Precision: ", round(cm_rf$byClass["Pos Pred Value"], 3)))
message(paste0("Random Forest F1-Score: ", round(cm_rf$byClass["F1"], 3)))


# ROC Curve and AUC
roc_glm <- roc(response = glm_pred$obs, predictor = glm_pred$ALS, levels = c("Control", "ALS"))
message(paste0("Logistic Regression AUC: ", round(auc(roc_glm), 3)))
pdf(file.path(fig.dir, "roc_curve_glm.pdf"), width = 6, height = 6) # Changed to PDF for better quality
plot(roc_glm, main = "ROC Curve - Logistic Regression", print.auc = TRUE, col = "blue")
dev.off()
message("Generated roc_curve_glm.pdf")

roc_rf <- roc(response = rf_pred$obs, predictor = rf_pred$ALS, levels = c("Control", "ALS"))
message(paste0("Random Forest AUC: ", round(auc(roc_rf), 3)))
pdf(file.path(fig.dir, "roc_curve_rf.pdf"), width = 6, height = 6) # Changed to PDF
plot(roc_rf, main = "ROC Curve - Random Forest", print.auc = TRUE, col = "red")
dev.off()
message("Generated roc_curve_rf.pdf")

pdf(file.path(fig.dir, "roc_curves_combined.pdf"), width = 7, height = 6) # Changed to PDF
plot(roc_glm, col = "blue", main = "ROC Curves: ALS vs Control")
lines(roc_rf, col = "red")
legend("bottomright", legend = c(paste0("Logistic Regression (AUC = ", round(auc(roc_glm), 2), ")"),
                                 paste0("Random Forest (AUC = ", round(auc(roc_rf), 2), ")")),
       col = c("blue", "red"), lty = 1, cex = 0.8)
dev.off()
message("Generated roc_curves_combined.pdf")


# Feature Importance (for Random Forest)
if (!is.null(varImp(model_rf)$importance)) {
  rf_importance <- varImp(model_rf)$importance
  plot_importance <- ggplot(rf_importance, aes(x = reorder(rownames(rf_importance), Overall), y = Overall)) +
    geom_bar(stat = "identity", fill = "steelblue") +
    coord_flip() +
    labs(title = "Random Forest Feature Importance",
         x = "Feature",
         y = "Importance (Overall)") +
    theme_minimal()
  ggsave(file.path(fig.dir, "rf_feature_importance.png"), plot = plot_importance, width = 8, height = 6, dpi = 300)
  message("Generated rf_feature_importance.png")
} else {
  message("Feature importance not available or could not be computed for Random Forest model.")
}

# Save metrics
glm_overall_metrics <- as.data.frame(t(cm_glm$overall))
glm_overall_metrics <- glm_overall_metrics %>%
  dplyr::select(Accuracy, Kappa) %>%
  rename(Overall_Accuracy = Accuracy, Overall_Kappa = Kappa)

glm_byclass_metrics <- as.data.frame(t(cm_glm$byClass))
glm_byclass_metrics <- glm_byclass_metrics %>%
  dplyr::select(`Sensitivity`, `Specificity`, `Pos Pred Value`, `Neg Pred Value`, `F1`) %>%
  rename(Sensitivity = `Sensitivity`, Specificity = `Specificity`,
         Precision = `Pos Pred Value`, Negative_Predictive_Value = `Neg Pred Value`, F1_Score = `F1`)

glm_auc <- data.frame(AUC = as.numeric(auc(roc_glm)))

glm_metrics_df <- bind_cols(glm_overall_metrics, glm_byclass_metrics, glm_auc) %>%
  mutate(Model = "Logistic Regression") %>%
  dplyr::select(Model, everything())


rf_overall_metrics <- as.data.frame(t(cm_rf$overall))
rf_overall_metrics <- rf_overall_metrics %>%
  dplyr::select(Accuracy, Kappa) %>%
  rename(Overall_Accuracy = Accuracy, Overall_Kappa = Kappa)

rf_byclass_metrics <- as.data.frame(t(cm_rf$byClass))
rf_byclass_metrics <- rf_byclass_metrics %>%
  dplyr::select(`Sensitivity`, `Specificity`, `Pos Pred Value`, `Neg Pred Value`, `F1`) %>%
  rename(Sensitivity = `Sensitivity`, Specificity = `Specificity`,
         Precision = `Pos Pred Value`, Negative_Predictive_Value = `Neg Pred Value`, F1_Score = `F1`)

rf_auc <- data.frame(AUC = as.numeric(auc(roc_rf)))

rf_metrics_df <- bind_cols(rf_overall_metrics, rf_byclass_metrics, rf_auc) %>%
  mutate(Model = "Random Forest") %>%
  dplyr::select(Model, everything())

all_classification_metrics <- bind_rows(glm_metrics_df, rf_metrics_df)

output_metrics_csv_file <- "classification_performance_metrics.csv"
write.csv(all_classification_metrics,
          file = file.path(fig.dir, output_metrics_csv_file), row.names = FALSE)
message(paste0("Generated ", output_metrics_csv_file))

message("cfDNA analysis complete. Check the 'Figures' directory for plots and CSVs.")