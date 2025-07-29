# install neccesary Packages within R if not already installed
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

#BiocManager::install("Rsamtools")
#install.packages("ggplot2")
#install.packages("dplyr")

library(Rsamtools)
library(ggplot2)
library(dplyr)
library(tidyr)
source("functions.R")

# Working directory should be Analysis folder
data.dir <- "../celfie"
fig.dir <- "../Figures"
#Read Metadata------
meta <- read.csv(file.path(data.dir,"celfie_cfDNA_ss.csv"))
ctrl.samples <- meta$Run[which(meta$disease_status=="ctrl")]
als.sample <- meta$Run[which(meta$disease_status=="als")]
bam_files_ctrl <- c()
for(name in ctrl.samples){
  bam_files_ctrl <- c(bam_files_ctrl,list.files(path = file.path(data.dir,"bam"),
                                                pattern = paste0("^",name,".*\\.bam$"), recursive = T))
}
bam_files_als <- c()
for(name in als.sample){
  bam_files_als <- c(bam_files_als,list.files(path = file.path(data.dir,"bam"),
                                              pattern = paste0("^",name,".*\\.bam$"), recursive = T))
}


# Combine for easier iteration and plotting
all_bam_files <- c(bam_files_ctrl, bam_files_als)
sample_groups <- c(rep("Control", length(bam_files_ctrl)),
                   rep("ALS", length(bam_files_als)))
sample_names <-c(ctrl.samples,als.sample)

# Insert Analysis -------------
# initial Explorations
max_insert_size <- 400 
min_insert_size <- 0 # Often start at 0 or a small positive value

# --- Extract insert sizes for all samples ---

all_insert_sizes_list <- list()

for (i in seq_along(all_bam_files)) {
  bam_path <- file.path(data.dir,"bam",all_bam_files[i])
  sample_name <- sample_names[i]
  group <- sample_groups[i]
  
  # Get insert sizes for the current BAM file
  sizes <- get_insert_sizes(bam_path, min_len = min_insert_size, max_len = max_insert_size)
  
  # Store in the list with sample name and group
  all_insert_sizes_list[[sample_name]] <- tibble(
    InsertSize = sizes,
    Sample = sample_name,
    Group = group
  )
}

# Combine all data frames into one large data frame for ggplot
all_insert_sizes_df <- bind_rows(all_insert_sizes_list)

# --- Generate Frequency Plot (Histogram) ---

# Plotting by sample
pdf(file.path(fig.dir,"InsertSizeDensity_bySample.pdf"),width = 10,height=10)
ggplot(all_insert_sizes_df, aes(x = InsertSize, fill = Group)) +
  geom_histogram(binwidth = 1, alpha = 0.7, position = "identity") +
  facet_wrap(~ Sample, scales = "free_y", ncol = 2) + 
  labs(title = "Insert Size Distribution per Sample",
       x = "Insert Size (bp)",
       y = "Frequency") +
  theme_minimal() +
  theme(legend.position = "bottom")
dev.off()
# Plotting by group 
pdf(file.path(fig.dir,"InsertSizeDensity_byGroup_300max.pdf"))
ggplot(all_insert_sizes_df, aes(x = InsertSize, color = Group, fill = Group)) +
  geom_density(alpha = 0.3) + 
  labs(title = "Combined Insert Size Distribution by Group",
       x = "Insert Size (bp)",
       y = "Density") +
  theme_minimal() +
  theme(legend.position = "bottom") +
  xlim(min_insert_size, max_insert_size) 
dev.off()


# differecns look apparent especially around  mono-nucleosomal and di-nucleosomal fragments

# lets test around those regions
mono_min <- 150
mono_max <- 180 # Or 190, 200 depending on where your first peak truly ends
di_min <- 300
di_max <- 380 # Or 380, 400 where your second peak truly ends

df_fragment_types <- all_insert_sizes_df %>%
  mutate(
    FragmentType = case_when(
      InsertSize >= mono_min & InsertSize <= mono_max ~ "Mono-nucleosomal",
      InsertSize >= di_min & InsertSize <= di_max ~ "Di-nucleosomal",
      TRUE ~ "Other_Fragment" # Catch all other sizes
    )
  )

fragment_counts_proportions <- df_fragment_types %>%
  group_by(Sample, Group, FragmentType) %>%
  summarise(Count = n(), .groups = 'drop_last') %>% # Count fragments for each type
  mutate(Proportion = Count / sum(Count)) %>%       # Calculate proportion within each sample
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
  # Handle potential Inf values if Di_Count is 0
  mutate(Mono_Di_Ratio = ifelse(is.infinite(Mono_Di_Ratio), NA, Mono_Di_Ratio))

# Test for Mono_Proportion difference
wilcox_mono_prop <- wilcox.test(Mono_Proportion ~ Group, data = mono_di_proportions)
print("\nWilcoxon Rank-Sum Test for Mono-nucleosomal Proportion difference:")
print(wilcox_mono_prop)
# Box plot for Mono_Proportion
pdf(file.path(fig.dir,"MononucleosomalProportion_boxplot.pdf"),width = 10,height=10)
ggplot(mono_di_proportions, aes(x = Group, y = Mono_Proportion, fill = Group)) +
  geom_boxplot(alpha = 0.7) +
  geom_point(position = position_jitter(width = 0.1), alpha = 0.8, color = "black") + 
  labs(title = "Mono-nucleosomal Fragment Proportion by Group",
       x = "Group",
       y = "Proportion of Fragments (150-180bp)") +
  theme_minimal() +
  theme(legend.position = "none")
dev.off()

#outlier?
mono_di_proportions_flagged <- identify_outliers_iqr(
  data = mono_di_proportions,
  variable_col = Mono_Proportion,
  group_col = Group,
  k = 1.5
)

# --- Plotting with Outliers Highlighted (Mono_Proportion) ---

ggplot(mono_di_proportions_flagged, aes(x = Group, y = Mono_Proportion, fill = Group)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) + 
  geom_point(aes(color = is_outlier), 
             position = position_jitter(width = 0.1), 
             size = 2, # Adjust point size
             alpha = 0.8) +
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "black"), 
                     labels = c("TRUE" = "Outlier", "FALSE" = "Not Outlier")) + 
  labs(title = "Mono-nucleosomal Fragment Proportion with Outliers Highlighted",
       x = "Group",
       y = "Proportion of Fragments (140-180bp)",
       color = "Data Point Type") + # Legend title for the 'color' aesthetic
  theme_minimal() +
  theme(legend.position = "bottom", # Place legend at the bottom
        plot.title = element_text(hjust = 0.5)) # Center the plot title
wilcox_mono_prop_cleaned <- wilcox.test(Mono_Proportion ~ Group, data = mono_di_proportions_cleaned)
print(wilcox_mono_prop_cleaned)

# Test for Mono_Di_Ratio difference 
wilcox_mono_di_ratio <- wilcox.test(Mono_Di_Ratio ~ Group, data = mono_di_proportions_cleaned)
print(wilcox_mono_di_ratio)

#Based on this Wilcoxon Rank-Sum test, there is no statistically significant difference

##Plot Nucleosome Box plots --------

mono_di_proportions_cleaned <- remove_outliers_iqr(
  data = mono_di_proportions,
  variable_col = Mono_Proportion,
  group_col = Group,
  k = 1.5 
)
pdf(file.path(fig.dir,"MononucleosomalProportionCLEAN_boxplot.pdf"),width = 10,height=10)
ggplot(mono_di_proportions_cleaned, aes(x = Group, y = Mono_Proportion, fill = Group)) +
  geom_boxplot(alpha = 0.7) +
  geom_point(position = position_jitter(width = 0.1), alpha = 0.8, color = "black") + # Show individual data points
  labs(title = "Mono-nucleosomal Fragment Proportion by Group",
       x = "Group",
       y = "Proportion of Fragments (150-180bp)") +
  theme_minimal() +
  theme(legend.position = "none")
dev.off()

# Box plot for Mono/Di-nucleosomal Ratio 
pdf(file.path(fig.dir,"Mono_Di_Ratio_boxplot.pdf"),width = 10,height=10)
ggplot(mono_di_proportions_cleaned, aes(x = Group, y = Mono_Di_Ratio, fill = Group)) +
  geom_boxplot(alpha = 0.7) +
  geom_point(position = position_jitter(width = 0.1), alpha = 0.8, color = "black") +
  labs(title = "Mono/Di-nucleosomal Fragment Ratio by Group",
       x = "Group",
       y = "Ratio of Fragments (Mono- / Di-nucleosomal)") +
  theme_minimal() +
  theme(legend.position = "none")
dev.off()


## Distribution tests -------
library(moments) # to extract skewness and kurtosis

sample_moments <- all_insert_sizes_df %>%
  group_by(Sample, Group) %>%
  summarise(
    Mean_InsertSize = mean(InsertSize, na.rm = TRUE),
    Median_InsertSize = median(InsertSize, na.rm = TRUE),
    SD_InsertSize = sd(InsertSize, na.rm = TRUE),
    Skewness_InsertSize = skewness(InsertSize, na.rm = TRUE),
    Kurtosis_InsertSize = kurtosis(InsertSize, na.rm = TRUE),
    .groups = 'drop'
  )

print(sample_moments)

dat <- pivot_longer(sample_moments,cols = c("Mean_InsertSize","Median_InsertSize",
                                            "SD_InsertSize","Skewness_InsertSize",
                                            "Kurtosis_InsertSize"))
# Wilcoxon tests
wilcox_skewness <- wilcox.test(Skewness_InsertSize ~ Group, data = sample_moments)
print(wilcox_skewness)
wilcox_Kurtosis <- wilcox.test(SD_InsertSize ~ Group, data = sample_moments)
print(wilcox_Kurtosis)


# Plotting (e.g., boxplots for each moment)
ggplot(dat, aes(x = Group, y = value, fill = Group)) +
  geom_boxplot() + geom_point(position = position_jitter(width = 0.1)) +
  facet_wrap(~ name,ncol=5,scales = "free_y")+ theme_minimal()+
  theme(strip.text = element_text(size=12, face='bold'),
        legend.position = "none")


## Short Fragments Test ---------
max_insert_size <- 120 
min_insert_size <- 0 

all_insert_sizes_list <- list()

for (i in seq_along(all_bam_files)) {
  bam_path <- file.path(data.dir,"bam",all_bam_files[i])
  sample_name <- sample_names[i]
  group <- sample_groups[i]
  
  # Get insert sizes for the current BAM file
  sizes <- get_insert_sizes(bam_path, min_len = min_insert_size, max_len = max_insert_size)
  
  # Store in the list with sample name and group
  all_insert_sizes_list[[sample_name]] <- tibble(
    InsertSize = sizes,
    Sample = sample_name,
    Group = group
  )
}

# Combine all data frames into one large data frame for ggplot
all_insert_sizes_df <- bind_rows(all_insert_sizes_list)


## Plot distribution of short frags
ggplot(all_insert_sizes_df, aes(x = InsertSize, color = Group, fill = Group)) +
  geom_density(alpha = 0.3) + 
  labs(title = "Combined Insert Size Distribution by Group",
       x = "Insert Size (bp)",
       y = "Density") +
  theme_minimal() +
  theme(legend.position = "bottom") +
  xlim(min_insert_size, max_insert_size) 

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

print(short_fragment_props)
wilcox_short_fragment_prop <- wilcox.test(ShortFragmentProportion ~ Group, data = short_fragment_props)
print(wilcox_short_fragment_prop)
ggplot(short_fragment_props,aes(x=Group,y=ShortFragmentProportion,fill=Group)) +
  geom_boxplot() + geom_point(position = position_jitter(width = 0.1))+
  ggtitle("ShortFragment Enrichment ")

## Fragment Coverage ------------
library(GenomicAlignments)
library(GenomicRanges)
library(rtracklayer)
library(GenomicFeatures)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)


all_fragments_list <- list()
for (i in seq_along(bam_files_ctrl)) {
  sample_name <- ctrl.samples[i]
  bam_path <- file.path(data.dir,"bam",bam_files_ctrl[i])
  all_fragments_list[[sample_name]] <- get_fragments_granges(bam_path)
  mcols(all_fragments_list[[sample_name]])$Sample <- sample_name
  mcols(all_fragments_list[[sample_name]])$Group <- "Control"
}
for (i in seq_along(bam_files_als)) {
  sample_name <- als.sample[i]
  bam_path <- file.path(data.dir,"bam",bam_files_als[i])
  all_fragments_list[[sample_name]] <- get_fragments_granges(bam_path)
  mcols(all_fragments_list[[sample_name]])$Sample <- sample_name
  mcols(all_fragments_list[[sample_name]])$Group <- "ALS"
}
###Overall Insert Size Distribution ----------
cpg_islands_chr21 <- import(file.path(data.dir,"cpgIslandExt.bed"))
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
promoters_full_genome <- promoters(txdb, upstream = 2000, downstream = 500)
enhancers_full <- 
  seqlevelsStyle(promoters_full_genome) <- "UCSC"
promoters_chr21 <- promoters_full_genome[seqnames(promoters_full_genome) == "chr21"]
promoters_chr21 <- keepSeqlevels(promoters_chr21, "chr21", pruning.mode="coarse")

all_feature_insert_sizes_df <- list()

for (sample_name in names(all_fragments_list)) {
  fragments_gr <- all_fragments_list[[sample_name]]
  
  # Analyze for CpG islands
  cpg_df <- get_feature_insert_sizes(fragments_gr, cpg_islands_chr21, "CpG_Island")
  if (!is.null(cpg_df)) {
    all_feature_insert_sizes_df[[paste0(sample_name, "_CpG")]] <- cpg_df
  }
  
  # Analyze for Promoters
  prom_df <- get_feature_insert_sizes(fragments_gr, promoters_chr21, "Promoter")
  if (!is.null(prom_df)) {
    all_feature_insert_sizes_df[[paste0(sample_name, "_Promoter")]] <- prom_df
  }
}

# Combine results
all_feature_insert_sizes_combined_df <- bind_rows(all_feature_insert_sizes_df)

# --- Plotting (Density plots per feature type, by group) ---
ggplot(all_feature_insert_sizes_combined_df, aes(x = InsertSize, color = Group, fill = Group)) +
  geom_density(alpha = 0.3) +
  facet_wrap(~ FeatureType, scales = "free_y") + 
  labs(title = "Insert Size Distribution by Group within Genomic Features",
       x = "Insert Size (bp)",
       y = "Density") +
  theme_minimal() +
  xlim(0, 400) # Adjust as needed

# Calculate median insert size per sample for each feature type
median_insert_sizes_per_feature <- all_feature_insert_sizes_combined_df %>%
  group_by(Sample, Group, FeatureType) %>%
  summarise(MedianInsertSize = median(InsertSize, na.rm = TRUE), .groups = 'drop')

print("Median Insert Sizes per Sample within Genomic Features:")
print(median_insert_sizes_per_feature)

# Perform Wilcoxon test for each feature type
feature_types_to_test <- unique(median_insert_sizes_per_feature$FeatureType)
for (ft in feature_types_to_test) {
  cat(paste0("\n--- Wilcoxon Test for Median Insert Size in ", ft, " ---\n"))
  data_subset <- filter(median_insert_sizes_per_feature, FeatureType == ft)
  if (nrow(data_subset) > 1 && length(unique(data_subset$Group)) > 1) {
    wilcox_result <- wilcox.test(MedianInsertSize ~ Group, data = data_subset)
    print(wilcox_result)
  } else {
    print("Not enough data to perform test for this feature type.")
  }
}


ggplot(median_insert_sizes_per_feature,aes(y=MedianInsertSize,x=Group,fill = Group)) +
  geom_boxplot(outlier.shape = NA)+geom_jitter(width = 0.4,na.rm = T,outlier.shape = NA) +
  facet_wrap(~FeatureType)

###Short/Long Fragment Ratio per Genomic Bin ------
library(BSgenome.Hsapiens.UCSC.hg38)
genome_info_chr21 <- seqinfo(BSgenome.Hsapiens.UCSC.hg38)["chr21"]

bin_size <- 1000000

genome_bins <- tileGenome(genome_info_chr21, tilewidth = bin_size, cut.last.tile.in.chrom = TRUE)

short_frag_min <- 100
short_frag_max <- 150
long_frag_min <- 151
long_frag_max <- 220

all_bin_ratios_list <- list()
for (sample_name in names(all_fragments_list)) {
  fragments_gr <- all_fragments_list[[sample_name]]
  group_name <- unique(fragments_gr$Group) 
  # Calculate ratios 
  sample_bin_ratios <- calculate_bin_ratios(fragments_gr, genome_bins, sample_name, group_name)
  all_bin_ratios_list[[sample_name]] <- sample_bin_ratios
}

all_bin_ratios_df <- bind_rows(all_bin_ratios_list)

ggplot(all_bin_ratios_df, aes(x = start / 1e6, y = Ratio, color = Group)) +
  geom_line(aes(group = Group)) + # Use group aesthetic for lines
  geom_smooth(method = "loess", se = FALSE) + # Add a smoothed trend line
  labs(title = paste0("Median Short/Long Fragment Ratio across Chr21"),
       x = "Genomic Position (Mb)",
       y = "Short/Long Fragment Ratio") +
  theme_minimal()

### Stats -------
wilcox.test(Ratio ~ Group, data = all_bin_ratios_df)




#Fragment End Analysis ------------------------------------------------

library(ggseqlogo) # For plotting sequence logos
motif_length <- 6
genome <- BSgenome.Hsapiens.UCSC.hg38
valid_chroms <- "chr21"
all_five_prime_motifs <- list()
all_three_prime_motifs <- list()

# Process Control samples
for (i in seq_along(bam_files_ctrl)) {
  bam_path <- file.path(data.dir,"bam",bam_files_ctrl[i])
  sample_name <- ctrl.samples[i]
  group_label <- "Control"
  
  motifs <- get_fragment_end_motifs(bam_path, ref_genome = genome, motif_len = motif_length, valid_chroms = valid_chroms)
  all_five_prime_motifs[[paste0(sample_name, "_", group_label)]] <- motifs$five_prime_motifs
  all_three_prime_motifs[[paste0(sample_name, "_", group_label)]] <- motifs$three_prime_motifs
}

# Process ALS samples
for (i in seq_along(bam_files_als)) {
  bam_path <- file.path(data.dir,"bam",bam_files_als[i])
  sample_name <- als.sample[i]
  group_label <- "ALS"
  
  motifs <- get_fragment_end_motifs(bam_path, ref_genome = genome, motif_len = motif_length, valid_chroms = valid_chroms)
  all_five_prime_motifs[[paste0(sample_name, "_", group_label)]] <- motifs$five_prime_motifs
  all_three_prime_motifs[[paste0(sample_name, "_", group_label)]] <- motifs$three_prime_motifs
}

## Motif Counting and Comparison ----

# --- Pool motifs by Group ---
### --- For 5' ends ----
filtered_5prime_control_list <- filter_empty_DNAStringSets(all_five_prime_motifs[grep("_Control", names(all_five_prime_motifs))])
pooled_5prime_control <- c(filtered_5prime_control_list[[1]],filtered_5prime_control_list[[2]],
                           filtered_5prime_control_list[[3]],filtered_5prime_control_list[[4]],
                           filtered_5prime_control_list[[5]],filtered_5prime_control_list[[6]],
                           filtered_5prime_control_list[[7]],filtered_5prime_control_list[[8]],
                           filtered_5prime_control_list[[9]],filtered_5prime_control_list[[10]])
filtered_5prime_als_list <- filter_empty_DNAStringSets(all_five_prime_motifs[grep("_ALS", names(all_five_prime_motifs))])

pooled_5prime_als <- c(filtered_5prime_als_list[[1]],filtered_5prime_als_list[[2]],
                       filtered_5prime_als_list[[3]],filtered_5prime_als_list[[4]],
                       filtered_5prime_als_list[[5]],filtered_5prime_als_list[[6]],
                       filtered_5prime_als_list[[7]],filtered_5prime_als_list[[8]],
                       filtered_5prime_als_list[[9]],filtered_5prime_als_list[[10]],
                       filtered_5prime_als_list[[11]],filtered_5prime_als_list[[12]])


### For 3' ends ----
# Filter control motifs
filtered_3prime_control_list <- filter_empty_DNAStringSets(all_three_prime_motifs[grep("_Control", names(all_three_prime_motifs))])

pooled_3prime_control <- c(filtered_3prime_control_list[[1]],filtered_3prime_control_list[[2]],
                           filtered_3prime_control_list[[3]],filtered_3prime_control_list[[4]],
                           filtered_3prime_control_list[[5]],filtered_3prime_control_list[[6]],
                           filtered_3prime_control_list[[7]],filtered_3prime_control_list[[8]],
                           filtered_3prime_control_list[[9]],filtered_3prime_control_list[[10]])
filtered_3prime_als_list <- filter_empty_DNAStringSets(all_three_prime_motifs[grep("_ALS", names(all_three_prime_motifs))])
pooled_3prime_als <- c(filtered_3prime_als_list[[1]],filtered_3prime_als_list[[2]],
                       filtered_3prime_als_list[[3]],filtered_3prime_als_list[[4]],
                       filtered_3prime_als_list[[5]],filtered_3prime_als_list[[6]],
                       filtered_3prime_als_list[[7]],filtered_3prime_als_list[[8]],
                       filtered_3prime_als_list[[9]],filtered_3prime_als_list[[10]],
                       filtered_3prime_als_list[[11]],filtered_3prime_als_list[[12]])


## Count k-mer  ----
# For 5' ends
counts_5prime_control <- oligonucleotideFrequency(pooled_5prime_control, width = motif_length, simplify.as="matrix")
counts_5prime_als <- oligonucleotideFrequency(pooled_5prime_als, width = motif_length, simplify.as="matrix")

# Convert to data frames for easier manipulation
df_5prime_control <- data.frame(Motif = colnames(counts_5prime_control), Count_Control = colSums(counts_5prime_control))
df_5prime_als <- data.frame(Motif = colnames(counts_5prime_als), Count_ALS = colSums(counts_5prime_als))

TRUE_Total_5prime_Control <- length(pooled_5prime_control)
TRUE_Total_5prime_ALS <- length(pooled_5prime_als)
TRUE_Total_3prime_Control <- length(pooled_3prime_control)
TRUE_Total_3prime_ALS <- length(pooled_3prime_als)
merged_5prime <- full_join(df_5prime_control, df_5prime_als, by = "Motif") %>%
  mutate(Count_Control = replace_na(Count_Control, 0),
         Count_ALS = replace_na(Count_ALS, 0),
         Total_Control = TRUE_Total_5prime_Control, # <--- CORRECTED
         Total_ALS = TRUE_Total_5prime_ALS,         # <--- CORRECTED
         Prop_Control = Count_Control / Total_Control,
         Prop_ALS = Count_ALS / Total_ALS) %>%
  rowwise() %>%
  mutate(Log2FC = log2((Prop_ALS + 1e-9) / (Prop_Control + 1e-9))) %>%
  ungroup() %>%
  arrange(desc(abs(Log2FC)))
# Merge and calculate proportions/log2FoldChange for 5' ends


cat("\nTop differentially abundant 5' end motifs (by absolute Log2FC):\n")
print(head(merged_5prime, 10))

# For 3' ends (similar process)
counts_3prime_control <- oligonucleotideFrequency(pooled_3prime_control, width = motif_length, simplify.as="matrix")
counts_3prime_als <- oligonucleotideFrequency(pooled_3prime_als, width = motif_length, simplify.as="matrix")

df_3prime_control <- data.frame(Motif = colnames(counts_3prime_control), Count_Control = colSums(counts_3prime_control))
df_3prime_als <- data.frame(Motif = colnames(counts_3prime_als), Count_ALS = colSums(counts_3prime_als))

merged_3prime <- full_join(df_3prime_control, df_3prime_als, by = "Motif") %>%
  mutate(Count_Control = replace_na(Count_Control, 0),
         Count_ALS = replace_na(Count_ALS, 0),
         Total_Control = TRUE_Total_3prime_Control,
         Total_ALS = TRUE_Total_3prime_ALS, 
         Prop_Control = Count_Control / Total_Control,
         Prop_ALS = Count_ALS / Total_ALS) %>%
  rowwise() %>%
  mutate(Log2FC = log2((Prop_ALS + 1e-9) / (Prop_Control + 1e-9))) %>%
  ungroup() %>%
  arrange(desc(abs(Log2FC)))
cat("\nTop differentially abundant 3' end motifs (by absolute Log2FC):\n")
print(head(merged_3prime, 10))

## Stats (Fisher's Exact Test) ----
# Perform Fisher's Exact Test for all motifs and apply multiple testing correction.
# This loop will run for all motifs.

# For 5' motifs:
fisher_results_5prime <- merged_5prime %>%
  rowwise() %>%
  mutate(
    p_value = tryCatch({
      test_table <- matrix(c(Count_Control, Total_Control - Count_Control,
                             Count_ALS, Total_ALS - Count_ALS), nrow = 2, byrow = TRUE)
      fisher.test(test_table)$p.value
    }, error = function(e) NA) # Handle potential errors if counts are problematic for test
  ) %>%
  ungroup() %>%
  mutate(p_adj = p.adjust(p_value, method = "BH")) %>% # Benjamini-Hochberg FDR correction
  arrange(p_adj)

cat("\nTop 5' end motifs by adjusted p-value:\n")
print(head(fisher_results_5prime, 10))
write.csv(x = fisher_results_5prime,file = file.path(fig.dir,"Motif_5primeDE.csv"))
## For 3' motifs:
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

cat("\nTop 3' end motifs by adjusted p-value:\n")
print(head(fisher_results_3prime, 10))
write.csv(x = fisher_results_3prime,file = file.path(fig.dir,"Motif_3primeDE.csv"))


##Plot End Distributions----


# Top 5' Enriched in ALS 
top_als_5prime_motifs <- fisher_results_5prime %>%
  filter(Log2FC > 0 ) %>% 
  arrange(desc(Log2FC)) %>%
  head(5) %>% pull(Motif)

if (length(top_als_5prime_motifs) > 0) {
  p_5prime_als <- ggplot() + geom_logo(top_als_5prime_motifs) +
    labs(title = paste0("Top 5' End Motifs Enriched in ALS (", motif_length, "bp)")) + theme_logo()
  pdf(file.path(fig.dir,"5prime_ALSEndMotifsEnriched.pdf"))
  print(p_5prime_als)
  dev.off()
} else {
  cat("No 5' motifs enriched in ALS found for plotting.\n")
}

# Top 5' Enriched in Control 
top_control_5prime_motifs <- fisher_results_5prime %>%
  filter(Log2FC < 0 ) %>% 
  arrange(Log2FC) %>%
  head(5) %>% pull(Motif)

if (length(top_control_5prime_motifs) > 0) {
  p_5prime_control <- ggplot() + geom_logo(top_control_5prime_motifs) +
    labs(title = paste0("Top 5' End Motifs Enriched in Control (", motif_length, "bp)")) + theme_logo()
  print(p_5prime_control)
  pdf(file.path(fig.dir,"5prime_CtrlEndMotifsEnriched.pdf"))
  print(p_5prime_control)
  dev.off()
} else {
  cat("No 5' motifs enriched in Control found for plotting.\n")
}

# Top 3' Enriched in ALS (positive Log2FC, filter by p_adj if desired)
top_als_3prime_motifs <- fisher_results_3prime %>%
  filter(Log2FC > 0 ) %>% 
  arrange(desc(Log2FC)) %>%
  head(5) %>% pull(Motif)

if (length(top_als_3prime_motifs) > 0) {
  p_3prime_als <- ggplot() + geom_logo(top_als_3prime_motifs) +
    labs(title = paste0("Top 3' End Motifs Enriched in ALS (", motif_length, "bp)")) + theme_logo()
  pdf(file.path(fig.dir,"3prime_ALSEndMotifsEnriched.pdf"))
  print(p_3prime_als)
  dev.off()
  
} else {
  cat("No 3' motifs enriched in ALS found for plotting.\n")
}

# Top 3' Enriched in Control (negative Log2FC, filter by p_adj if desired)
top_control_3prime_motifs <- fisher_results_3prime %>%
  filter(Log2FC < 0 & p_adj < 0.05) %>%
  arrange(Log2FC) %>%
  head(5) %>% pull(Motif)

if (length(top_control_3prime_motifs) > 0) {
  p_3prime_control <- ggplot() + geom_logo(top_control_3prime_motifs) +
    labs(title = paste0("Top 3' End Motifs Enriched in Control (", motif_length, "bp)")) + theme_logo()
  pdf(file.path(fig.dir,"3prime_CtrlEndMotifsEnriched.pdf"))
  print(p_3prime_control)
  dev.off()
} else {
  cat("No 3' motifs enriched in Control found for plotting.\n")
}

cat("\nFragment End Motif Analysis Complete!\n")
cat("Check your working directory for generated sequence logo plots (.png files).\n")


#Start end site Dist ------
flanking_window <- 1000
bin_size_bp <- 6

tss_gr <- transcripts(txdb)
tss_gr <- resize(tss_gr, width = 1, fix = "start") # Get 1bp at TSS
tss_gr <- tss_gr[seqnames(tss_gr) == valid_chroms]
tss_gr <- keepSeqlevels(tss_gr, valid_chroms, pruning.mode="coarse")


promoter_centers <- promoters_chr21
start(promoter_centers) <- start(promoters_chr21) + floor(width(promoters_chr21) / 2)
end(promoter_centers) <- start(promoter_centers)

seqinfo(promoter_centers) <- seqinfo(genome)[seqlevels(promoter_centers)]

all_end_pos_dfs <- list()

for (sample_name in names(all_fragments_list)) {
  fragments_gr_sample <- all_fragments_list[[sample_name]]
  
  # Process for Promoters (TSS centers)
  prom_end_pos_df <- calculate_relative_distances(
    fragments_gr_sample, promoter_centers, "Promoter_TSS", window_size = flanking_window
  )
  if (!is.null(prom_end_pos_df)) {
    all_end_pos_dfs[[paste0(sample_name, "_Promoter")]] <- prom_end_pos_df
  }
  
  # You can add other features here if you have appropriate 'center' GRanges for them
  # E.g., for CpG Islands, you might use their centers:
  # cpg_centers <- cpg_islands_chr21
  # start(cpg_centers) <- start(cpg_islands_chr21) + floor(width(cpg_islands_chr21) / 2)
  # end(cpg_centers) <- start(cpg_centers)
  # cpg_end_pos_df <- calculate_relative_distances(fragments_gr_sample, cpg_centers, "CpG_Island", window_size = flanking_window)
  # if (!is.null(cpg_end_pos_df)) { all_end_pos_dfs[[paste0(sample_name, "_CpG")]] <- cpg_end_pos_df }
}

final_end_pos_df <- bind_rows(all_end_pos_dfs)

