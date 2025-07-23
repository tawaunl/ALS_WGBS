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
max_insert_size <- 1000 
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
pdf(file.path(fig.dir,"InsertSizeDensity_byGroup.pdf"))
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

write.csv(sample_moments,
          file = file.path(fig.dir,"fragment_length_summary_statistics_chr21.csv"),
          row.names = FALSE)

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
filtered_5prime_control_list <- filter_empty_DNAStringSets(all_five_prime_motifs[grep("_Control",
                                                                                      names(all_five_prime_motifs))])

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

# --- For 3' ends ---
# Filter control motifs
filtered_3prime_control_list <- filter_empty_DNAStringSets(all_three_prime_motifs[grep("_Control", names(all_three_prime_motifs))])

pooled_3prime_control <- c(filtered_3prime_control_list[[1]],filtered_3prime_control_list[[2]],
                           filtered_3prime_control_list[[3]],filtered_3prime_control_list[[4]],
                           filtered_3prime_control_list[[5]],filtered_3prime_control_list[[6]],
                           filtered_3prime_control_list[[7]],filtered_3prime_control_list[[8]],
                           filtered_3prime_control_list[[9]],filtered_3prime_control_list[[10]])


# Filter ALS motifs
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


# Merge and calculate proportions/log2FoldChange for 5' ends
merged_5prime <- full_join(df_5prime_control, df_5prime_als, by = "Motif") %>%
  mutate(Count_Control = replace_na(Count_Control, 0), # Replace NA counts with 0 (motif not found in one group)
         Count_ALS = replace_na(Count_ALS, 0),
         Total_Control = sum(Count_Control), 
         Total_ALS = sum(Count_ALS),
         Prop_Control = Count_Control / Total_Control, 
         Prop_ALS = Count_ALS / Total_ALS) %>%        
  rowwise() %>% 
  mutate(Log2FC = log2((Prop_ALS + 1e-9) / (Prop_Control + 1e-9))) %>% # Add small pseudocount to avoid log(0) issues
  ungroup() %>%
  arrange(desc(abs(Log2FC))) # Sort 

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
         Total_Control = sum(Count_Control),
         Total_ALS = sum(Count_ALS), 
         Prop_Control = Count_Control / Total_Control,
         Prop_ALS = Count_ALS / Total_ALS) %>%
  rowwise() %>%
  mutate(Log2FC = log2((Prop_ALS + 1e-9) / (Prop_Control + 1e-9))) %>%
  ungroup() %>%
  arrange(desc(abs(Log2FC)))


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
write.csv(fisher_results_5prime,
          file = file.path(fig.dir,"Motif_5primeDE.csv"),row.names = FALSE)
# For 3' motifs:
cat("\nRunning Fisher's Exact Tests for all 3' end motifs and adjusting p-values...\n")
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
write.csv(fisher_results_3prime,
          file = file.path(fig.dir,"Motif_3primeDE.csv"),row.names = FALSE)

##Plot----


# Top 5' Enriched in ALS 
top_als_5prime_motifs <- fisher_results_5prime %>%
  filter(Log2FC > 0 ) %>% 
  arrange(desc(Log2FC)) %>%
  head(5) %>% pull(Motif)

if (length(top_als_5prime_motifs) > 0) {
  p_5prime_als <- ggplot() + geom_logo(top_als_5prime_motifs) +
    labs(title = paste0("Top 5' End Motifs Enriched in ALS (", motif_length, "bp)")) + theme_logo()
  print(p_5prime_als)
  ggsave(file.path(fig.dir,paste0("top_5prime_als_motifs_", motif_length, "bp.png")), plot = p_5prime_als, width = 8, height = 4, dpi = 300)
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
  ggsave(file.path(fig.dir,paste0("top_5prime_control_motifs_", motif_length, "bp.png")), plot = p_5prime_control, width = 8, height = 4, dpi = 300)
} else {
  cat("No 5' motifs enriched in Control found for plotting.\n")
}

# Top 3' Enriched in ALS (positive Log2FC, filter by p_adj if desired)
top_als_3prime_motifs <- fisher_results_3prime %>%
  filter(Log2FC > 0) %>% 
  arrange(desc(Log2FC)) %>%
  head(5) %>% pull(Motif)

if (length(top_als_3prime_motifs) > 0) {
  p_3prime_als <- ggplot() + geom_logo(top_als_3prime_motifs) +
    labs(title = paste0("Top 3' End Motifs Enriched in ALS (", motif_length, "bp)")) + theme_logo()
  print(p_3prime_als)
  ggsave(file.path(fig.dir,paste0("top_3prime_als_motifs_", motif_length, "bp.png")), plot = p_3prime_als, width = 8, height = 4, dpi = 300)
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
  print(p_3prime_control)
  ggsave(file.path(fig.dir,paste0("top_3prime_control_motifs_", motif_length, "bp.png")), plot = p_3prime_control, width = 8, height = 4, dpi = 300)
} else {
  cat("No 3' motifs enriched in Control found for plotting.\n")
}

cat("\nFragment End Motif Analysis Complete!\n")
cat("Check your working directory for generated sequence logo plots (.png files).\n")

# Start End Analysis ------------
flanking_window <- 1000
bin_size_bp <- 5
tss_gr <- transcripts(txdb)
tss_gr <- resize(tss_gr, width = 1, fix = "start") # Get 1bp at TSS
tss_gr <- tss_gr[seqnames(tss_gr) == valid_chroms]
seqlevelsStyle(tss_gr) <- "UCSC"
tss_gr <- keepSeqlevels(tss_gr, valid_chroms, pruning.mode="coarse")


promoter_centers <- promoters_chr21
start(promoter_centers) <- start(promoters_chr21) + floor(width(promoters_chr21) / 2)
end(promoter_centers) <- start(promoter_centers) # Make it a 1bp point feature
seqinfo(promoter_centers) <- seqinfo(genome)[seqlevels(promoter_centers)]

cpg_centers <- cpg_islands_chr21
start(cpg_centers) <- start(cpg_islands_chr21) + floor(width(cpg_islands_chr21) / 2)
end(cpg_centers) <- start(cpg_centers)
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
  
}

final_end_pos_df <- bind_rows(all_end_pos_dfs)
## Plot ------
if (nrow(final_end_pos_df) > 0) {
  ggplot(final_end_pos_df, aes(x = RelativePosition, color = Group, fill = Group)) +
    geom_histogram(aes(y = after_stat(density)), binwidth = bin_size_bp, alpha = 0.3, position = "identity") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") + # Mark the center of the feature
    facet_grid(EndType ~ FeatureType, scales = "free_y") + # Separate by 5'/3' and Feature Type
    labs(title = paste0("Fragment End Distribution around Genomic Features (", flanking_window, "bp window)"),
         x = paste0("Distance from Feature Center (bp)"),
         y = "Density") +
    theme_minimal() +
    theme(legend.position = "bottom") +
    xlim(-flanking_window, flanking_window) # Ensure consistent x-axis limits
  ggsave(file.path(fig.dir,"fragment_end_distributions_chr21.png"), width = 10, height = 8, dpi = 300)
} else {
  cat("No fragment ends found overlapping specified features for this analysis.\n")
}
### Stats --------
breaks <- seq(-flanking_window, flanking_window, by = bin_size_bp)

labels <- seq(-flanking_window + bin_size_bp/2, flanking_window - bin_size_bp/2, by = bin_size_bp)

# Bin the data and count fragments per bin per group/end_type/feature_type
binned_counts <- final_end_pos_df %>%
  mutate(Bin = cut(RelativePosition,
                   breaks = breaks,
                   labels = labels, # Assign center of bin as label
                   include.lowest = TRUE, # Include the lowest value in the first bin
                   right = TRUE)) %>% # Bins are (start, end], adjust if needed
  filter(!is.na(Bin)) %>% # Remove any NA bins (fragments outside defined window)
  group_by(Group, FeatureType, EndType, Bin) %>%
  summarise(Count = n(), .groups = 'drop_last') %>% # Count fragments per bin
  ungroup()

all_combinations <- expand_grid(
  Group = unique(binned_counts$Group),
  FeatureType = unique(binned_counts$FeatureType),
  EndType = unique(binned_counts$EndType),
  Bin = unique(binned_counts$Bin)
)

# Left join to fill in missing bins with 0 counts
binned_counts_full <- left_join(all_combinations, binned_counts, by = c("Group", "FeatureType", "EndType", "Bin")) %>%
  mutate(Count = replace_na(Count, 0)) # Replace NA counts (for empty bins) with 0

# Calculate total fragments per group, feature type, and end type for proportion calculation
# This TotalCount is the sum of all fragments in a given Group/FeatureType/EndType combination
total_counts_per_group_end_feature <- binned_counts_full %>%
  group_by(Group, FeatureType, EndType) %>%
  summarise(TotalCount = sum(Count), .groups = 'drop')

# Join the total counts back to the binned data to calculate proportions
binned_proportions <- left_join(binned_counts_full, total_counts_per_group_end_feature, by = c("Group", "FeatureType", "EndType")) %>%
  mutate(Proportion = Count / TotalCount)


fisher_results_bins <- data.frame(binned_proportions) %>%
  dplyr::select(FeatureType, EndType, Bin, Group, Count, TotalCount) %>%
  pivot_wider(names_from = Group, values_from = c(Count, TotalCount)) %>%
  rowwise() %>% # Process row by row for the Fisher's test
  mutate(
    Count_Control = replace_na(Count_Control, 0), # Ensure counts are 0 if NA
    Count_ALS = replace_na(Count_ALS, 0),
    TotalCount_Control = replace_na(TotalCount_Control, 0), # Ensure totals are 0 if NA (e.g., if a group has no data for that FeatureType/EndType)
    TotalCount_ALS = replace_na(TotalCount_ALS, 0),
    p_value = tryCatch({
      if (TotalCount_Control == 0 && TotalCount_ALS == 0) {
        return(NA_real_)
      }
      # Construct the 2x2 contingency table for the current bin
      test_table <- matrix(c(Count_Control, TotalCount_Control - Count_Control, # Count of motif IN bin vs. OUT of bin for Control
                             Count_ALS, TotalCount_ALS - Count_ALS),             # Count of motif IN bin vs. OUT of bin for ALS
                           nrow = 2, byrow = TRUE)
      # Perform Fisher's Exact Test
      fisher.test(test_table)$p.value
    }, error = function(e) NA_real_) # Return NA if any error occurs during the test
  ) %>%
  ungroup() %>% # Ungroup after row-wise operations
  mutate(p_adj = p.adjust(p_value, method = "BH")) 


cat("\nTop bins with statistically significant differences (p_adj < 0.05):\n")
significant_bins <- fisher_results_bins %>% filter(p_value < 0.05)
print(significant_bins)

if (nrow(final_end_pos_df) > 0 && nrow(significant_bins) > 0) {
  
  # 1. Prepare data for plotting rectangles for significant bins
  # Ensure 'Bin' in significant_bins is treated numerically for xmin/xmax calculation.
  # It might be a factor, so converting to character then numeric is safest.
  highlight_rects <- significant_bins %>%
    mutate(
      xmin = as.numeric(as.character(Bin)) - bin_size_bp / 2, # Calculate start of the bin
      xmax = as.numeric(as.character(Bin)) + bin_size_bp / 2, # Calculate end of the bin
      ymin = -Inf, # Rectangles should span the full y-axis range
      ymax = Inf   # (will be correctly scaled by facet_grid's free_y)
    ) %>%
    # Select only the necessary columns for geom_rect to avoid aesthetic inheritance issues
    dplyr::select(FeatureType, EndType, xmin, xmax, ymin, ymax) %>%
    distinct() # Remove any duplicate rectangles if they somehow arose
  
  # 2. Create the main ggplot
  p <- ggplot(final_end_pos_df, aes(x = RelativePosition, color = Group, fill = Group)) +
    # Add highlight rectangles BEFORE geom_histogram so they are drawn *behind* the density curves
    geom_rect(data = highlight_rects,
              aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
              fill = "red", # Color for the highlight
              alpha = 0.5,        # Transparency of the highlight
              inherit.aes = FALSE) + # IMPORTANT: Prevents geom_rect from trying to inherit Group, color, fill from main aes()
    
    # Add the density profiles using histograms (with density normalization)
    geom_histogram(aes(y = after_stat(density)),
                   binwidth = bin_size_bp,
                   alpha = 0.3,
                   position = "identity") + # 'identity' means no stacking/dodging for fill
    
    # Add a dashed vertical line at the feature center (0 bp)
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
    
    # Facet the plot by EndType (5' vs 3') and FeatureType (Promoter, CpG, etc.)
    facet_grid(EndType ~ FeatureType, scales = "free_y") + # 'free_y' allows y-axes to scale independently per facet
    
    # Add labels and theme
    labs(title = paste0("Fragment End Distribution around Genomic Features (", flanking_window, "bp window)"),
         x = "Distance from Feature Center (bp)",
         y = "Density") +
    theme_minimal() +
    theme(legend.position = "bottom",
          plot.title = element_text(hjust = 0.5)) + # Center plot title
    xlim(-flanking_window, flanking_window) # Ensure consistent x-axis limits
  
  # Print and save the plot
  print(p)
  ggsave(file.path(fig.dir,"fragment_end_distributions_significant_bins_chr21.png"), plot = p, width = 12, height = 8, dpi = 300)
  
} else if (nrow(final_end_pos_df) > 0) {
  # If no significant bins were found, plot the profiles without highlights
  cat("No significant bins found to highlight. Plotting original density profiles.\n")
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
  print(p)
  ggsave(file.path(fig.dir,"fragment_end_distributions_chr21_no_highlights.png"), plot = p, width = 10, height = 8, dpi = 300)
} else {
  # If no data at all, print a message
  cat("No data available in final_end_pos_df to plot fragment end distributions.\n")
}


# Classification --------------
library(caret)   
library(pROC)

feature_data <- sample_moments %>%
  dplyr::select(Sample, Group, Mean_InsertSize, Median_InsertSize,
         SD_InsertSize, Skewness_InsertSize, Kurtosis_InsertSize)

# Join with mono_di_proportions
feature_data <- feature_data %>%
  left_join(mono_di_proportions %>% dplyr::select(Sample, Mono_Proportion, Mono_Di_Ratio), by = "Sample")

# Join with short_fragment_analysis
feature_data <- feature_data %>%
  left_join(short_fragment_props %>% dplyr::select(Sample,ShortFragmentCount,TotalFragments, ShortFragmentProportion), by = "Sample")

# Handle potential NA values if some samples missed a metric (e.g., Mono_Di_Ratio if Di_Count was 0)
# A simple approach is to impute with the mean/median, or remove.
# For simplicity, we'll impute with the column mean for now.
# In a real scenario, carefully consider NA handling (e.g., k-NN imputation, specific to each feature).
for(col_name in names(feature_data)) {
  if(is.numeric(feature_data[[col_name]]) && any(is.na(feature_data[[col_name]]))) {
    feature_data[[col_name]][is.na(feature_data[[col_name]])] <- mean(feature_data[[col_name]], na.rm = TRUE)
  }
}

# Ensure the 'Group' variable is a factor for classification
# And rename levels to be R-friendly (no spaces or special characters)
feature_data$Group <- factor(feature_data$Group)

# Remove the 'Sample' identifier column from the features for training
# We keep it separate for now for clarity
classification_df <- feature_data %>% dplyr::select(-Sample)

cat("Prepared classification data (first 5 rows):\n")
print(head(classification_df, 5))
cat(paste0("Total samples for classification: ", nrow(classification_df), "\n"))
cat(paste0("Number of features: ", ncol(classification_df) - 1, "\n"))


# Use repeated cross-validation for more robust estimates
# method = "repeatedcv": Repeated Cross-Validation
# number = 10: 10-fold CV
# repeats = 5: Repeat 5 times
# classProbs = TRUE: Needed for AUC-ROC calculation
# summaryFunction = twoClassSummary: Calculates AUC, Sensitivity, Specificity
train_control <- trainControl(
  method = "repeatedcv",
  number = 10,
  repeats = 50,
  classProbs = TRUE,
  summaryFunction = twoClassSummary,
  savePredictions = "final" # Save predictions for ROC curve plotting
)

# Train Classification Models ---

cat("\n--- Training Logistic Regression Model ---\n")
# Logistic Regression
model_glm <- train(
  Group ~ ., # Predict 'Group' based on all other features
  data = classification_df,
  method = "glm",
  family = "binomial", # For binary classification
  trControl = train_control,
  metric = "ROC" # Optimize for AUC-ROC
)
cat("Logistic Regression Model Training Complete.\n")
print(model_glm)

cat("\n--- Training Random Forest Model ---\n")
# Random Forest
# tuneLength = 3 means it will try 3 values for mtry (number of variables to randomly sample at each split)
# Alternatively, tuneGrid can be used for specific mtry values.
model_rf <- train(
  Group ~ .,
  data = classification_df,
  method = "rf",
  trControl = train_control,
  metric = "ROC", # Optimize for AUC-ROC
  tuneLength = 4 # Adjust for more thorough tuning if needed (longer runtime)
)
cat("Random Forest Model Training Complete.\n")
print(model_rf)

# Evaluate and Present Metrics ---
print(model_glm$results) # Shows performance for different tuning parameters (if any)
cat(paste0("Best Logistic Regression ROC: ", round(max(model_glm$results$ROC), 3), "\n"))
cat(paste0("Best Logistic Regression Sensitivity: ", round(model_glm$results$Sens[which.max(model_glm$results$ROC)], 3), "\n"))
cat(paste0("Best Logistic Regression Specificity: ", round(model_glm$results$Spec[which.max(model_glm$results$ROC)], 3), "\n"))


print(model_rf$results)
cat(paste0("Best Random Forest ROC: ", round(max(model_rf$results$ROC), 3), "\n"))
cat(paste0("Best Random Forest Sensitivity: ", round(model_rf$results$Sens[which.max(model_rf$results$ROC)], 3), "\n"))
cat(paste0("Best Random Forest Specificity: ", round(model_rf$results$Spec[which.max(model_rf$results$ROC)], 3), "\n"))

# Logistic Regression Confusion Matrix
if ("mtry" %in% names(model_glm$bestTune)) {
  glm_pred <- model_glm$pred %>% filter(mtry == model_glm$bestTune$mtry)
} else { # For glm, there is no mtry, so just use all predictions
  glm_pred <- model_glm$pred
}

cat("\nLogistic Regression Confusion Matrix:\n")
cm_glm <- confusionMatrix(data = glm_pred$pred, reference = glm_pred$obs, positive = "ALS")
print(cm_glm)
cat(paste0("Logistic Regression Precision: ", round(cm_glm$byClass["Pos Pred Value"], 3), "\n"))
cat(paste0("Logistic Regression F1-Score: ", round(cm_glm$byClass["F1"], 3), "\n"))


# Random Forest Confusion Matrix
rf_pred <- model_rf$pred %>% filter(mtry == model_rf$bestTune$mtry)

cat("\nRandom Forest Confusion Matrix:\n")
cm_rf <- confusionMatrix(data = rf_pred$pred, reference = rf_pred$obs, positive = "ALS")
print(cm_rf)
cat(paste0("Random Forest Precision: ", round(cm_rf$byClass["Pos Pred Value"], 3), "\n"))
cat(paste0("Random Forest F1-Score: ", round(cm_rf$byClass["F1"], 3), "\n"))


## ROC Curve and AUC -----


# For Logistic Regression
roc_glm <- roc(response = glm_pred$obs, predictor = glm_pred$ALS, levels = c("Control", "ALS"))
cat(paste0("Logistic Regression AUC: ", round(auc(roc_glm), 3), "\n"))
plot(roc_glm, main = "ROC Curve - Logistic Regression", print.auc = TRUE, col = "blue")
ggsave(file.path(fig.dir,"roc_curve_glm.png"), plot = last_plot(), width = 6, height = 6, dpi = 300)

# For Random Forest
roc_rf <- roc(response = rf_pred$obs, predictor = rf_pred$ALS, levels = c("Control", "ALS"))
cat(paste0("Random Forest AUC: ", round(auc(roc_rf), 3), "\n"))
plot(roc_rf, main = "ROC Curve - Random Forest", print.auc = TRUE, col = "red")
ggsave(file.path(fig.dir,"roc_curve_rf.png"), plot = last_plot(), width = 6, height = 6, dpi = 300)

# Compare ROC curves on one plot
plot(roc_glm, col = "blue", main = "ROC Curves: ALS vs Control")
lines(roc_rf, col = "red")
legend("bottomright", legend = c(paste0("Logistic Regression (AUC = ", round(auc(roc_glm), 2), ")"),
                                 paste0("Random Forest (AUC = ", round(auc(roc_rf), 2), ")")),
       col = c("blue", "red"), lty = 1, cex = 0.8)
ggsave(file.path(fig.dir,"roc_curves_combined.png"), plot = last_plot(), width = 7, height = 6, dpi = 300)


# --- 6. Feature Importance (for Random Forest) ---
if (!is.null(varImp(model_rf)$importance)) {
  rf_importance <- varImp(model_rf)$importance
  # For classification, importance is often based on MeanDecreaseAccuracy or MeanDecreaseGini
  # For "ROC" metric, caret's varImp typically uses a weighted average or permutations
  print(rf_importance)
  
  # Plotting importance
  plot_importance <- ggplot(rf_importance, aes(x = reorder(rownames(rf_importance), Overall), y = Overall)) +
    geom_bar(stat = "identity", fill = "steelblue") +
    coord_flip() + # Flip coordinates to make labels readable
    labs(title = "Random Forest Feature Importance",
         x = "Feature",
         y = "Importance (Overall)") +
    theme_minimal()
  print(plot_importance)
  ggsave(file.path(fig.dir,"rf_feature_importance.png"), plot = plot_importance, width = 8, height = 6, dpi = 300)
} else {
  cat("Feature importance not available or could not be computed for Random Forest model.\n")
}

### Save metrics ---------
glm_overall_metrics <- as.data.frame(t(cm_glm$overall))
glm_overall_metrics <- glm_overall_metrics %>%
  dplyr::select(Accuracy, Kappa) %>%
  rename(Overall_Accuracy = Accuracy, Overall_Kappa = Kappa)

# By-class metrics (for the positive class "ALS")
glm_byclass_metrics <- as.data.frame(t(cm_glm$byClass))
glm_byclass_metrics <- glm_byclass_metrics %>%
  dplyr::select(`Sensitivity`, `Specificity`, `Pos Pred Value`, `Neg Pred Value`, `F1`) %>%
  rename(Sensitivity = `Sensitivity`, Specificity = `Specificity`,
         Precision = `Pos Pred Value`, Negative_Predictive_Value = `Neg Pred Value`, F1_Score = `F1`)

# AUC-ROC
glm_auc <- data.frame(AUC = as.numeric(auc(roc_glm)))

glm_metrics_df <- bind_cols(glm_overall_metrics, glm_byclass_metrics, glm_auc) %>%
  mutate(Model = "Logistic Regression") %>%
  dplyr::select(Model, everything()) # Reorder to put Model first


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
  dplyr::select(Model, everything()) # Reorder to put Model first


all_classification_metrics <- bind_rows(glm_metrics_df, rf_metrics_df)

output_metrics_csv_file <- "classification_performance_metrics.csv"
write.csv(all_classification_metrics,
          file =file.path(fig.dir, output_metrics_csv_file), row.names = FALSE)

# Methelaytion analysis --------------------------------------
# # need to install sam tools from command line using homebrew
# 
# sample_groups_map <- feature_data$Group
# names(sample_groups_map) <-feature_data$Sample
# bismark_extractor_path <- "../Bismark-0.22.3/bismark_methylation_extractor"
# base_methylation_output_dir <- "methylation_reports_chr21"
# if (!dir.exists(file.path(fig.dir,base_methylation_output_dir))) {
#   dir.create(file.path(fig.dir,base_methylation_output_dir), recursive = TRUE)
# }
# methylation_report_files <- list()
# 
# for (bam_path in all_bam_files) {
#   sample_id <- strsplit(bam_path,".d")[[1]][1]
#   output_dir_for_sample <- file.path(file.path(fig.dir,base_methylation_output_dir),
#                                      tools::file_path_sans_ext(sample_id, compress = FALSE))
#   
#   # Create output directory for the current sample
#   if (!dir.exists(output_dir_for_sample)) {
#     dir.create(output_dir_for_sample, recursive = TRUE)
#     
#   }
#   # Construct the Bismark command arguments
#   # Ensure all paths are absolute to avoid issues with R's working directory
#   command_args <- c(
#     "--paired-end",
#     "--comprehensive",
#     "--no_header", # Often useful if you just want data
#     "-o", shQuote(output_dir_for_sample), # shQuote handles spaces in paths
#     shQuote(file.path(data.dir,"bam",bam_path))
#   )
#   
#   cat(paste0("Processing: ", sample_id, "\n"))
#   cat(paste0("  Command: ", bismark_extractor_path, " ", paste(command_args, collapse = " "), "\n"))
#   
#   # Execute the command
#   # intern=TRUE captures stdout/stderr, wait=TRUE ensures R waits for command to finish
#   # stdout/stderr are useful for debugging.
#   result <- system2(
#     command = bismark_extractor_path,
#     args = command_args,
#     stdout = TRUE,
#     stderr = TRUE,
#     wait = TRUE
#   )
#   
#   # Check for successful execution (status code 0 usually means success)
#   # Note: system2 doesn't return the exit status directly if stdout/stderr are captured.
#   # You'd need to check the 'result' for error messages.
#   if (!any(grepl("Error|Fail", result, ignore.case = TRUE))) { # Simple check for common errors in output
#     cat(paste0("  Successfully processed ", sample_id, "\n"))
#     # Store the path to the CpG report file
#     # Bismark names them consistently: [original_bam_name]_CpG_report.txt.gz
#     cpg_report_file <- file.path(output_dir_for_sample, paste0(tools::file_path_sans_ext(sample_id, compress = FALSE), "_CpG_report.txt.gz"))
#     
#     # Check if the file actually exists and is not empty (for robustness)
#     if (file.exists(cpg_report_file) && file.size(cpg_report_file) > 0) {
#       methylation_report_files[[sample_id]] <- cpg_report_file
#     } else {
#       warning(paste0("CpG report file not found or empty for ", sample_id, ": ", cpg_report_file))
#     }
#     
#   } else {
#     stop(paste0("Error processing ", sample_id, ". Bismark output:\n", paste(result, collapse = "\n")))
#   }
# }
# 
# cat("\nAll Bismark methylation extraction commands executed.\n")
# cat("\nList of generated CpG report files:\n")
# print(unlist(methylation_report_files))
# 
