# useful functions
get_insert_sizes <- function(bam_file_path, min_len = 0, max_len = 1000) {
  require(Rsamtools)
  print(paste("Processing:", bam_file_path))
  
  # Define what fields to read from the BAM file
  # "qname": Query Name
  # "flag": SAM flags
  # "rname": Reference Name
  # "mpos": 1-based leftmost mapping position
  # "isize": Observed Template Length (insert size)
  param <- ScanBamParam(what = c("qname", "flag", "rname", "mpos", "isize"),
                        flag = scanBamFlag(isPaired = TRUE,
                                           isProperPair = TRUE, # Ensure proper pairs
                                           isDuplicate = FALSE)) # Optional: filter out duplicates
  
  # Read the BAM file
  bam_data <- scanBam(bam_file_path, param = param)[[1]]
  
  # Filter for reads where tlen is not NA and is within the desired range
  # tlen (TLEN) can be negative if mate is to the left, take absolute value
  insert_sizes <- abs(bam_data$isize)
  
  # Filter by min/max length
  insert_sizes <- insert_sizes[!is.na(insert_sizes) &
                                 insert_sizes >= min_len &
                                 insert_sizes <= max_len]
  
  return(insert_sizes)
}


remove_outliers_iqr <- function(data, variable_col, group_col, k = 1.5) {
  variable_name <- deparse(substitute(variable_col)) # Get the name of the column
  
  data_cleaned <- data %>%
    group_by({{ group_col }}) %>% # Group by the specified grouping column
    mutate(
      Q1 = quantile({{ variable_col }}, 0.25, na.rm = TRUE),
      Q3 = quantile({{ variable_col }}, 0.75, na.rm = TRUE),
      IQR = Q3 - Q1,
      lower_bound = Q1 - k * IQR,
      upper_bound = Q3 + k * IQR,
      is_outlier = !({{ variable_col }} >= lower_bound & {{ variable_col }} <= upper_bound)
    ) %>%
    ungroup()
  
  # Report identified outliers
  outliers_info <- data_cleaned %>%
    filter(is_outlier) %>%
    dplyr::select({{ group_col }}, Sample, {{ variable_col }}) # Select relevant columns for reporting
  if (nrow(outliers_info) > 0) {
    cat(paste0("\n--- Outliers Identified for '", variable_name, "' (k=", k, ") ---\n"))
    print(outliers_info)
    cat(paste0("Total outliers removed: ", nrow(outliers_info), "\n"))
  } else {
    cat(paste0("\n--- No outliers identified for '", variable_name, "' (k=", k, ") ---\n"))
  }
  
  # Return data without outliers
  data_cleaned_filtered <- data_cleaned %>%
    filter(!is_outlier) %>%
    dplyr::select(-Q1, -Q3, -IQR, -lower_bound, -upper_bound, -is_outlier) # Remove helper columns
  
  return(data_cleaned_filtered)
}

identify_outliers_iqr <- function(data, variable_col, group_col, k = 1.5) {
  variable_name <- deparse(substitute(variable_col)) # Get the name of the column
  
  data_with_outliers_flag <- data %>%
    group_by({{ group_col }}) %>% # Group by the specified grouping column
    mutate(
      Q1 = quantile({{ variable_col }}, 0.25, na.rm = TRUE),
      Q3 = quantile({{ variable_col }}, 0.75, na.rm = TRUE),
      IQR = Q3 - Q1,
      lower_bound = Q1 - k * IQR,
      upper_bound = Q3 + k * IQR,
      is_outlier = !({{ variable_col }} >= lower_bound & {{ variable_col }} <= upper_bound)
    ) %>%
    ungroup() %>%
    # Remove helper columns if you don't need them in the final df, but keep is_outlier
    dplyr::select(-Q1, -Q3, -IQR, -lower_bound, -upper_bound)
  
  return(data_with_outliers_flag)
}

calculate_relative_distances <- function(fragments_gr, feature_centers_gr, feature_type_name, window_size = 1000) {
  cat(paste0("  Calculating distances for ", unique(fragments_gr$Sample), " around ", feature_type_name, "...\n"))
  
  # Find overlaps between fragments and the expanded window around features
  # Expand features by window_size on both sides
  expanded_features <- GenomicRanges::flank(feature_centers_gr, width = window_size, start = TRUE)
  expanded_features <- c(expanded_features, GenomicRanges::flank(feature_centers_gr, width = window_size, start = FALSE))
  expanded_features <- reduce(expanded_features) # Merge overlapping expanded regions for efficiency
  
  # Find overlaps of fragment ends with expanded feature regions
  # We are interested in fragment START and END positions
  # 5' ends of fragments
  five_prime_points <- resize(fragments_gr, width = 1, fix = "start")
  # 3' ends of fragments
  three_prime_points <- resize(fragments_gr, width = 1, fix = "end")
  
  # Only consider points that overlap the expanded feature regions
  five_prime_overlaps <- findOverlaps(five_prime_points, expanded_features)
  three_prime_overlaps <- findOverlaps(three_prime_points, expanded_features)
  
  # Get the coordinates of the overlapping fragment ends
  five_prime_coords <- start(five_prime_points[queryHits(five_prime_overlaps)])
  three_prime_coords <- start(three_prime_points[queryHits(three_prime_overlaps)])
  
  # Get the coordinates of the feature centers
  feature_center_coords_for_5prime <- start(feature_centers_gr[subjectHits(five_prime_overlaps)])
  feature_center_coords_for_3prime <- start(feature_centers_gr[subjectHits(three_prime_overlaps)])
  
  # Calculate relative distances
  relative_distances_5prime <- five_prime_coords - feature_center_coords_for_5prime
  relative_distances_3prime <- three_prime_coords - feature_center_coords_for_3prime
  
  # Filter to be within the plotting window
  relative_distances_5prime <- relative_distances_5prime[
    relative_distances_5prime >= -window_size & relative_distances_5prime <= window_size
  ]
  relative_distances_3prime <- relative_distances_3prime[
    relative_distances_3prime >= -window_size & relative_distances_3prime <= window_size
  ]
  
  # Combine into a data frame
  if (length(relative_distances_5prime) > 0 || length(relative_distances_3prime) > 0) {
    df_5prime <- data.frame(
      RelativePosition = relative_distances_5prime,
      EndType = "5' End",
      Sample = unique(fragments_gr$Sample),
      Group = unique(fragments_gr$Group),
      FeatureType = feature_type_name
    )
    df_3prime <- data.frame(
      RelativePosition = relative_distances_3prime,
      EndType = "3' End",
      Sample = unique(fragments_gr$Sample),
      Group = unique(fragments_gr$Group),
      FeatureType = feature_type_name
    )
    return(bind_rows(df_5prime, df_3prime))
  } else {
    return(NULL)
  }
}

get_fragments_granges <- function(bam_file_path, max_len = 1000) {
  print(paste("Extracting fragments from:", bam_file_path))
  
  param <- ScanBamParam(
    what = c("qname", "flag", "rname", "pos", "isize"),
    flag = scanBamFlag(isPaired = TRUE, isProperPair = TRUE, isDuplicate = FALSE)
  )
  
  gal_pairs <- readGAlignmentPairs(bam_file_path, param = param)
  
  fragments_gr <- granges(gal_pairs)
  fragments_gr$insert_size <- width(fragments_gr)
  fragments_gr <- fragments_gr[fragments_gr$insert_size <= max_len]
  
  # Ensure fragments also only contain chr21 levels (if they don't already)
  # This step might not be strictly necessary if BAM is already chr21-only
  # But it's good for consistency with feature files.
  fragments_gr <- keepSeqlevels(fragments_gr, "chr21", pruning.mode="coarse")
  
  
  return(fragments_gr)
}

get_feature_insert_sizes <- function(fragments_gr, feature_gr, feature_type_name) {
  # Find overlaps
  overlaps <- findOverlaps(fragments_gr, feature_gr)
  
  # Get fragments that overlap the feature
  fragments_in_feature <- fragments_gr[queryHits(overlaps)]
  
  # Extract their insert sizes
  insert_sizes <- fragments_in_feature$insert_size
  
  # Return a data frame for plotting
  if (length(insert_sizes) > 0) {
    return(data.frame(
      InsertSize = insert_sizes,
      Sample = unique(fragments_in_feature$Sample), # Should be only one sample
      Group = unique(fragments_in_feature$Group),   # Should be only one group
      FeatureType = feature_type_name
    ))
  } else {
    return(NULL) # No fragments found
  }
}


calculate_bin_ratios <- function(fragments_gr, genome_bins_gr, sample_name, group_name) {
  # Prepare fragments by length category
  fragments_short <- fragments_gr[fragments_gr$insert_size >= short_frag_min & fragments_gr$insert_size <= short_frag_max]
  fragments_long <- fragments_gr[fragments_gr$insert_size >= long_frag_min & fragments_gr$insert_size <= long_frag_max]
  
  # Count overlaps per bin
  short_counts <- countOverlaps(genome_bins_gr, fragments_short)
  long_counts <- countOverlaps(genome_bins_gr, fragments_long)
  
  # Calculate ratio, handle division by zero
  ratio <- ifelse(long_counts > 0, short_counts / long_counts, ifelse(short_counts > 0, Inf, 0))
  
  # Create a data frame for results
  bin_results <- data.frame(
    chrom = as.character(seqnames(genome_bins_gr)),
    start = start(genome_bins_gr),
    end = end(genome_bins_gr),
    ShortCount = short_counts,
    LongCount = long_counts,
    Ratio = ratio,
    Sample = sample_name,
    Group = group_name
  )
  return(bin_results)
}

get_fragment_end_motifs <- function(bam_file_path, ref_genome, motif_len = 6, valid_chroms = "chr21") {
  print(paste("Extracting fragment ends from:", basename(bam_file_path)))
  
  param <- ScanBamParam(
    what = c("qname", "flag", "rname", "pos", "strand", "qwidth", "cigar", "mrnm", "mpos"),
    flag = scanBamFlag(isPaired = TRUE, isProperPair = TRUE, isDuplicate = FALSE)
  )
  
  gal_pairs <- readGAlignmentPairs(bam_file_path, param = param)
  
  # Initial filtering for valid chromosomes
  gal_pairs <- gal_pairs[as.character(seqnames(gal_pairs)) %in% valid_chroms]
  if (length(gal_pairs) == 0) {
    message("No valid fragments found for processing in ", basename(bam_file_path), ". Returning empty DNAStringSet.")
    return(list(five_prime_motifs = DNAStringSet(), three_prime_motifs = DNAStringSet()))
  }
  
  fragments_gr <- granges(gal_pairs)
  
  # --- Vectorized creation of 5' and 3' end GRanges objects ---
  
  # 5' End GRanges
  # Start of motif is the start of the fragment
  five_prime_starts <- start(fragments_gr)
  five_prime_ends <- five_prime_starts + motif_len - 1
  
  # Create a GRanges object for all 5' end regions
  five_prime_ranges <- GRanges(
    seqnames = seqnames(fragments_gr),
    ranges = IRanges(five_prime_starts, five_prime_ends),
    strand = "*" # Strand irrelevant for fetching from genome
  )
  
  # 3' End GRanges
  # End of motif is the end of the fragment
  three_prime_ends_coord <- end(fragments_gr)
  three_prime_starts_coord <- three_prime_ends_coord - motif_len + 1
  
  # Create a GRanges object for all 3' end regions
  three_prime_ranges <- GRanges(
    seqnames = seqnames(fragments_gr),
    ranges = IRanges(three_prime_starts_coord, three_prime_ends_coord),
    strand = "*"
  )
  
  # --- Handle out-of-bounds ranges efficiently ---
  # This uses GenomicRanges' built-in capabilities to handle ranges that extend
  # beyond chromosome boundaries.
  # We'll first make sure the seqinfo in our ranges matches the ref_genome,
  # then reduce them to valid ranges.
  
  # Apply sequence information from the reference genome to the fragment ranges
  seqinfo(five_prime_ranges) <- seqinfo(ref_genome)[seqlevels(five_prime_ranges)]
  seqinfo(three_prime_ranges) <- seqinfo(ref_genome)[seqlevels(three_prime_ranges)]
  
  # Trim ranges that go beyond chromosome boundaries (this automatically handles both start and end)
  five_prime_ranges_trimmed <- trim(five_prime_ranges)
  three_prime_ranges_trimmed <- trim(three_prime_ranges)
  
  # Filter out ranges that became too short after trimming (i.e., less than motif_len)
  valid_5prime_idx <- width(five_prime_ranges_trimmed) == motif_len
  valid_3prime_idx <- width(three_prime_ranges_trimmed) == motif_len
  
  five_prime_ranges_final <- five_prime_ranges_trimmed[valid_5prime_idx]
  three_prime_ranges_final <- three_prime_ranges_trimmed[valid_3prime_idx]
  
  # --- Vectorized fetching of sequences from reference genome ---
  five_prime_motifs <- getSeq(ref_genome, five_prime_ranges_final)
  three_prime_motifs <- getSeq(ref_genome, three_prime_ranges_final)
  
  return(list(five_prime_motifs = five_prime_motifs, three_prime_motifs = three_prime_motifs))
}

filter_empty_DNAStringSets <- function(list_of_dna_stringsets) {
  filtered_list <- list()
  for (name in names(list_of_dna_stringsets)) {
    if (length(list_of_dna_stringsets[[name]]) > 0) {
      filtered_list[[name]] <- list_of_dna_stringsets[[name]]
    } else {
      message(paste("Skipping empty DNAStringSet for:", name))
    }
  }
  return(filtered_list)
}