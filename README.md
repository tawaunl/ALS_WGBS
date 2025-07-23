# ALS WGBS

This repository contains R scripts for analyzing cell-free DNA (cfDNA) samples, specifically focusing on identifying features that can distinguish between ALS (Amyotrophic Lateral Sclerosis) and Control samples. The analysis involves insert size distribution, fragment end motif analysis, fragment coverage over genomic features, and a classification model to predict disease status.

The provided data is a downsampled subset of chr21 from a larger dataset, processed for reduced file size to facilitate quick analysis.

## Project Structure

```
.
├── Analysis/
│   ├── functions.R
│   ├── InitialAnalysis.R (Provided - not directly run, its logic is integrated into run_analysis.R)
│   ├── run_analysis.R
│   └── setup.R
├── celfie/ (Data directory - to be placed at the same level as Analysis)
│   ├── bam/
│   │   ├── 
```
The celfie data folder can be downloaded from the data provided [here.](https://drive.google.com/file/d/1Xt4SPbL8sYa9WAnUP_j0rZ_T9AJl7oTr/view)
The data is a subset of data taken from [this paper](https://www.nature.com/articles/s41467-021-22901-x#Sec11) and contains 12 ALS samples and 10 Control samples. Data is processed and downsampled, and only selected for chr21 to reduce file size. 

## Setup Instructions

To run this analysis, you need R installed. The necessary R packages can be installed using the `setup.R` script.

1.  **Clone the Repository (or download files):**
    Ensure you have the `Analysis/` folder containing `functions.R`, `run.R`, and `setup.R`.

2.  **Place Data:**
    Create a `celfie/` directory at the same level as your `Analysis/` directory. Inside `celfie/`, place:
    * `bam/` directory containing all `.bam` files.
    * `celfie_cfDNA_ss.csv` (metadata file).
    * `cpgIslandExt.bed` (CpG island annotations).

    Your directory structure should look like this:

    ```
    .
    ├── Analysis/
    │   ├── functions.R
    │   ├── run_analysis.R
    │   └── setup.R
    └── celfie/
        ├── bam/
        ├── celfie_cfDNA_ss.csv
        └── cpgIslandExt.bed
    ```

3.  **Install R Packages:**
    Open an R console or RStudio, navigate to the `Analysis/` directory (or set your working directory to `Analysis/`), and run the setup script:

    ```R
    setwd("/path/to/your/Analysis/directory")
    source("setup.R")
    ```
    This script will install all required Bioconductor and CRAN packages. This step may take some time depending on your internet connection and if packages are already installed.

## Running the Analysis

After setting up the environment and data, you can run the main analysis script:

1.  **Open R/RStudio:**
    Ensure your working directory is set to the `Analysis/` folder.

2.  **Execute the Analysis Script:**

    ```R
    source("run_analysis.R")
    ```

The `run_analysis.R` script will:
* Read metadata and BAM file paths.
* Perform insert size analysis, generating distribution plots and statistics.
* Analyze mono-nucleosomal and di-nucleosomal fragment proportions.
* Calculate statistical moments (skewness, kurtosis) of insert sizes.
* Analyze short fragment enrichment.
* Extract fragments overlapping genomic features (CpG islands, promoters) and analyze their insert sizes.
* Calculate short/long fragment ratios across genomic bins on chr21.
* Extract and compare fragment end motifs (5' and 3' ends) between groups, generating sequence logos.
* Analyze fragment end distributions around genomic features.
* Perform binary classification (ALS vs. Control) using Logistic Regression and Random Forest models based on derived features, and report classification metrics (precision, sensitivity, F1-score, AUC-ROC).

## Outputs

All generated plots (PDF and PNG) and summary CSV files will be saved in the `Figures/` directory.

Key output files include:
* `InsertSizeDensity_bySample.pdf`: Histogram of insert sizes per sample.
* `InsertSizeDensity_byGroup.pdf`: Density plot of insert sizes by group.
* `MononucleosomalProportionCLEAN_boxplot.pdf`: Boxplot of mono-nucleosomal fragment proportions.
* `Mono_Di_Ratio_boxplot.pdf`: Boxplot of mono/di-nucleosomal fragment ratios.
* `fragment_length_summary_statistics_chr21.csv`: CSV with mean, median, SD, skewness, and kurtosis of insert sizes.
* `FragmentLengthMoments_boxplots.pdf`: Boxplots of statistical moments.
* `ShortInsertSizeDensity_byGroup.pdf`: Density plot of short insert sizes by group.
* `ShortFragmentEnrichment_boxplot.pdf`: Boxplot of short fragment proportions.
* `InsertSizeDistribution_Features_chr21.pdf`: Insert size distributions within CpG islands and promoters.
* `ShortLongRatio_Chr21_Bins.pdf`: Short/long fragment ratio across chr21 in 1Mb bins.
* `Motif_5primeDE.csv` & `Motif_3primeDE.csv`: Differentially enriched fragment end motifs.
* `top_*_motifs_6bp.png`: Sequence logo plots for top enriched motifs.
* `fragment_end_distributions_chr21.png`: Fragment end position distributions around features.
* `fragment_end_distributions_significant_bins_chr21.png`: Same as above, with statistically significant bins highlighted.
* `roc_curve_glm.pdf`, `roc_curve_rf.pdf`, `roc_curves_combined.pdf`: ROC curves for classification models.
* `rf_feature_importance.png`: Feature importance plot from Random Forest.
* `classification_performance_metrics.csv`: CSV summarizing classification performance metrics.

## Scalability Notes

The provided data is downsampled to chr21 for faster execution. For full genomic datasets (e.g., 50GB BAM files), consider:
* **Parallel Processing:** Utilizing `BiocParallel` for `scanBam` or `readGAlignmentPairs` operations across multiple BAM files or even within single large BAM files.
* **Streaming/Chunking:** Processing BAM files in chunks rather than loading the entire file into memory (though `Rsamtools` functions are often optimized for this).
* **High-Performance Computing (HPC):** Running analyses on clusters with sufficient RAM and CPU cores.
* **Database Integration:** For extremely large datasets, consider loading processed data into a genomic database for querying.

## Contact

For any questions or issues, please refer to the R script comments.

---

# Analysis Logic and Implementation Details for cfDNA Data

This document provides a high-level overview of the analytical steps performed on the cell-free DNA (cfDNA) data to identify potential markers differentiating ALS from Control samples. The analysis leverages various aspects of cfDNA fragmentomics to uncover distinguishing patterns.

## 1. Data Loading and Setup

The analysis begins by organizing the input data. This involves identifying the paths to the BAM files (which contain the aligned sequencing reads) and loading a metadata file that categorizes each sample as either ALS or Control. This initial step ensures that samples are correctly grouped for comparative analysis.

## 2. Insert Size Analysis

Insert size, or fragment length, is a crucial characteristic of cfDNA. This section explores how the lengths of DNA fragments differ between the ALS and Control groups.

* **Distribution Analysis:** The distribution of fragment lengths is examined for all samples, both individually and aggregated by disease group. This helps visualize overall patterns.
* **Nucleosomal Fragment Proportions:** Specific attention is paid to fragments corresponding to mono-nucleosomal (around 150-180 base pairs) and di-nucleosomal (around 300-380 base pairs) lengths. The proportions of these fragment types are calculated and statistically compared between groups to identify potential differences related to nucleosome packaging.
* **Statistical Moments:** Beyond simple averages, statistical measures like skewness and kurtosis are calculated for fragment length distributions. These metrics describe the shape of the distribution, which can reveal subtle but significant differences between the groups.
* **Short Fragment Enrichment:** The proportion of very short fragments (e.g., 85-105 base pairs) is specifically investigated, as these can be indicative of particular biological processes or disease states.

## 3. Fragment Coverage and Feature Overlap Analysis

This part of the analysis investigates how cfDNA fragments are distributed across specific functional regions of the genome, such as CpG islands and gene promoters on chromosome 21. This can provide insights into chromatin accessibility and gene regulation.

* **Fragment Overlap with Features:** Fragments that overlap with known genomic features are identified, and their insert size distributions within these regions are analyzed and compared between ALS and Control samples.
* **Regional Fragment Ratios:** The ratio of short fragments to long fragments is calculated across large genomic bins (e.g., 1 megabase) along chromosome 21. This helps to detect larger-scale regional variations in fragmentation patterns that might be associated with disease.

## 4. Fragment End Motif Analysis

The specific DNA sequences at the very ends of cfDNA fragments can offer clues about the enzymes that cut the DNA in the body. This section aims to identify and compare these "fragment end motifs."

* **Motif Extraction:** Short DNA sequences (e.g., 6 base pairs) at both the 5' and 3' ends of the fragments are extracted for all samples.
* **Differential Abundance:** The frequencies of these motifs are compared between the ALS and Control groups. Statistical tests are performed to identify motifs that are significantly more or less abundant in one group compared to the other.
* **Sequence Logos:** Visual representations called "sequence logos" are generated for the most differentially abundant motifs, highlighting conserved nucleotides and providing insights into potential cleavage preferences.

## 5. Fragment Start and End Position Distributions around Features

This analysis delves into the precise positioning of fragment ends relative to key genomic landmarks, such as transcription start sites (TSS) of genes and the centers of CpG islands. This can reveal patterns related to nucleosome positioning or transcription factor binding.

* **Relative Positioning:** The distances of fragment ends (both 5' and 3') from the centers of specified genomic features are calculated.
* **Distribution Comparison:** The distributions of these relative positions are plotted and compared between the ALS and Control groups. Regions where fragment ends are significantly enriched or depleted are identified.
* **Statistical Significance:** Statistical tests are applied to determine if differences in fragment end distributions within specific genomic "bins" around features are statistically significant between the disease and control groups. Significant bins are highlighted to draw attention to regions of interest.

## 6. Binary Classification (ALS vs. Control)

The ultimate goal is to determine if the various cfDNA features identified in the previous steps can be used to accurately distinguish between ALS and Control samples. This section builds and evaluates predictive models.

* **Feature Compilation:** Key statistical summaries derived from the insert size analysis (e.g., mean, median, skewness, kurtosis, nucleosomal proportions, short fragment enrichment) are combined into a single dataset.
* **Model Training:** Two common machine learning models, Logistic Regression and Random Forest, are trained using this compiled feature set. Repeated cross-validation is used to ensure robust model evaluation.
* **Performance Evaluation:** The models' ability to classify samples is assessed using standard metrics such as:
    * **ROC AUC (Area Under the Receiver Operating Characteristic Curve):** A measure of overall discriminatory power.
    * **Sensitivity:** The model's ability to correctly identify ALS samples.
    * **Specificity:** The model's ability to correctly identify Control samples.
    * **Precision:** The accuracy of positive (ALS) predictions.
    * **F1-Score:** A balance between precision and sensitivity.
* **ROC Curves and Feature Importance:** ROC curves are plotted to visualize the trade-off between sensitivity and specificity, and feature importance is extracted from the Random Forest model to highlight which cfDNA characteristics are most influential in predicting disease status.

## 7. Methylation Analysis (Planned/Future Work)

The provided data is bisulfite-treated, which allows for the detection of DNA methylation. While not explicitly implemented in the current scripts, extracting and analyzing methylation patterns from the `XM` tag in the BAM files is identified as a valuable area for future research. This would involve identifying methylated cytosines and comparing their distribution and frequency between ALS and Control samples.
Additionally, DMA analysis and  methylation patterens at fragment ends.

## 8. Expanded Genomic Context Analysis
The current analysis looks at promoters and CpG islands. You could broaden this to other functional genomic elements because they were the easiest to extract.
Future work can look at Enhancers an other regulatory elements, repeat elemetns and maybe disease specific regions.

## 9. Expanded nuclesome Profiling
We can expand on the profiling done to look at phasing as well as nuclesome free regions.

## 10. Machine Learning Model Enhancements
instead of using basic features we can create more complex features that may better than the current ones used. This will increase our F1- scores
we can also explore more Advanced models like SVMs and Neural networks.

## Conclusion

This comprehensive analysis pipeline provides a robust framework for exploring cfDNA fragmentation patterns and their potential as biomarkers for distinguishing disease states. The current implementation focuses on various aspects of fragment length, genomic distribution, and end motifs, culminating in a classification module. The results, including statistical summaries, plots, and classification metrics, are output to the `Figures` directory, offering valuable insights into the differences between ALS and Control cfDNA profiles. The modular design, with separate functions and analysis scripts, ensures maintainability and scalability for future investigations with larger datasets or additional analyses like methylation.