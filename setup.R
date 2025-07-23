# setup.R
# Script to install and load all necessary R packages for the cfDNA analysis.

# List of Bioconductor packages
bioc_packages <- c(
  "Rsamtools",
  "GenomicAlignments",
  "GenomicRanges",
  "rtracklayer",
  "GenomicFeatures",
  "TxDb.Hsapiens.UCSC.hg38.knownGene",
  "BSgenome.Hsapiens.UCSC.hg38"
)

# List of CRAN packages
cran_packages <- c(
  "ggplot2",
  "dplyr",
  "tidyr",
  "moments",
  "ggseqlogo",
  "caret",
  "pROC"
)

# Function to install and load Bioconductor packages
install_bioc_packages <- function(packages) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
  }
  for (pkg in packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      BiocManager::install(pkg)
    }
    library(pkg, character.only = TRUE)
  }
}

# Function to install and load CRAN packages
install_cran_packages <- function(packages) {
  for (pkg in packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      install.packages(pkg)
    }
    library(pkg, character.only = TRUE)
  }
}

# Run installation and loading
message("Installing and loading Bioconductor packages...")
install_bioc_packages(bioc_packages)

message("Installing and loading CRAN packages...")
install_cran_packages(cran_packages)

message("All packages installed and loaded successfully!")