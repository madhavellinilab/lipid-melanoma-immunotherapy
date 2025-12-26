# Install all required packages for the project
# Run this once after installing R and RStudio

cat("========================================\n")
cat("Installing Required R Packages\n")
cat("========================================\n\n")

# CRAN packages
packages_to_install <- c(
  # Data handling
  "tidyverse", "data.table", "readxl",
  # Genomics
  "fgsea", "msigdbr",
  # Survival analysis
  "survival", "survminer", "ggsurvfit", "cmprsk",
  # Statistics
  "broom", "corrr",
  # Visualization
  "ggplot2", "gridExtra", "cowplot",
  # Data manipulation
  "reshape2", "plyr"
)

# Install CRAN packages
cat("Installing CRAN packages...\n")
for (pkg in packages_to_install) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    cat(sprintf("  Installing %s...\n", pkg))
    install.packages(pkg, dependencies = TRUE)
  } else {
    cat(sprintf("  ✓ %s already installed\n", pkg))
  }
}

# Install Bioconductor
cat("\nInstalling BiocManager...\n")
if (!require("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# Install Bioconductor packages
bioc_packages <- c("TCGAbiolinks", "fgsea", "msigdbr", "immunedeconv")
cat("\nInstalling Bioconductor packages...\n")
for (pkg in bioc_packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    cat(sprintf("  Installing %s...\n", pkg))
    BiocManager::install(pkg, update = FALSE)
  } else {
    cat(sprintf("  ✓ %s already installed\n", pkg))
  }
}

cat("\n========================================\n")
cat("All packages installed successfully!\n")
cat("========================================\n")
