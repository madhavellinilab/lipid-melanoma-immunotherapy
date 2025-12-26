# R Project-Specific Settings
# Loaded automatically when opening this project

# Set working directory
setwd(getwd())

# Create data directories (if not exist)
dir.create("data/TCGA", showWarnings = FALSE, recursive = TRUE)
dir.create("data/processed", showWarnings = FALSE, recursive = TRUE)
dir.create("figures/main_figures", showWarnings = FALSE, recursive = TRUE)
dir.create("results", showWarnings = FALSE, recursive = TRUE)

# Load common packages
suppressWarnings(library(tidyverse))
suppressWarnings(library(data.table))

# Custom function for saving figures publication-ready
save_figure <- function(plot, filename, width = 8, height = 6, dpi = 300) {
  ggsave(filename, plot, width = width, height = height, dpi = dpi)
  cat(sprintf("Saved: %s\n", filename))
}

# Print startup message
cat("\n=== Lipid Metabolism Project Loaded ===\n")
cat("Working directory:", getwd(), "\n")
cat("Ready for analysis!\n\n")
