@echo off
echo ========================================
echo Creating Project Folder Structure
echo ========================================

REM Create data folders
mkdir data\TCGA 2>nul
mkdir data\processed 2>nul
echo Created data folders

REM Create scripts folder
mkdir scripts 2>nul
echo Created scripts folder

REM Create figures folders
mkdir figures\main_figures 2>nul
mkdir figures\supplementary_figures 2>nul
echo Created figures folders

REM Create results folder
mkdir results 2>nul
echo Created results folder

REM Create manuscript folder
mkdir manuscript 2>nul
echo Created manuscript folder

REM Create placeholder files
echo # Data Codebook > data\metadata.txt
echo # Statistical Results > results\statistical_results.csv
echo # Model Summaries > results\model_summaries.txt
echo Created placeholder files

REM Create empty R script files
echo # Download TCGA Data > scripts\01_download_tcga.R
echo # Lipid Signature Analysis > scripts\02_lipid_signature.R
echo # Immune Deconvolution > scripts\03_immune_deconvolution.R
echo # Merge Data > scripts\04_merge_data.R
echo # Correlation Analysis > scripts\05_correlation_analysis.R
echo # Survival Analysis > scripts\06_survival_analysis.R
echo # Generate Figures > scripts\07_figures.R
echo Created R script templates

echo.
echo ========================================
echo Project Structure Created Successfully!
echo ========================================
echo.
echo Next steps:
echo 1. Install R and RStudio
echo 2. Create .Rprofile in project root
echo 3. Install required packages
echo.
pause
