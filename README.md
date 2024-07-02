# MMP_GWAS_Formatting_38
# GWAS_Formatting_Pipeline

This repository contains scripts and data for formatting GWAS (Genome-Wide Association Study) summary statistics, preparing them for further analysis such as Mendelian Randomization (MR) studies.

## Overview

The goal of this project is to format GWAS summary statistics, compute necessary metrics like Z-scores and allele frequencies, and lift over genomic coordinates to a consistent build (GRCh38). This prepares the data for downstream analysis, including MR studies.

## Data Source

The input data for this project is a GWAS summary statistics file. Ensure that your GWAS data is in a tab-delimited format with columns such as chromosome (CHR), base pair position (BP), SNP ID, effect allele (A1), other allele (A2), beta, and standard error (SE).

## Prerequisites

The following R packages are required to run the scripts in this repository:

- `MungeSumstats`
- `SNPlocs.Hsapiens.dbSNP144.GRCh37`
- `BSgenome.Hsapiens.1000genomes.hs37d5`
- `SNPlocs.Hsapiens.dbSNP144.GRCh38`
- `BSgenome.Hsapiens.NCBI.GRCh38`
- `SNPlocs.Hsapiens.dbSNP155.GRCh38`
- `BSgenome.Hsapiens.UCSC.hg19.masked`
- `data.table`
- `readr`
- `tidyr`
- `dplyr`
- `ieugwasr`
- `GwasDataImport`

You can install these packages using the following command in R:

```r
install.packages(c("MungeSumstats", "SNPlocs.Hsapiens.dbSNP144.GRCh37", "BSgenome.Hsapiens.1000genomes.hs37d5", "SNPlocs.Hsapiens.dbSNP144.GRCh38", "BSgenome.Hsapiens.NCBI.GRCh38", "SNPlocs.Hsapiens.dbSNP155.GRCh38", "BSgenome.Hsapiens.UCSC.hg19.masked", "data.table", "readr", "tidyr", "dplyr", "ieugwasr", "GwasDataImport"))
```

## Directory Structure

- `scripts/`: Contains the main analysis script.
- `data/`: Directory to store input data files.
- `output/`: Directory to store output files.

## Usage

1. **Download the Data**: Ensure your GWAS summary statistics file is ready and place it in the `data/` directory.

2. **Run the Script**: Execute the R script located in the `scripts/` directory to perform the analysis.

## Step-by-Step Guide

1. **Load Necessary Libraries**: The script begins by loading the required R packages.
2. **Set Up Working Directory**: The working directory is set to the location of the script.
3. **Define Relative Paths**: The paths for input and output files are defined relative to the working directory.
4. **Create Output Directory**: The output directory is created if it doesn't exist.
5. **Format Summary Statistics**: The summary statistics file is formatted using the `MungeSumstats` package.
6. **Calculate Allele Frequency**: Allele frequencies are calculated using the `BETA` values.
7. **Liftover to Build 38**: The GWAS data is lifted over from its current build to build 38 using the `liftover_gwas` function.
8. **Save the Output**: The lifted over GWAS data is saved to a new file in the output directory.

## Example

```r
# Load necessary libraries
library(MungeSumstats)
library(SNPlocs.Hsapiens.dbSNP144.GRCh37)
library(BSgenome.Hsapiens.1000genomes.hs37d5)
library(SNPlocs.Hsapiens.dbSNP144.GRCh38)
library(BSgenome.Hsapiens.NCBI.GRCh38)
library(SNPlocs.Hsapiens.dbSNP155.GRCh38)
library(BSgenome.Hsapiens.UCSC.hg19.masked)
library(data.table)
library(readr)
library(tidyr)
library(dplyr)
library(ieugwasr)
library(GwasDataImport)

# Set up the working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Define relative paths for input and output files
input_file <- "data/fileb8493ce9bac7.tsv.gz"
output_dir <- "output"
formatted_file <- file.path(output_dir, "formatted_sumstats.tsv.gz")
output_file <- file.path(output_dir, "HF_HRC_GWAS_UKBB_EUR_Genome38.txt")

# Create the output directory if it doesn't exist
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

# Format the summary statistics
reformatted <- MungeSumstats::format_sumstats(
  path = input_file,
  compute_z = TRUE,
  impute_se = TRUE,
  impute_FRQ = TRUE,
  frq_is_maf = FALSE
)

# Save the reformatted summary statistics
write.table(reformatted, formatted_file, row.names = FALSE, col.names = TRUE, sep = "\t")

# Read the formatted file
Open_Formatted_file <- fread(formatted_file)

# Calculate the allele frequency from beta values
Open_Formatted_file$FRQ <- (exp(Open_Formatted_file$BETA)) / (1 + exp(Open_Formatted_file$BETA))

# Liftover to build 38
Gwaslifting <- liftover_gwas(
  Open_Formatted_file,
  build = c(37, 38, 36),
  to = 38,
  chr_col = "CHR",
  pos_col = "BP",
  snp_col = "SNP",
  ea_col = "A1",
  oa_col = "A2",
  build_fallback = "position"
)

# Save the lifted over GWAS data
write.table(Gwaslifting, output_file, row.names = FALSE, col.names = TRUE, sep = "\t")
```

