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

# Set up the working directory to the script's location
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

# Step 1: Format the summary statistics using MungeSumstats
reformatted <- MungeSumstats::format_sumstats(
  path = input_file,
  compute_z = TRUE,       # Compute Z-scores if not present
  impute_se = TRUE,       # Impute standard errors if not present
  impute_FRQ = TRUE,      # Impute allele frequencies if not present
  frq_is_maf = FALSE      # Indicate if the frequency column is not MAF
)

# Save the reformatted summary statistics file
write.table(reformatted, formatted_file, row.names = FALSE, col.names = TRUE, sep = "\t")

# Step 2: Read the formatted file
Open_Formatted_file <- fread(formatted_file)

# Step 4: GWAS liftover to build 38 (as eQTL data is in build 38)
Gwaslifting <- liftover_gwas(
  Open_Formatted_file,
  build = c(37, 38, 36),  # Specify the builds of the input data
  to = 38,                # Target build
  chr_col = "CHR",        # Chromosome column
  pos_col = "BP",         # Position column
  snp_col = "SNP",        # SNP ID column
  ea_col = "A1",          # Effect allele column
  oa_col = "A2",          # Other allele column
  build_fallback = "position"  # Fallback method if exact match not found
)

# Step 5: Save the lifted over GWAS data to a new file
write.table(Gwaslifting, output_file, row.names = FALSE, col.names = TRUE, sep = "\t")
