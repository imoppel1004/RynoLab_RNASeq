# Isabella Moppel
# 2023/06/30

# Create new install location (used when running in hpc)
#dir.create(Sys.getenv("R_LIBS_USER"), recursive = TRUE)
#.libPaths(Sys.getenv("R_LIBS_USER"))

# Set CRAN mirror (used when running from terminal)
#r = getOption("repos")
#r["CRAN"] = "http://cran.us.r-project.org"
#options(repos = r)

# Install packages
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install('Rsubread')
BiocManager::install("DESeq2")
BiocManager::install("tximport")
BiocManager::install("rhdf5")
BiocManager::install("BSgenome")
BiocManager::install("GenomicFeatures")
BiocManager::install('apeglm')
BiocManager::install('EnhancedVolcano')
BiocManager::install('ensembldb')


library(GenomicFeatures)
library(Rsubread)
library(DESeq2)
library(tximport)
library(rhdf5)
library(EnhancedVolcano)
library(apeglm)
library(ensembldb)

# build index
buildindex(basename = "eColiRefIndex",
           reference = "refSeq.fna",
           gappedIndex = FALSE,
           memory = 8000,
           TH_subread = 100,
           colorspace = FALSE)

# align paired files
# 28-0
align(index = "eColiRefIndex",
      readfile1 = "LMR-B1-28-0_S19_L002_R1_001.fastq",
      readfile2 = "LMR-B1-28-0_S19_L002_R2_001.fastq",
      output_format = "SAM")
align(index = "eColiRefIndex",
      readfile1 = "LMR-B2-28-0_S20_L002_R1_001.fastq",
      readfile2 = "LMR-B2-28-0_S20_L002_R2_001.fastq",
      output_format = "SAM")
align(index = "eColiRefIndex",
      readfile1 = "LMR-B3-28-0_S21_L002_R1_001.fastq",
      readfile2 = "LMR-B3-28-0_S21_L002_R2_001.fastq",
      output_format = "SAM")
# 28-5
align(index = "eColiRefIndex",
      readfile1 = "LMR-B1-28-5_S22_L002_R1_001.fastq",
      readfile2 = "LMR-B1-28-5_S22_L002_R2_001.fastq",
      output_format = "SAM")
align(index = "eColiRefIndex",
      readfile1 = "LMR-B2-28-5_S23_L002_R1_001.fastq",
      readfile2 = "LMR-B2-28-5_S23_L002_R2_001.fastq",
      output_format = "SAM")
align(index = "eColiRefIndex",
      readfile1 = "LMR-B3-28-5_S24_L002_R1_001.fastq",
      readfile2 = "LMR-B3-28-5_S24_L002_R2_001.fastq",
      output_format = "SAM")
# 37-0
align(index = "eColiRefIndex",
      readfile1 = "LMR-B1-37-0_S25_L002_R1_001.fastq",
      readfile2 = "LMR-B1-37-0_S25_L002_R2_001.fastq",
      output_format = "SAM")
align(index = "eColiRefIndex",
      readfile1 = "LMR-B2-37-0_S26_L002_R1_001.fastq",
      readfile2 = "LMR-B2-37-0_S26_L002_R2_001.fastq",
      output_format = "SAM")
align(index = "eColiRefIndex",
      readfile1 = "LMR-B3-37-0_S27_L002_R1_001.fastq",
      readfile2 = "LMR-B3-37-0_S27_L002_R2_001.fastq",
      output_format = "SAM")
# 37-5
align(index = "eColiRefIndex",
      readfile1 = "LMR-B1-37-5_S28_L002_R1_001.fastq",
      readfile2 = "LMR-B1-37-5_S28_L002_R2_001.fastq",
      output_format = "SAM")
align(index = "eColiRefIndex",
      readfile1 = "LMR-B2-37-5_S29_L002_R1_001.fastq",
      readfile2 = "LMR-B2-37-5_S29_L002_R2_001.fastq",
      output_format = "SAM")
align(index = "eColiRefIndex",
      readfile1 = "LMR-B3-37-5_S30_L002_R1_001.fastq",
      readfile2 = "LMR-B3-37-5_S30_L002_R2_001.fastq",
      output_format = "SAM")