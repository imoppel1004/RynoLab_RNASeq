# Isabella Moppel
# 2023/06/29

#Install Packages
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

# align paired files for 2022 samples
# 28-0
align(index = "eColiRefIndex",
      readfile1 = "LMR-P1-28-0_S31_L001_R1_001.fastq",
      readfile2 = "LMR-P1-28-0_S31_L001_R2_001.fastq",
      output_format = "SAM")
align(index = "eColiRefIndex",
      readfile1 = "LMR-P2-28-0_S32_L001_R1_001.fastq",
      readfile2 = "LMR-P2-28-0_S32_L001_R2_001.fastq",
      output_format = "SAM")
align(index = "eColiRefIndex",
      readfile1 = "LMR-P3-28-0_S33_L001_R1_001.fastq",
      readfile2 = "LMR-P3-28-0_S33_L001_R2_001.fastq",
      output_format = "SAM")
# 28-5
align(index = "eColiRefIndex",
      readfile1 = "LMR-P1-28-5_S34_L001_R1_001.fastq",
      readfile2 = "LMR-P1-28-5_S34_L001_R2_001.fastq",
      output_format = "SAM")
align(index = "eColiRefIndex",
      readfile1 = "LMR-P2-28-5_S35_L001_R1_001.fastq",
      readfile2 = "LMR-P2-28-5_S35_L001_R2_001.fastq",
      output_format = "SAM")
align(index = "eColiRefIndex",
      readfile1 = "LMR-P3-28-5_S36_L001_R1_001.fastq",
      readfile2 = "LMR-P3-28-5_S36_L001_R2_001.fastq",
      output_format = "SAM")
# 37-0
align(index = "eColiRefIndex",
      readfile1 = "LMR-P1-37-0_S37_L001_R1_001.fastq",
      readfile2 = "LMR-P1-37-0_S37_L001_R2_001.fastq",
      output_format = "SAM")
align(index = "eColiRefIndex",
      readfile1 = "LMR-P2-37-0_S38_L001_R1_001.fastq",
      readfile2 = "LMR-P2-37-0_S38_L001_R2_001.fastq",
      output_format = "SAM")
align(index = "eColiRefIndex",
      readfile1 = "LMR-P3-37-0_S39_L001_R1_001.fastq",
      readfile2 = "LMR-P3-37-0_S39_L001_R2_001.fastq",
      output_format = "SAM")
# 37-5
align(index = "eColiRefIndex",
      readfile1 = "LMR-P1-37-5_S40_L001_R1_001.fastq",
      readfile2 = "LMR-P1-37-5_S40_L001_R2_001.fastq",
      output_format = "SAM")
align(index = "eColiRefIndex",
      readfile1 = "LMR-P2-37-5_S41_L001_R1_001.fastq",
      readfile2 = "LMR-P2-37-5_S41_L001_R2_001.fastq",
      output_format = "SAM")
align(index = "eColiRefIndex",
      readfile1 = "LMR-P3-37-5_S42_L001_R1_001.fastq",
      readfile2 = "LMR-P3-37-5_S42_L001_R2_001.fastq",
      output_format = "SAM")