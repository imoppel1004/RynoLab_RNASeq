# EdgeR Differential Gene Expression Analysis Workflow
#Compiled by Lisa Ryno, Isabella Moppel, and Jason Kuchtey June 2024

# install BiocManager
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# change the library path as needed, to the first value returned by typing
#   .libPaths()
#   in the console
my_lib_path = .libPaths()[1]

# install packages
BiocManager::install(c('Rsubread'), lib=my_lib_path)
BiocManager::install(c('DESeq2'), lib=my_lib_path)
BiocManager::install(c('tximport'), lib=my_lib_path)
BiocManager::install(c("rhdf5"), lib=my_lib_path)
BiocManager::install(c("BSgenome"), lib=my_lib_path)
BiocManager::install(c("GenomicFeatures"), lib=my_lib_path)
BiocManager::install(c('apeglm'), lib=my_lib_path)
BiocManager::install(c('EnhancedVolcano'), lib=my_lib_path)
BiocManager::install(c("Matrix"), lib=my_lib_path)
BiocManager::install(c('ensembldb'), lib=my_lib_path)
BiocManager::install(c("AnnotationDbi"), lib=my_lib_path)
BiocManager::install(c("org.EcK12.eg.db"), lib=my_lib_path)
install.packages(c("readxl"), lib=my_lib_path)
install.packages(c("openxlsx"), lib=my_lib_path)
install.packages(c("stringr"), lib=my_lib_path)
install.packages(c("edgeR"), lib=my_lib_path)


# load packages
library(GenomicFeatures)
library(Rsubread)
library(DESeq2)
library(tximport)
library(rhdf5)
library(EnhancedVolcano)
library(apeglm)
library(Matrix)
library(ensembldb)
library(AnnotationDbi)
library(org.EcK12.eg.db)
library(readxl)
library(openxlsx)
library(stringr)
library(edgeR)

# create list of quantification file names
# your working directory should be in the directory with the LMR-P1-28-0 etc folders
# can set this in the console with setwd("<dirname>")
# these .txt files need to be filled manually by you :)
# have a Filename column (ex LMR-P1-28-0/quant/quant.sf), a Sample column
#    (P1-28), and a Sugar column (treated vs untreated) 
samples <- read.table(file="samples.txt", header = TRUE)
filenames <- paste(samples$Filename)
files <- file.path(filenames)
names(files) <- paste0(1:48) # NOTE: change this between 12, 24 & 48 depending on how many samples you're working with
files
all(file.exists(files))

# create tx2gene 
# can get the gff3.gz file (renamed) from:
# https://ftp.ensemblgenomes.ebi.ac.uk/pub/bacteria/release-56/gff3/bacteria_0_collection/escherichia_coli_str_k_12_substr_mg1655_gca_000005845/Escherichia_coli_str_k_12_substr_mg1655_gca_000005845.ASM584v2.49.gff3.gz
txdb <- makeTxDbFromGFF("EColi_k12.gff3.gz")
k <- keys( txdb, keytype = "TXNAME" )
tx2gene <- select(txdb, k, "GENEID", "TXNAME") # this is "EcoGene ID, accession number"; change as you want to!
head(tx2gene)

#import quantification data 

txi_all <- tximport(files = files, type = "salmon", tx2gene = tx2gene)
names(txi_all)
head(txi_all$counts)


#check comparability

all(rownames(samples) %in% colnames(txi_all$counts))


#Get necessary columns from sample info files
# Sample is the column containing the descriptors, like P1-28, P2-37, etc
# Temp refers to the temperatures, either 28 or 37
# Media refers to M9 or LB
# Sugar is the column detailing whether the sample has the sugar or not.
#   Should be either "treated" or "untreated", or you'll run into problems further down.
all_info <- samples[c("Sample", "BorP", "Temp", "Media", "Sugar")]

#Write the counts to a file that can be read by edgeR
write.csv(txi_all$counts, "allcounts.csv")

# Check comparability
if (!all(rownames(samples) %in% colnames(txi_all$counts))) {
  stop("Sample row names do not match counts column names.")
}
print("Sample row names match counts column names.")

# Assign counts and normalization matrix
cts <- txi_all$counts
normMat <- txi_all$length
normMat <- normMat / exp(rowMeans(log(normMat)))

# Normalize counts
normCts <- cts / normMat
eff.lib <- calcNormFactors(normCts) * colSums(normCts)
normMat <- sweep(normMat, 2, eff.lib, "*")
normMat <- log(normMat)

# Create DGEList
y <- DGEList(counts = cts)

# Assign sample information
y$samples <- all_info

# Check dimensions
print(paste("Dimensions of cts:", dim(cts)))
print(paste("Dimensions of y$samples:", dim(y$samples)))
if (ncol(cts) != nrow(all_info)) {
  stop("Number of columns in counts does not match number of rows in sample info.")
}
print("Number of columns in counts matches number of rows in sample info.")

# Check structure of y$samples
print(str(y$samples))

# Explicitly calculate and set library sizes
y$samples$lib.size <- colSums(y$counts)
print(paste("Library sizes after setting:", y$samples$lib.size))

# Normalize factors
y <- calcNormFactors(y)
print(paste("Normalized factors:", y$samples$norm.factors))

# Calculate normalized counts (log-transformed)
norm_counts <- cpm(y, log = TRUE)

# Check dimensions of normalized counts
print(paste("Dimensions of norm_counts:", dim(norm_counts)))

# Perform PCA
pca <- prcomp(t(norm_counts), scale = TRUE)

# Print PCA summary
summary(pca)

# Get the principal component scores
pc_scores <- as.data.frame(pca$x)

# Get the explained variance ratios
explained_variance <- as.data.frame(pca$sdev^2 / sum(pca$sdev^2))

# Save the principal component scores to a CSV file
write.csv(pc_scores, "PCA_scores.csv", row.names = TRUE)

# Save the explained variance ratios to a CSV file
write.csv(explained_variance, "PCA_explained_variance.csv", row.names = FALSE)
