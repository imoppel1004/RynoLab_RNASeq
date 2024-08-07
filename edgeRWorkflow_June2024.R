# EdgeR Differential Gene Expression Analysis Workflow
# Compiled by Lisa Ryno, Jason Kuchtey, Isabella Moppel. June 2024

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
# Ryno Lab members: your working directory should be in the directory with the LMR-P1-28-0 etc folders, can set this in the console with setwd("<dirname>")
# there should be TWO of these potential working directories, one for planktonic and one for biofilm, run them separately
# your samples.txt, 28_samples.txt, 37_samples.txt should be in this directory
# these .txt files need to be filled manually by you :)
# have a Filename column (ex LMR-P1-28-0/quant/quant.sf), a Sample column
#    (P1-28), and a Sugar column (treated vs untreated) 
samples <- read.table(file="samples.txt", header = TRUE)
filenames <- paste(samples$Filename)
files <- file.path(filenames)
names(files) <- paste0(1:12) # NOTE: change this between 12 & 24 depending on how many samples you're working with
files
all(file.exists(files))



# make the file list for 28 degree samples
cold_samples <- read.table(file="28_samples.txt", header = TRUE)
cold_filenames <- paste(cold_samples$Filename)
cold_files <- file.path(cold_filenames)
names(cold_files) <- paste0(1:6) # NOTE: change this between 6 and 12 depending on how many samples you're working with

# and for the 37 degree samples
hot_samples <- read.table(file="37_samples.txt", header = TRUE)
hot_filenames <- paste(hot_samples$Filename)
hot_files <- file.path(hot_filenames)
names(hot_files) <- paste0(1:6) # NOTE: change this between 6 and 12 depending on how many samples you're working with

# create tx2gene 
# can get the gff3.gz file (renamed) from:
# https://ftp.ensemblgenomes.ebi.ac.uk/pub/bacteria/release-56/gff3/bacteria_0_collection/escherichia_coli_str_k_12_substr_mg1655_gca_000005845/Escherichia_coli_str_k_12_substr_mg1655_gca_000005845.ASM584v2.49.gff3.gz
txdb <- makeTxDbFromGFF("EColi_k12.gff3.gz")
k <- keys( txdb, keytype = "TXNAME" )
tx2gene <- select(txdb, k, "GENEID", "TXNAME") # this is "EcoGene ID, accession number"; change as you want to!
head(tx2gene)



#import quantification data to DESEQ2
#cold
txi_cold <- tximport(files = cold_files, type = "salmon", tx2gene = tx2gene)
names(txi_cold)
head(txi_cold$counts)
#hot
txi_hot <- tximport(files = hot_files, type = "salmon", tx2gene = tx2gene)
names(txi_hot)
head(txi_hot$counts)



#check comparability
all(rownames(cold_samples) %in% colnames(txi_cold$counts))
all(rownames(hot_samples) %in% colnames(txi_hot$counts))



#Get necessary columns from sample info files
# Sample is the column containing the descriptors, like P1-28, P2-37, etc
# Sugar is the column detailing whether the sample has the sugar or not.
#   Should be either "treated" or "untreated", or you'll run into problems a lil further down.
cold_info <- cold_samples[c("Sample", "Sugar")]
hot_info <- hot_samples[c("Sample", "Sugar")]

#Write the counts to a file that can be read by edgeR (note, this is done with a file txi_hot that is written to a file hotBM9counts.csv)
write.csv(txi_hot$counts, "hotBM9counts.csv")

#Note: You must change the first column header to read as "Symbol"
cts <- txi_hot$counts
normMat <- txi_hot$length
normMat <- normMat/exp(rowMeans(log(normMat)))
normCts <- cts/normMat

eff.lib <- calcNormFactors(normCts) * colSums(normCts)
normMat <- sweep(normMat, 2, eff.lib, "*")
normMat <- log(normMat)
y <- DGEList(cts)
y <- scaleOffset(y, normMat)
design <- model.matrix(~Sugar, data = hot_info)
keep <- filterByExpr(y, design)
y <- y[keep, ]

x <- read.csv("hotBM9counts.csv",row.names="Symbol")
x <- read.csv("hotBM9counts.csv", row.names = 1)
group <- factor(c("control", "treatment", "control", "treatment", "control", "treatment"))
dge <- DGEList(counts = x, group = group)
dge <- calcNormFactors(dge)
dge <- estimateDisp(dge)
fit <- glmFit(dge)
lrt <- glmLRT(fit)
topTags(lrt)
results <- topTags(lrt, n = nrow(lrt$table))
results_df <- as.data.frame(results)

gff <- read.delim("EColi_k12.gff3", header=FALSE, comment.char="#") # read in gff3 file
gff_length <- nrow(gff)
lookup_table <- data.frame(EcoGeneID = character(0), GeneSymbol = character(0)) # data frame to be filled
gff_line_no <- 1
while (gff_line_no < gff_length) {
  gff_type <- gff$V3[gff_line_no]
  gff_info <- gff$V9[gff_line_no]
  if (gff_type == "gene") {
    ecogene_id <- str_extract(gff_info, "(?<=gene_id=)[^;]+")
    gene_symbol <- str_extract(gff_info, "(?<=Name=)[^;]+")
    lookup_table <- rbind(lookup_table, data.frame(EcoGeneID = ecogene_id, GeneSymbol = gene_symbol))
  }
  gff_line_no <- gff_line_no + 1
}
b_numbers_in_res <- rownames(results_df)
res_length <- length(b_numbers_in_res)
iter <- 1
while (iter <= res_length) {
  if (b_numbers_in_res[iter] %in% lookup_table$EcoGeneID) {
    b_number <- subset(lookup_table, EcoGeneID == b_numbers_in_res[iter], select = GeneSymbol)[1,]
    rownames(results_df)[iter] <- b_number
  }
  print(iter)
  iter <- iter + 1
}
write.csv(results_df, file = "differential_expression_results_37BM9.csv")
