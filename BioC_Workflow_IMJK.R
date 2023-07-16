# Jason Kuchtey & Isabella Moppel
# Ryno Lab 2023


# install packages
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
BiocManager::install("Matrix")
BiocManager::install('ensembldb')
BiocManager::install("AnnotationDbi")
BiocManager::install("org.EcK12.eg.db")
install.packages("readxl")



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


# create list of quantification file names
# your working directory should be in the directory with the LMR-P1-28-0 etc folders
# can set this in the console with setwd("<dirname>")
# your samples.txt, 28_samples.txt etc should be in this directory
# these .txt files need to be filled manually by you :)
# have a Filename column (ex LMR-P1-28-0/quant/quant.sf), a Sample column
#    (P1-28), and a Rhamnose column (treated vs untreated) 
samples <- read.table(file="samples.txt", header = TRUE)
filenames <- paste(samples$Filename)
files <- file.path(filenames,  "quant/quant.sf")
names(files) <- paste0(1:12)
files
all(file.exists(files))



# make the file list for 28 degree samples
cold_samples <- read.table(file="28_samples.txt", header = TRUE)
cold_filenames <- paste(cold_samples$Filename)
cold_files <- file.path(cold_filenames, "quant/quant.sf")
names(cold_files) <- paste0(1:6)

# and for the 37 degree samples
hot_samples <- read.table(file="37_samples.txt", header = TRUE)
hot_filenames <- paste(hot_samples$Filename)
hot_files <- file.path(hot_filenames, "quant/quant.sf")
names(hot_files) <- paste0(1:6)



# create tx2gene 
# can get the gff3.gz file (renamed) from:
# https://ftp.ensemblgenomes.ebi.ac.uk/pub/bacteria/release-56/gff3/bacteria_0_collection/escherichia_coli_str_k_12_substr_mg1655_gca_000005845/Escherichia_coli_str_k_12_substr_mg1655_gca_000005845.ASM584v2.49.gff3.gz
txdb <- makeTxDbFromGFF("EColi_k12.gff3.gz")
k <- keys( txdb, keytype = "TXNAME" )
tx2gene <- select(txdb, k, "GENEID", "TXNAME") # this is "EcoGene ID, accession number"; change as you want to!
head(tx2gene)


# convert all supported GENEIDs (b0001, b0002, etc) to Entrez Gene IDs (944742, 945803, etc)
acc_nos_in_txdb <- tx2gene$TXNAME
total_acc_nos_in_txdb <- length(acc_nos_in_txdb)
total_keys <- length(mappedRkeys(org.EcK12.egACCNUM2EG)) # the right keys are the accession numbers
how_many_processed <- 0
iter <- 1
while (iter < total_acc_nos_in_txdb) {
  acc_no <- acc_nos_in_txdb[iter]
  key_iter <- 1
  while (key_iter < total_keys + 1) {
    if (acc_no == mappedRkeys(org.EcK12.egACCNUM2EG)[key_iter]) {
      tx2gene$GENEID[iter] <- mappedLkeys(org.EcK12.egACCNUM2EG[key_iter])
      key_iter <- total_keys
    }
    key_iter <- key_iter + 1
    if (key_iter == total_keys + 1) {
      how_many_processed <- how_many_processed + 1
      print(how_many_processed)
    }
  }
  iter <- iter + 1
}
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



# now, convert all supported Entrez Gene IDs (944742, 945803, etc) to Official Gene Symbols - cold
entrez <- rownames(txi_cold$counts)
total_entrez <- length(entrez)
supported_entrez <- mappedkeys(org.EcK12.egSYMBOL) # gives Entrez that have a gene_name
entrez_to_gene_name_map <- as.list(org.EcK12.egSYMBOL[supported_entrez]) # maps Entrez to gene_name
total_keys <- length(mappedLkeys(org.EcK12.egSYMBOL))
how_many_processed <- 0
iter <- 1
while (iter < total_entrez) {
  e <- entrez[iter]
  print(e)
  key_iter <- 1
  while (key_iter < total_keys + 1) {
    if (e == mappedLkeys(org.EcK12.egSYMBOL)[key_iter]) {
      rownames(txi_cold$counts)[iter] <- mappedRkeys(org.EcK12.egSYMBOL[key_iter])
      key_iter <- total_keys
    }
    key_iter <- key_iter + 1
    if (key_iter == total_keys + 1) {
      how_many_processed <- how_many_processed + 1
      print(how_many_processed)
    }
  }
  iter <- iter + 1
}
head(rownames(txi_cold$counts))

# and the same but hot
entrez <- rownames(txi_hot$counts)
total_entrez <- length(entrez)
supported_entrez <- mappedkeys(org.EcK12.egSYMBOL) # gives Entrez that have a gene_name
entrez_to_gene_name_map <- as.list(org.EcK12.egSYMBOL[supported_entrez]) # maps Entrez to gene_name
total_keys <- length(mappedLkeys(org.EcK12.egSYMBOL))
how_many_processed <- 0
iter <- 1
while (iter < total_entrez) {
  e <- entrez[iter]
  print(e)
  key_iter <- 1
  while (key_iter < total_keys + 1) {
    if (e == mappedLkeys(org.EcK12.egSYMBOL)[key_iter]) {
      rownames(txi_hot$counts)[iter] <- mappedRkeys(org.EcK12.egSYMBOL[key_iter])
      key_iter <- total_keys
    }
    key_iter <- key_iter + 1
    if (key_iter == total_keys + 1) {
      how_many_processed <- how_many_processed + 1
      print(how_many_processed)
    }
  }
  iter <- iter + 1
}
head(rownames(txi_hot$counts))



#check comparability
all(rownames(cold_samples) %in% colnames(txi_cold$counts))
all(rownames(hot_samples) %in% colnames(txi_hot$counts))


#Get necessary columns from sample info files
# Sample is the column containing the descriptors, like P1-28, P2-37, etc
# Rhamnose is the column detailing whether the sample has rhamnose or not.
#   Should be either "treated" or "untreated", or you'll run into problems a lil further down.
cold_info <- cold_samples[c("Sample", "Rhamnose")]
hot_info <- hot_samples[c("Sample", "Rhamnose")]


# make DESEQ2 objects
dds_cold <- DESeqDataSetFromTximport(txi_cold, colData=cold_info, design=~Rhamnose)
dds_hot <- DESeqDataSetFromTximport(txi_hot, colData=hot_info, design=~Rhamnose)


#pre filtering out low count obs
filt <- rowSums(counts(dds_cold)) >= 10
dds_cold <- dds_cold[filt,]
filt <- rowSums(counts(dds_hot)) >= 10
dds_hot <- dds_hot[filt,]


# set factor level - this is detailing where the baseline is
# this is where it's important that you have the Rhamnose column as "treated" and "untreated"
# if you want to change that above, just make sure you change it here too
dds_cold$Rhamnose <- relevel(dds_cold$Rhamnose, ref = "untreated")
dds_hot$Rhamnose <- relevel(dds_hot$Rhamnose, ref = "untreated")


# DESeq
dds_cold <- DESeq(dds_cold)
dds_hot <- DESeq(dds_hot)

# can change the alpha value as desired
res_cold <- results(dds_cold, alpha = 0.05)
res_hot <- results(dds_hot, alpha = 0.05)



# everything after this is OPTIONAL analysis! the hard part is over! :)



# view!
View(data.frame(res_cold))
View(data.frame(res_hot))



# order by smallest p-value
resOrdered_cold <- res_cold[order(res_cold$padj),]
resOrdered_hot <- res_hot[order(res_hot$padj),]
  
# how many adjusted p-values are less than 0.001?
sum(res_cold$padj < 0.001, na.rm=TRUE)
sum(res_hot$padj < 0.001, na.rm=TRUE)



# subsets with padj < 0.01 && log2FoldChange >= 1
subset_cold <- subset(data.frame(res_cold), res_cold$padj < 0.01 & res_cold$log2FoldChange >= 1)
subset_cold_ordered_by_FC <- subset_cold[order(subset_cold$log2FoldChange, decreasing=TRUE),]
subset_hot <- subset(data.frame(res_hot), res_hot$padj < 0.01 & res_hot$log2FoldChange >= 1)
subset_hot_ordered_by_FC <- subset_hot[order(subset_hot$log2FoldChange, decreasing=TRUE),]

# print the above subsets to file
write.xlsx(subset_cold_ordered_by_FC, file = "subset_cold_by_FC.xlsx", rowNames=TRUE)
write.xlsx(subset_hot_ordered_by_FC, file = "subset_hot_by_FC.xlsx", rowNames=TRUE)



# compare genes in the hot and the cold subsets for either planktonic or biofilm
# (depending on where your working directory is currently, move working directory as needed)
hot_cold_comp <- data.frame(only_28 = character(), only_28_FC = numeric(), only_28_padj = numeric(),
                            both = character(), both_28_FC = numeric(), both_28_padj = numeric(), both_37_FC = numeric(), both_37_padj = numeric(),
                            only_37 = character(), only_37_FC = numeric(), only_37_padj = numeric(),
                            stringsAsFactors = FALSE)
cold_iter <- 1
cold_length <- length(rownames(data.frame(subset_cold_ordered_by_FC)))
while (cold_iter <= cold_length) {
  cold_name <- rownames(data.frame(subset_cold_ordered_by_FC))[cold_iter]
  if (!(cold_name %in% rownames(data.frame(subset_hot)))) { # cold but not hot
    hot_cold_comp <- rbind(hot_cold_comp, list(only_28 = cold_name,
                                               only_28_FC = subset_cold_ordered_by_FC$log2FoldChange[cold_iter],
                                               only_28_padj = subset_cold_ordered_by_FC$padj[cold_iter],
                                               both = "", both_28_FC = "", both_28_padj = "", both_37_FC = "", both_37_padj = "",
                                               only_37 = "", only_37_FC = "", only_37_padj = ""))
  }
  cold_iter <- cold_iter + 1
}  
cold_iter <- 1
cold_length <- length(rownames(data.frame(subset_cold_ordered_by_FC)))
while (cold_iter <= cold_length) {
  cold_name <- rownames(data.frame(subset_cold_ordered_by_FC))[cold_iter]
  if (cold_name %in% rownames(data.frame(subset_hot))) { # cold and hot
    hot_iter <- 1
    hot_length <- length(rownames(data.frame(subset_hot_ordered_by_FC)))
    while (hot_iter <= hot_length) {
      hot_name <- rownames(data.frame(subset_hot_ordered_by_FC))[hot_iter]
      if (hot_name == cold_name) {
        hot_cold_comp <- rbind(hot_cold_comp, list(only_28 = "", only_28_FC = "", only_28_padj = "",
                                                   both = cold_name,
                                                   both_28_FC = subset_cold_ordered_by_FC$log2FoldChange[cold_iter],
                                                   both_28_padj = subset_cold_ordered_by_FC$padj[cold_iter],
                                                   both_37_FC = subset_hot_ordered_by_FC$log2FoldChange[hot_iter],
                                                   both_37_padj = subset_hot_ordered_by_FC$padj[hot_iter],
                                                   only_37 = "", only_37_FC = "", only_37_padj = ""))
      }
      hot_iter <- hot_iter + 1
    }
  }
  cold_iter <- cold_iter + 1
}
hot_iter <- 1
hot_length <- length(rownames(data.frame(subset_hot_ordered_by_FC)))
while (hot_iter <= hot_length) {
  hot_name <- rownames(data.frame(subset_hot_ordered_by_FC))[hot_iter]
  if (!(hot_name %in% rownames(data.frame(subset_cold)))) { # hot but not cold
    hot_cold_comp <- rbind(hot_cold_comp, list(only_28 = "", only_28_FC = "", only_28_padj = "",
                                               both = "", both_28_FC = "", both_28_padj = "", both_37_FC = "", both_37_padj = "",
                                               only_37 = hot_name,
                                               only_37_FC = subset_hot_ordered_by_FC$log2FoldChange[hot_iter],
                                               only_37_padj = subset_hot_ordered_by_FC$padj[hot_iter]))
  }
  hot_iter <- hot_iter + 1
}
# view comparison
View(hot_cold_comp)
# print comparison to file
write.xlsx(hot_cold_comp, file = "hot_cold_comp.xlsx", rowNames=FALSE)



# now, compare the planktonic and biofilm subsets for the cold genes; and the hot genes
# first, STOP!! and please move the "subset_cold_by_FC.xlsx" and "subset_hot_by_FC.xlsx"
#   files into your working directory, and make sure to rename them with "_planktonic"
#   or "_biofilm" on the end (as an example, one of your finished files should be
#   "subset_cold_by_FC_planktonic.xlsx"). then you may proceed :)
# cold genes first
cold_planktonic <- read_excel("subset_cold_by_FC_planktonic.xlsx")
colnames(cold_planktonic)[colnames(cold_planktonic) == "...1"] <- "names" # renaming after that pesky auto-rename
cold_biofilm <- read_excel("subset_cold_by_FC_biofilm.xlsx")
colnames(cold_biofilm)[colnames(cold_biofilm) == "...1"] <- "names" # renaming after that pesky auto-rename
planktonic_vs_biofilm_comp_cold <- data.frame(only_planktonic = character(), only_planktonic_FC = numeric(), only_planktonic_padj = numeric(),
                                              both = character(), both_planktonic_FC = numeric(), both_planktonic_padj = numeric(), both_biofilm_FC = numeric(), both_biofilm_padj = numeric(),
                                              only_biofilm = character(), only_biofilm_FC = numeric(), only_biofilm_padj = numeric(),
                                              stringsAsFactors = FALSE)
plank_iter <- 1
plank_length <- length(cold_planktonic$names)
while (plank_iter <= plank_length) {
  plank_name <- cold_planktonic$names[plank_iter]
  if (!(plank_name %in% cold_biofilm$names)) { # planktonic but not biofilm
    planktonic_vs_biofilm_comp_cold <- rbind(planktonic_vs_biofilm_comp_cold, list(
      only_planktonic = plank_name, only_planktonic_FC = cold_planktonic$log2FoldChange[plank_iter], only_planktonic_padj = cold_planktonic$padj[plank_iter],
      both = "", both_planktonic_FC = "", both_planktonic_padj = "", both_biofilm_FC = "", both_biofilm_padj = "",
      only_biofilm = "", only_biofilm_FC = "", only_biofilm_padj = ""))
  }
  plank_iter <- plank_iter + 1
}  
plank_iter <- 1
plank_length <- length(cold_planktonic$names)
while (plank_iter <= plank_length) {
  plank_name <- cold_planktonic$names[plank_iter]
  if (plank_name %in% cold_biofilm$names) { # planktonic and biofilm
    biofilm_iter <- 1
    biofilm_length <- length(cold_biofilm$names)
    while (biofilm_iter <= biofilm_length) {
      biofilm_name <- cold_biofilm$names[biofilm_iter]
      if (plank_name == biofilm_name) { # both planktonic and biofilm
        planktonic_vs_biofilm_comp_cold <- rbind(planktonic_vs_biofilm_comp_cold, list(
          only_planktonic = "", only_planktonic_FC = "", only_planktonic_padj = "",
          both = plank_name, both_planktonic_FC = cold_planktonic$log2FoldChange[plank_iter], both_planktonic_padj = cold_planktonic$padj[plank_iter],
          both_biofilm_FC = cold_biofilm$log2FoldChange[biofilm_iter], both_biofilm_padj = cold_biofilm$padj[biofilm_iter],
          only_biofilm = "", only_biofilm_FC = "", only_biofilm_padj = ""))
      }
      biofilm_iter <- biofilm_iter + 1
    }
  }
  plank_iter <- plank_iter + 1
}
biofilm_iter <- 1
biofilm_length <- length(cold_biofilm$names)
while (biofilm_iter <= biofilm_length) {
  biofilm_name <- biofilm_name <- cold_biofilm$names[biofilm_iter]
  if (!(biofilm_name %in% cold_planktonic$names)) { # biofilm but not planktonic
    planktonic_vs_biofilm_comp_cold <- rbind(planktonic_vs_biofilm_comp_cold, list(
      only_planktonic = "", only_planktonic_FC = "", only_planktonic_padj = "",
      both = "", both_planktonic_FC = "", both_planktonic_padj = "", both_biofilm_FC = "", both_biofilm_padj = "",
      only_biofilm = biofilm_name, only_biofilm_FC = cold_biofilm$log2FoldChange[biofilm_iter], only_biofilm_padj = cold_biofilm$padj[biofilm_iter]))
  }
  biofilm_iter <- biofilm_iter + 1
}
# view comparison
View(planktonic_vs_biofilm_comp_cold)
# print comparison to file
write.xlsx(planktonic_vs_biofilm_comp_cold, file = "planktonic_vs_biofilm_comp_cold.xlsx", rowNames=FALSE)



# next, compare the planktonic and biofilm subsets for the hot genes
hot_planktonic <- read_excel("subset_hot_by_FC_planktonic.xlsx")
colnames(hot_planktonic)[colnames(hot_planktonic) == "...1"] <- "names" # renaming after that pesky auto-rename
hot_biofilm <- read_excel("subset_hot_by_FC_biofilm.xlsx")
colnames(hot_biofilm)[colnames(hot_biofilm) == "...1"] <- "names" # renaming after that pesky auto-rename
planktonic_vs_biofilm_comp_hot <- data.frame(only_planktonic = character(), only_planktonic_FC = numeric(), only_planktonic_padj = numeric(),
                                              both = character(), both_planktonic_FC = numeric(), both_planktonic_padj = numeric(), both_biofilm_FC = numeric(), both_biofilm_padj = numeric(),
                                              only_biofilm = character(), only_biofilm_FC = numeric(), only_biofilm_padj = numeric(),
                                              stringsAsFactors = FALSE)
plank_iter <- 1
plank_length <- length(hot_planktonic$names)
while (plank_iter <= plank_length) {
  plank_name <- hot_planktonic$names[plank_iter]
  if (!(plank_name %in% hot_biofilm$names)) { # planktonic but not biofilm
    planktonic_vs_biofilm_comp_hot <- rbind(planktonic_vs_biofilm_comp_hot, list(
      only_planktonic = plank_name, only_planktonic_FC = hot_planktonic$log2FoldChange[plank_iter], only_planktonic_padj = hot_planktonic$padj[plank_iter],
      both = "", both_planktonic_FC = "", both_planktonic_padj = "", both_biofilm_FC = "", both_biofilm_padj = "",
      only_biofilm = "", only_biofilm_FC = "", only_biofilm_padj = ""))
  }
  plank_iter <- plank_iter + 1
}  
plank_iter <- 1
plank_length <- length(hot_planktonic$names)
while (plank_iter <= plank_length) {
  plank_name <- hot_planktonic$names[plank_iter]
  if (plank_name %in% hot_biofilm$names) { # planktonic and biofilm
    biofilm_iter <- 1
    biofilm_length <- length(hot_biofilm$names)
    while (biofilm_iter <= biofilm_length) {
      biofilm_name <- hot_biofilm$names[biofilm_iter]
      if (plank_name == biofilm_name) { # both planktonic and biofilm
        planktonic_vs_biofilm_comp_hot <- rbind(planktonic_vs_biofilm_comp_hot, list(
          only_planktonic = "", only_planktonic_FC = "", only_planktonic_padj = "",
          both = plank_name, both_planktonic_FC = hot_planktonic$log2FoldChange[plank_iter], both_planktonic_padj = hot_planktonic$padj[plank_iter],
          both_biofilm_FC = hot_biofilm$log2FoldChange[biofilm_iter], both_biofilm_padj = hot_biofilm$padj[biofilm_iter],
          only_biofilm = "", only_biofilm_FC = "", only_biofilm_padj = ""))
      }
      biofilm_iter <- biofilm_iter + 1
    }
  }
  plank_iter <- plank_iter + 1
}
biofilm_iter <- 1
biofilm_length <- length(hot_biofilm$names)
while (biofilm_iter <= biofilm_length) {
  biofilm_name <- biofilm_name <- hot_biofilm$names[biofilm_iter]
  if (!(biofilm_name %in% hot_planktonic$names)) { # biofilm but not planktonic
    planktonic_vs_biofilm_comp_hot <- rbind(planktonic_vs_biofilm_comp_hot, list(
      only_planktonic = "", only_planktonic_FC = "", only_planktonic_padj = "",
      both = "", both_planktonic_FC = "", both_planktonic_padj = "", both_biofilm_FC = "", both_biofilm_padj = "",
      only_biofilm = biofilm_name, only_biofilm_FC = hot_biofilm$log2FoldChange[biofilm_iter], only_biofilm_padj = hot_biofilm$padj[biofilm_iter]))
  }
  biofilm_iter <- biofilm_iter + 1
}
# view comparison
View(planktonic_vs_biofilm_comp_hot)
# print comparison to file
write.xlsx(planktonic_vs_biofilm_comp_hot, file = "planktonic_vs_biofilm_comp_hot.xlsx", rowNames=FALSE)



#plotting
plotMA(res_cold, ylim=c(-8,8))
plotMA(res_hot, ylim=c(-8,8))
plotCounts(dds_cold, gene=which.min(res_cold$padj), intgroup="Rhamnose")
plotCounts(dds_hot, gene=which.min(res_hot$padj), intgroup="Rhamnose")

# volcano plots
# watch the coef, make sure your Rhamnose column entries match the coef
resLFC <- lfcShrink(dds_cold, coef="Rhamnose_treated_vs_untreated", type="apeglm")
EnhancedVolcano(res_cold,
                lab = rownames(res_cold),
                x = 'log2FoldChange',
                y = 'pvalue',
                title = "28º Samples"
)
resLFC <- lfcShrink(dds_hot, coef="Rhamnose_treated_vs_untreated", type="apeglm")
EnhancedVolcano(res_hot,
                lab = rownames(res_hot),
                x = 'log2FoldChange',
                y = 'pvalue',
                title = "37º Samples"
)