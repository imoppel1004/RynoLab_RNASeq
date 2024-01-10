# Jason Kuchtey & Isabella Moppel
# Ryno Lab 2023


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



# create list of quantification file names
# your working directory should be in the directory with the LMR-P1-28-0 etc folders
# can set this in the console with setwd("<dirname>")
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



# make DESEQ2 objects
dds_cold <- DESeqDataSetFromTximport(txi_cold, colData=cold_info, design=~Sugar)
dds_hot <- DESeqDataSetFromTximport(txi_hot, colData=hot_info, design=~Sugar)



#pre filtering out low count obs
filt <- rowSums(counts(dds_cold)) >= 10
dds_cold <- dds_cold[filt,]
filt <- rowSums(counts(dds_hot)) >= 10
dds_hot <- dds_hot[filt,]



# set factor level - this is detailing where the baseline is
# this is where it's important that you have the Sugar column as "treated" and "untreated"
# if you want to change that above, just make sure you change it here too
dds_cold$Sugar <- relevel(dds_cold$Sugar, ref = "untreated")
dds_hot$Sugar <- relevel(dds_hot$Sugar, ref = "untreated")



# DESeq
dds_cold <- DESeq(dds_cold)
dds_hot <- DESeq(dds_hot)



# can change the alpha value as desired
res_cold <- results(dds_cold, alpha = 0.05)
res_hot <- results(dds_hot, alpha = 0.05)



# convert EcoGene IDs to Official Gene Symbols using the gff3 file

# build the lookup_table of EcoGene IDs : Official Gene Symbols
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

# convert all cold EcoGeneIDs (b0001, b0002, etc) to Official Gene Symbols using the lookup table
b_numbers_in_res <- rownames(res_cold)
res_length <- length(b_numbers_in_res)
iter <- 1
while (iter <= res_length) {
  if (b_numbers_in_res[iter] %in% lookup_table$EcoGeneID) {
    b_number <- subset(lookup_table, EcoGeneID == b_numbers_in_res[iter], select = GeneSymbol)[1,]
    rownames(res_cold)[iter] <- b_number
  }
  print(iter)
  iter <- iter + 1
}
# and the same for hot
b_numbers_in_res <- rownames(res_hot)
res_length <- length(b_numbers_in_res)
iter <- 1
while (iter <= res_length) {
  if (b_numbers_in_res[iter] %in% lookup_table$EcoGeneID) {
    b_number <- subset(lookup_table, EcoGeneID == b_numbers_in_res[iter], select = GeneSymbol)[1,]
    rownames(res_hot)[iter] <- b_number
  }
  print(iter)
  iter <- iter + 1
}



# everything after this is OPTIONAL analysis! the hard part is over! :)



# view!
View(data.frame(res_cold))
View(data.frame(res_hot))

# print out these unfiltered data frames to file, ordered by FC
unfiltered_cold_by_FC <- res_cold[order(res_cold$log2FoldChange),]
unfiltered_hot_by_FC <- res_hot[order(res_hot$log2FoldChange),]
write.xlsx(data.frame(unfiltered_cold_by_FC), file = "unfiltered_cold_by_FC.xlsx", rowNames=TRUE)
write.xlsx(data.frame(unfiltered_hot_by_FC), file = "unfiltered_hot_by_FC.xlsx", rowNames=TRUE)


# order by smallest p-value
resOrdered_cold <- res_cold[order(res_cold$padj),]
resOrdered_hot <- res_hot[order(res_hot$padj),]
  
# how many adjusted p-values are less than 0.001?
sum(res_cold$padj < 0.001, na.rm=TRUE)
sum(res_hot$padj < 0.001, na.rm=TRUE)



# subsets with padj < 0.01 && log2FoldChange >= 1
subset_cold <- subset(data.frame(res_cold), res_cold$padj < 0.01 & (res_cold$log2FoldChange >= 1 | res_cold$log2FoldChange <= -1))
subset_cold_ordered_by_FC <- subset_cold[order(subset_cold$log2FoldChange, decreasing=TRUE),]
subset_hot <- subset(data.frame(res_hot), res_hot$padj < 0.01 & (res_hot$log2FoldChange >= 1 | res_hot$log2FoldChange <= -1))
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


# compare genes in the equivalent subsets between LB and M9
# for examples, the planktonic @28C in LB vs the planktonic @28C in M9
LB_data <- data.frame(read_excel("subset_cold_by_FC_planktonic_LB.xlsx")) # change this depending on which subset you're looking at
M9_data <- data.frame(read_excel("subset_cold_by_FC_planktonic_M9.xlsx")) # and change this to match the above but the M9 version
LB_M9_comp <- data.frame(only_LB = character(), only_LB_FC = numeric(), only_LB_padj = numeric(),
                            both = character(), both_LB_FC = numeric(), both_LB_padj = numeric(), both_M9_FC = numeric(), both_M9_padj = numeric(),
                            only_M9 = character(), only_M9_FC = numeric(), only_M9_padj = numeric(),
                            stringsAsFactors = FALSE)
LB_iter <- 1
LB_length <- length(rownames(LB_data))
while (LB_iter <= LB_length) {
  LB_name <- LB_data$...1[LB_iter]
  if (!(LB_name %in% M9_data$...1)) { # LB but not M9
    LB_M9_comp <- rbind(LB_M9_comp, list(only_LB = LB_name,
                                               only_LB_FC = LB_data$log2FoldChange[LB_iter],
                                               only_LB_padj = LB_data$padj[LB_iter],
                                               both = "", both_LB_FC = "", both_LB_padj = "", both_M9_FC = "", both_M9_padj = "",
                                               only_M9 = "", only_M9_FC = "", only_M9_padj = ""))
  }
  LB_iter <- LB_iter + 1
}  
LB_iter <- 1
LB_length <- length(rownames(LB_data))
while (LB_iter <= LB_length) {
  LB_name <- LB_data$...1[LB_iter]
  if (LB_name %in% M9_data$...1) { # LB and M9
    M9_iter <- 1
    M9_length <- length(rownames(M9_data))
    while (M9_iter <= M9_length) {
      M9_name <- M9_data$...1[M9_iter]
      if (M9_name == LB_name) {
        LB_M9_comp <- rbind(LB_M9_comp, list(only_LB = "", only_LB_FC = "", only_LB_padj = "",
                                                   both = LB_name,
                                                   both_LB_FC = LB_data$log2FoldChange[LB_iter],
                                                   both_LB_padj = LB_data$padj[LB_iter],
                                                   both_M9_FC = M9_data$log2FoldChange[M9_iter],
                                                   both_M9_padj = M9_data$padj[M9_iter],
                                                   only_M9 = "", only_M9_FC = "", only_M9_padj = ""))
      }
      M9_iter <- M9_iter + 1
    }
  }
  LB_iter <- LB_iter + 1
}
M9_iter <- 1
M9_length <- length(rownames(M9_data))
while (M9_iter <= M9_length) {
  M9_name <- M9_data$...1[M9_iter]
  if (!(M9_name %in% LB_data$...1)) { # M9 but not LB
    LB_M9_comp <- rbind(LB_M9_comp, list(only_LB = "", only_LB_FC = "", only_LB_padj = "",
                                               both = "", both_LB_FC = "", both_LB_padj = "", both_M9_FC = "", both_M9_padj = "",
                                               only_M9 = M9_name,
                                               only_M9_FC = M9_data$log2FoldChange[M9_iter],
                                               only_M9_padj = M9_data$padj[M9_iter]))
  }
  M9_iter <- M9_iter + 1
}
# view comparison
View(LB_M9_comp)
# print comparison to file
write.xlsx(LB_M9_comp, file = "LB_M9_comp_planktonic_28.xlsx", rowNames=FALSE) #change this file output name as needed


#plotting
plotMA(res_cold, ylim=c(-8,8))
plotMA(res_hot, ylim=c(-8,8))
plotCounts(dds_cold, gene=which.min(res_cold$padj), intgroup="Sugar")
plotCounts(dds_hot, gene=which.min(res_hot$padj), intgroup="Sugar")

# volcano plots
# watch the coef, make sure your Sugar column entries match the coef
resLFC <- lfcShrink(dds_cold, coef="Sugar_treated_vs_untreated", type="apeglm")
EnhancedVolcano(res_cold,
                lab = rownames(res_cold),
                x = 'log2FoldChange',
                y = 'pvalue',
                title = "28º Samples"
)
resLFC <- lfcShrink(dds_hot, coef="Sugar_treated_vs_untreated", type="apeglm")
EnhancedVolcano(res_hot,
                lab = rownames(res_hot),
                x = 'log2FoldChange',
                y = 'pvalue',
                title = "37º Samples"
)