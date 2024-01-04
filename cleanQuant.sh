# Isabella Moppel
# 2023/07/06


# trim off adapters with scythe - don't need to do this for paired-end reads
#for i in $(ls *_001.fastq.gz | sed -r 's/_001.fastq.gz//' | uniq)
#do
#	suffixOut="_adapt"
#	suffixIn="_001"
#	scythe -a adapters.fasta -o "${i}${suffixOut}.fastq" "${i}${suffixIn}.fastq.gz";
#done


# trim low-quality ends with sickle
for i in $(ls *_R1_001.fastq | sed -r 's/_R1_001.fastq//' | uniq)
do
	suffixR1In="_R1_001"
	suffixR2In="_R2_001"
	suffixR1Out="_R1_trimmed"
	suffixR2Out="_R2_trimmed"
	sickle pe -f "${i}${suffixR1In}.fastq" -r "${i}${suffixR2In}.fastq" -t sanger -o "${i}${suffixR1Out}.fastq" -p "${i}${suffixR2Out}.fastq" -s /dev/null;
done


# align reads with hisat2
#hisat2-build refSeq1.fna hs2_ecoli_idx # this aligns to the genome
#for i in $(ls *_R1_trimmed.fastq | sed -r 's/_R1_trimmed.fastq//' | uniq)
#do
#	suffixR1="_R1_trimmed"
#	suffixR2="_R2_trimmed"
#	hisat2 -x hs2_ecoli_idx -1 "${i}${suffixR1}.fastq" -2 "${i}${suffixR2}.fastq" -S ${i}.sam;
#done

# quantify with salmon

# ----- version that aligns with salmon: -----
# get eColi_transcriptome_fasta.fa from:
# https://ftp.ensemblgenomes.ebi.ac.uk/pub/bacteria/release-56/fasta/bacteria_0_collection/escherichia_coli_str_k_12_substr_mg1655_gca_000005845/cdna/Escherichia_coli_str_k_12_substr_mg1655_gca_000005845.ASM584v2.cdna.all.fa.gz

# build decoys - this isn't working right now!
#grep "^>" <(gunzip -c eColi_genomic_fasta.fna.gz) | cut -d " " -f 1 > decoys.txt
#sed -i.bak -e 's/>//g' decoys.txt
#cat eColi_transcriptome_fasta.fa.gz eColi_genomic_fasta.fna.gz > gentrome.fa.gz
# build index with decoys
#salmon index -t gentrome.fa.gz -d decoys.txt -i salmon_ecoli_index --gencode

# build index without decoys
salmon index -t eColi_transcriptome_fasta.fa -i salmon_ecoli_index
# quantify, mapping-based mode
for i in $(ls *_R1_trimmed.fastq | sed -r 's/_R1_trimmed.fastq//' | uniq)
do
	suffixR1="_R1_trimmed"
	suffixR2="_R2_trimmed"
	salmon quant -i salmon_ecoli_index -l IU -1 "${i}${suffixR1}.fastq" -2 "${i}${suffixR2}.fastq" --validateMappings -o quant;
done 

# ----- version that uses the alignments from HISAT2: ----- NOT CURRENTLY WORKING
#for i in $(ls *_L001.sam | sed -r 's/_L001.sam//' | uniq)
#do
#	suffix="_L001"
#	salmon quant -t refSeq1.fna -l IU -a "${i}${suffix}.sam" -o quant;
#done
