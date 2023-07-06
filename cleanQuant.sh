# Isabella Moppel
# 2023/07/06


# trim off adapters with scythe
for i in $(ls *_001.fastq.gz | sed -r 's/_001.fastq.gz//' | uniq)
do
	suffixOut="_adapt"
	suffixIn="_001"
	scythe -a adapters.fasta -o "${i}${suffixOut}.fastq" "${i}${suffixIn}.fastq.gz";
done


# trim low-quality ends with sickle
for i in $(ls *_R1_adapt.fastq | sed -r 's/_R1_adapt.fastq//' | uniq)
do
	suffixR1In="_R1_adapt"
	suffixR2In="_R2_adapt"
	suffixR1Out="_R1_trimmed"
	suffixR2Out="_R2_trimmed"
	sickle pe -f "${i}${suffixR1In}.fastq" -r "${i}${suffixR2In}.fastq" -t sanger -o "${i}${suffixR1Out}.fastq" -p "${i}${suffixR2Out}.fastq" -s /dev/null;
done


# align reads with hisat2
hisat2-build refSeq1.fna hs2_ecoli_idx
for i in $(ls *_R1_trimmed.fastq | sed -r 's/_R1_trimmed.fastq//' | uniq)
do
	suffixR1="_R1_trimmed"
	suffixR2="_R2_trimmed"
	hisat2 -x hs2_ecoli_idx -1 "${i}${suffixR1}.fastq" -2 "${i}${suffixR2}.fastq" -S ${i}.sam;
done

# quantify with salmon
# get EColi_k12_cdna.all.fa from:
# https://ftp.ensemblgenomes.ebi.ac.uk/pub/bacteria/release-56/fasta/bacteria_0_collection/escherichia_coli_str_k_12_substr_mg1655_gca_000005845/cdna/Escherichia_coli_str_k_12_substr_mg1655_gca_000005845.ASM584v2.cdna.all.fa.gz
salmon index -t EColi_k12_cdna.all.fa -i salmon_ecoli_index
for i in $(ls *_R1_trimmed.fastq | sed -r 's/_R1_trimmed.fastq//' | uniq)
do
	suffixR1="_R1_trimmed"
	suffixR2="_R2_trimmed"
	salmon quant -i salmon_ecoli_index -l IU -1 "${i}${suffixR1}.fastq" -2 "${i}${suffixR2}.fastq" --validateMappings -o quant;
done 
