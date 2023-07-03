# Isabella Moppel
# 2023/07/03

# requires that index is already built, named ref_index
# otherwise, build index from refSeq.fna by uncommenting the below
#salmon index -t refSeq.fna -i ref_index

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
for i in $(ls *_R1_trimmed.fastq | sed -r 's/_R1_trimmed.fastq//' | uniq)
do
	suffixR1="_R1_trimmed"
	suffixR2="_R2_trimmed"
	hisat2 -x ref_index -1 "${i}${suffixR1}.fastq" -2 "${i}${suffixR2}.fastq" -S ${i}.sam;
done

# quantify with salmon
for i in $(ls *_R1_trimmed.fastq | sed -r 's/_R1_trimmed.fastq//' | uniq)
do
	suffixR1="_R1_trimmed"
	suffixR2="_R2_trimmed"
	salmon quant -i ref_index -l A -1 "${i}${suffixR1}.fastq" -2 "${i}${suffixR2}.fastq" -o quant;
done 
