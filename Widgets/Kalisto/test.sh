kallisto index -i transcripts.idx transcripts.fasta.gz
kallisto quant -i transcripts.idx -o output -b 100 reads_1.fastq.gz reads_2.fastq.gz
cd output
cat abundance.tsv
