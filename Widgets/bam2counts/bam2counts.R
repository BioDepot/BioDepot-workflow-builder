
cmd_args <- commandArgs(TRUE)
if(length(cmd_args) < 3){
    cat(paste("Usage: Rscript bam2counts.R [GTF File] [Sample Table File] [BAM file Directory] [Workers number]\n"))
    cat("\n")
    cat("Arguments:\n\n")
    cat("\t[GTF File] a file specifying the genomic features.\n\n")
    cat("\t[Sample Table File] is a table with detailed information for each of samples that links samples to the associated FASTQ and BAM files.\n\n")
    cat("\t[BAM file Directory] is a directory contains BAM files.\n\n")
    cat("\t[Workers number] is a parameter defined number of threads, default = 1")
    stop()
}

start.time <- Sys.time()

gtf_file <- cmd_args[1]
sample_table <- cmd_args[2]
bam_filepath <- cmd_args[3]
workers_number <- 1
if(length(cmd_args)>=4){
    workers_number <- cmd_args[4]
}

print(paste0("Working threads: ", workers_number))

print("Loading sample table....")
sampleTable <- read.csv(sample_table, row.names = 1)

print("Loading GTF file....")
library("GenomicFeatures")
txdb <- makeTxDbFromGFF(gtf_file, format = "gtf", circ_seqs = character())
ebg <- exonsBy(txdb, by="gene")
print("GTF file loaded.")

print("Scanning BAM files......")
filenames <-list.files(bam_filepath, '*.bam')
bamfilenames <- unlist(lapply(filenames, function(x) {file.path(bam_filepath,x)}))
print(paste0("BAM files: ", filenames))

library("Rsamtools")
bamfiles <- BamFileList(bamfilenames)

library("GenomicAlignments")
library("BiocParallel")

if(workers_number == 1){
    register(SerialParam())
}else{
    register(MulticoreParam(workers = workers_number))
}

print("Summerize overlaps......")
se <- summarizeOverlaps(features=ebg, reads=bamfiles,
                        mode="Union",
                        singleEnd=FALSE,
                        ignore.strand=TRUE,
                        fragments=TRUE )


se_file <- file.path(bam_filepath, "se.RDATA")
print(paste0("Saving se data to: ", se_file))
save(se, file=se_file)

library("DESeq2")
print("DESeq......")
colData(se) <- DataFrame(sampleTable)

dds <- DESeqDataSet(se, design = ~ cell + dex)

dds_file <- file.path(bam_filepath, "dds.RDATA")

print(paste0("Saving dds to: ", dds_file))

save(dds, file=dds_file)

result_file <- file.path(bam_filepath, "deseq_results.csv")
print(paste0("Generating and saving result table to: ", result_file))
dds <- DESeq(dds)
res <- results(dds)
resOrdered <- res[order(res$pvalue),]
write.csv(as.data.frame(resOrdered), file=result_file)

end.time <- Sys.time()
time.taken <- end.time - start.time

print(paste0("Time used: ", time.taken))
