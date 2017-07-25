
cmd_args <- commandArgs(TRUE)
if(length(cmd_args) < 2){
    cat(paste("Usage: Rscript bam2counts.R [GTF File] [BAM file Directory] [Workers number]\n"))
    cat("\n")
    cat("Arguments:\n\n")
    cat("\t[GTF File] a file specifying the genomic features.\n\n")
    cat("\t[BAM file Directory] is a directory contains BAM files.\n\n")
    cat("\t[Workers number] is a parameter defined number of threads, default = 1")
    stop()
}

start.time <- Sys.time()

gtf_file <- cmd_args[1]
bam_filepath <- cmd_args[2]
workers_number <- 1
if(length(cmd_args)>=3){
    workers_number <- cmd_args[3]
}

print(paste0("Working threads: ", workers_number))

print(paste0("Loading GTF file....", gtf_file))
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

bam_counts = assay(se)
counts_file <- file.path(bam_filepath, "bamcounts.csv")
print(paste0("Saving counts table to: ", counts_file))
write.table(bam_counts, file=counts_file, sep = ",")


end.time <- Sys.time()
time.taken <- end.time - start.time

print(paste0("Time used: ", time.taken))
