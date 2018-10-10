#!/bin/bash

#Env variables passed from widget

#descriptionFile: table with at least a sample name column and a condition/category column for identifying the baseline and treatment conditions for differential expression eg.
# sample	condition	path
#SRR493366	scramble	/data/SRR493366
#SRR493367	scramble	/data/SRR493367
#SRR493368	scramble	/data/SRR493368
#SRR493369	HOXA1KD	/data/SRR493369
#SRR493370	HOXA1KD	/data/SRR493370
#SRR493371	HOXA1KD	/data/SRR493371
#
#countsFile: tab delimited file with table with gene names as rows names, counts as rows and sample names as columns eg.
#             SRR493366 SRR493368 SRR493367 SRR493369 SRR493370 SRR493371
#ENST0000001       209       164       143       162        80       151
#ENST0000002       150       102        45        99        68       250
#    .              .         .         .         .         .         .
#control: The value of the condition entry for the baseline samples to compare the treatment samples against
#treatment: The value of the condition entry for the treatment samples
#nGenes: The number of top genes to display - not currently used but could be used to generate smaller results file
#outputDir: 

#countsFile='deseqCountsTable.tsv'
#descriptionFile='exp_desc.tsv'
#control='scramble'
#treatment='HOXA1KD'
#nGenes=50
#outputDir='/data'
#nLines=$((nGenes+1))

function createRscript {
	echo '#!/usr/bin/env Rscript'
	echo 'library(DESeq2)'
	echo 'descriptionTable <- read.table('"'""${descriptionFile}""'"', header = TRUE, stringsAsFactors = FALSE)'
	echo 'sampleNames <- descriptionTable$sample'
	echo 'condition <- descriptionTable$condition'
	echo 'dataPath<-'"'""${countsFile}""'"
	echo "countData = read.table(file = dataPath, header = TRUE, row.names = 1, sep = '\t')"
	echo 'colData <- data.frame(row.names=colnames(countData), condition=factor(condition,levels=c('"'""${control}""','""${treatment}""')))"
	echo 'dataset <- DESeqDataSetFromMatrix(countData = countData,colData = colData,design = ~condition)'
	echo 'DESeq(dataset)'
	echo 'dds <-DESeq(dataset)'
	echo "result <- results(dds, contrast=c('condition','""${treatment}','""${control}'))"
	echo 'result <- result[complete.cases(result),]'
	echo 'resOrdered <- result[order(result$padj),]'
	echo "write.table(resOrdered,'${outputDir}/results.tsv',sep='\t',row.names=TRUE)"
}
createRscript > runRscript.sh
chmod +x runRscript.sh
./runRscript.sh
nLines=$((nGenes+1))
echo "sed '1s/^/Gene\t/' ${outputDir}/results.tsv | head -${nLines} | sed 's/\\"'"//g'"' >${outputDir}/top${nGenes}Genes.tsv"
eval "sed '1s/^/Gene\t/' ${outputDir}/results.tsv | head -${nLines} | sed 's/\\"'"//g'"' >${outputDir}/top${nGenes}Genes.tsv"
