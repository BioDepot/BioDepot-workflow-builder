#!/bin/bash
#saves name and counts column depending on the type of read
#The parameters passed by bwb are 3 environment variables and 2 sets of of command line parameters
#The env variable inputFile defines the STAR filename (not including path)
#The env variable outputFile defines the DESEQ output filename
#The env variable column defines the column used for the counts and depends on whether the reads are single or paired ends
#The -i flag is followed by all the directories where the inputFiles are found usually i.e dir1/ReadsPerGene.out.tab dir2/ReadsPerGene.out.tab
#if there is no outputDir then we use the same directory as the inputDirectory to output the file 

#
#gather the input files
trap "exit" INT TERM
trap "kill 0" EXIT
#inputFile='ReadsPerGene.out.tab'
#outputFile='deseqCountsTable.tsv'
column=3
function getList()
{
    list=()
    foundFlag='False'
    for arg in "${args[@]}"; do
        if [ ${arg} = $1 ]; then
            foundFlag='True'
        elif [ ${arg:0:1} = '-' ]; then
            if [ ${foundFlag} = 'True' ]; then
                break
            fi
        elif [ ${foundFlag} = 'True' ]; then
            list+=(${arg})
        fi 
    done
}    
function pasteColumns()
{
	#echo "tail -n +5 $1/${inputFile} | cut -d$'\t' -f""${column} | paste -d $'\t' ${outputFile} - > ${outputFile}_tmp"
 eval "tail -n +5 $1/${inputFile} | cut -d$'\t' -f""${column} | paste -d $'\t' ${outputFile} - > ${outputFile}_tmp"
 mv "${outputFile}"_tmp "${outputFile}"
}

args=($@)
list=()
getList '-i'
inputDirs=("${list[@]}")

#get the gene list from first file
tail -n +5 "${inputDirs[0]}/${inputFile}" | cut -d$'\t' -f1 > "${outputFile}"

#join columns from counts files
for inputDir in "${inputDirs[@]}"; do 
    pasteColumns "${inputDir}"
done

#add header
header=$( IFS=$'\t'; echo "${inputDirs[*]}" )
sed -i "1i${header}" "${outputFile}"

