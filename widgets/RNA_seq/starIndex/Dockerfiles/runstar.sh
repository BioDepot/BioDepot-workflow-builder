#!/bin/bash

#This wrapper around STAR does 3 things
#1 Append the outputDirectory to the outputPrefix - STAR only provides for a prefix but it is necessary to monge the directory path - so the form asks for the path separately
#2 Rearrange the fastq files for paired ends - STAR recognizes paired end reads by having 2 tokens after the readFilesIn flag. If there are multiple pairs, then the tokens are files joined by commas. If there are multiple single end reads then there is one token joined by commas
#3 Schedule multiple samples. While multiple files are possible for STAR - these are aggregated into one set of Counts. For multiple samples we need to run STAR separately. The MULTISAMPLE flag indicates that the multiple files should be parsed as multiple samples. Multiple files can also be part of a sample if they are separated by columns

#make sure everything gets killed upon exit
trap "exit" INT TERM
trap "kill 0" EXIT

function makeOutFileNamePrefix()
{
 if [[ -z ${outputDirectory} && -z ${outputPrefix} ]]; then
    outFileStr=""
 elif [[ -z ${outputDirectory} ]]; then
    outFileNamePrefix="${outputPrefix}"
    outFileStr="--outFileNamePrefix ${outFileNamePrefix} "
 else
    outFileNamePrefix="${outputDirectory}/${outputPrefix}"
    outFileStr="--outFileNamePrefix ${outFileNamePrefix} "
 fi
}    

function findOutputSubDir()
{
    local fullFile=$1
    #echo 'name for modification is ' $1
    #if it is comma separated get first element
    
    firstFile=$(echo $1 | cut -d "," -f1) 
    #echo 'first file is ' $firstFile
    
    #keep characters to right of last slash to get outputDir
    outputSubDir="${firstFile##*/}"

    #remove .gz if there 
    if [[ ${outputSubDir: -3} == '.gz' ]]; then
        outputSubDir="${outputSubDir%.*}"
    fi
    #remove .fq or fastq if there
    if [[ ${outputSubDir: -3} == '.fq' ]]; then
        outputSubDir="${outputSubDir%.*}"
    elif [[ ${outputSubDir: -6} == '.fastq' ]]; then
        outputSubDir="${outputSubDir%.*}"
    fi
    #if paired end remove the _/.R1/2 _/.1/2    
    if [[ -z "${APELIST}" && -z "${SPELIST}" ]]; then
        #single end
        echo 'single-end files do not look for R1/R1 1/2 endings' 
    else
        if [[ ${outputSubDir: -3} == '_R1' ]] || [[ ${outputSubDir: -3} == '_R2' ]]; then
            outputSubDir="${outputSubDir%_*}"
        elif [[ ${outputSubDir: -3} == '.R1' ]] || [[ ${outputSubDir: -3} == '.R2' ]]; then
            outputSubDir="${outputSubDir%.*}"
        elif [[ ${outputSubDir: -2} == '_1' ]] || [[ ${outputSubDir: -3} == '_2' ]]; then
         outputSubDir="${outputSubDir%_*}"
        elif [[ ${outputSubDir: -2} == '.1' ]] || [[ ${outputSubDir: -3} == '.2' ]]; then
         outputSubDir="${outputSubDir%.*}"
        fi  
    fi
    #add existing prefix if there is a prefix
    if [[ -z ${outFileNamePrefix} ]]; then
        return
    else
        outputSubDir="${outFileNamePrefix}${outputSubDir}/"
    fi
}   

 
#wrapper for STAR to handle multiple paired end files
makefilestr ()
{
    if [[ -z "${APELIST}" && -z "${SPELIST}" ]]; then
        filestr=$(IFS=$','; echo "${files[*]}")
    elif [ -z "${APELIST}" ]; then
        #sequential list #put commas in inputfiles        
        R1files=(${files[@]:0:$(( nFiles / 2 ))})
        R2files=(${files[@]:$(( nFiles / 2 ))})
        filestr=$(IFS=$','; echo "${R1files[*]}")' '$( IFS=$','; echo "${R2files[*]}")
    else
        R1str=""
        R2str=""    
        for ((i = 0 ; i < nFiles/2 ; i++ )); do
            R1files+=(${files[i*2]})
            R2files+=(${files[i*2+1]})
            R1str+=${files[i*2]},
            R2str+=${files[i*2+1]},
        done
        filestr+=${R1str::-1}' '${R2str::-1}
    fi
    cmdStr=${baseCmd}${outFileStr}$( IFS=' '; echo "${cmd[*]}" )' '${filestr}' '${remainingFlagsStr}
}

##Main routine starts here
#globals
outFileStr=""
baseCmd='STAR --runMode alignReads '
cmdStr=""
cmd=(${cmdStr})
remainingFlagsStr=""
remainingFlags=()
R1Files=()
R2Files=()
files=()
filestr=""
#save args
args=( "${@}" )

#will be of the form --prevFlags prevflag --inputfiles file1.. fileN --nextFlag nextflag

#APELIST SPELIST are two env variables that are passed to indicate alternating aR1 aR2 bR1 bR2 or serial aR1 aR1 bR2 bR2 list of pairs of files
#If MULTISAMPLE is not selected the lists will be converted to STAR's native format which is two tokens consisting of commas separating R1 files followed by a second token of commas separating R2 files i.e. aR1,bR1 aR2,bR2
#If neither APELIST or SPELIST are set and MULTISAMPLE is not set then the command is passed through as is - this allows for users to enter the files using the native STAR format if they so desire
#MULTISAMPLE indicates the number of samples involved
#If MULTISAMPLE is set the tokens are reordered depending on whether APELIST or SPELIST are set and pairs sent to individual STAR align commands
#If MULTISAMPLE is set and neither APELIST or SPELIST are set then tokens are sent to individual STAR align commands one at a time (single-end files)
#NOTE that tokens can still be comma joined lists of files 

makeOutFileNamePrefix 

#gather the input files
foundFiles='False'
while [ $# -gt 0 ]; do
    arg=$1
    if [ ${foundFiles} = 'False' ]; then
        #we are not in the section with the input files
        cmd+=(${arg})
        if [ ${arg} = '--readFilesIn' ]; then
            foundFiles='True'
        fi
        shift
    else
        #we are parsing now parsing filenames until we hit a flag
        if [ ${arg:0:1} = '-' ]; then
            #found flag  -no longer parsing filenames
            #process the foundfiles
            remainingFlagsStr="$*"
            remainingFlags=( "$@" )
            break
        else
            files+=(${arg})
            shift
        fi
    fi
done
#check if there is one file
nFiles=${#files[@]}
echo "number of files is ${nFiles}"
if [ ${nFiles} -le 1 ]; then
    cmdStr="${baseCmd}$outFileStr"$( IFS=' '; echo "${args[*]}" )
    echo "${cmdStr}"
    eval "${cmdStr}"
    exit 0
fi
#makefilestr is run here because it also creates the R1files and R2files array
makefilestr
if [[ -z "${MULTISAMPLE}" ]]; then
    echo 'single sample'
    echo "${cmdStr}"
    eval "${cmdStr}"
    exit 0
else
    #find and remove outFileNamePrefix from the command line
    #Now find if we are dealing with paired ends
    if [[ -z "${APELIST}" && -z "${SPELIST}" ]]; then
        #single ends - feed files one at a time
        for file in "${files[@]}"; do
            echo "working on single-end $file"
            findOutputSubDir $file
            echo "mkdir -p $outputSubDir"
            mkdir -p "${outputSubDir}"
            echo "${baseCmd} --outFileNamePrefix $outputSubDir ""$( IFS=' '; echo "${cmd[*]}" )" "$file ${remainingFlags[@]}" 
            eval "${baseCmd} --outFileNamePrefix $outputSubDir ""$( IFS=' '; echo "${cmd[*]}" )" "$file ${remainingFlags[@]}"
        done
    else
        echo "working on paired ends"
        i=0
        for R1file in "${R1files[@]}"; do
            echo "$R1file"
            R2file=${R2files[i]}
            findOutputSubDir $R1file
            echo "mkdir -p $outputSubDir"
            mkdir -p "${outputSubDir}"
            echo "${baseCmd} --outFileNamePrefix $outputSubDir ""$( IFS=' '; echo "${cmd[*]}" )" "$R1file $R2file ${remainingFlags[@]}"
            eval "${baseCmd} --outFileNamePrefix $outputSubDir ""$( IFS=' '; echo "${cmd[*]}" )" "$R1file $R2file ${remainingFlags[@]}"
            i=$((i + 1))
        done
   fi
fi
