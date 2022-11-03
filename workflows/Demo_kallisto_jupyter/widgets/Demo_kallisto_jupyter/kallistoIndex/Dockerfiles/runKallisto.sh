#!/bin/bash
#wrap kallisto to check for multiple samples

#make sure everything gets killed upon exit
trap "exit" INT TERM
trap "kill 0" EXIT

function getName(){
	local fullFile=$1
	#keep characters to right of last slash to get filename
	filename="${fullFile##*/}"
	#remove .gz if there 
	if [ ${filename: -3} == '.gz' ]; then
		filename="${filename%.*}"
	fi
	#remove .fq or fastq if there
	if [ ${filename: -3} == '.fq' ]; then
		filename="${filename%.*}"
	elif [ ${filename: -6} == '.fastq' ]; then
		filename="${filename%.*}"
	fi
	#if paired end remove the _/.R1/2 _/.1/2	
	if [ -z ${singleFlag} ]; then
		if [[ ${filename: -3} == '_R1' ]] || [[ ${filename: -3} == '_R2' ]]; then
			filename="${filename%_*}"
		elif [[ ${filename: -3} == '.R1' ]] || [[ ${filename: -3} == '.R2' ]]; then
			filename="${filename%.*}"
		elif [[ ${filename: -2} == '_1' ]] || [[ ${filename: -3} == '_2' ]]; then
		 filename="${filename%_*}"
		elif [[ ${filename: -2} == '.1' ]] || [[ ${filename: -3} == '.2' ]]; then
		 filename="${filename%.*}"
		fi	
 fi
}	

if [ -z ${MULTISAMPLE} ]; then
	echo kallisto quant $@
	kallisto quant $@
else
	#default is paired end reads - check if the single read flag is set 
	nargs=$((MULTISAMPLE + MULTISAMPLE))
	params=("$@")
	
	for ((i = 0; i < ${#params[@]}; i++ )); do
	 if [ "${params[$i]}" == '--single' ]; then
			nargs=$MULTISAMPLE
			singleFlag=1
			break
	 fi
	done
	
	#find the -o parameter and remove it
	for ((i = 0; i < ${#params[@]}; i++ )); do
	 if [ "${params[$i]}" == '-o' ]; then
	  j=$((i + 1))
	  outputDir=${params[j]}
	  params=(${params[@]:0:$i} ${params[@]:$((j + 1))})
	  break
	 fi
	done
	
	nparms=${#params[@]}
	cmdSize=$((nparms - nargs))
	cmdSlice=("${params[@]:0:$cmdSize}")
 
 args=()
 outputPaths=()
 #get the arguments and paths
	for (( i = $cmdSize ; i < $nparms; i++ )); do
  if [ -z ${singleFlag} ]; then
	   file1=${params[i]}
	   i=$((i + 1))
	   args+=("$file1 ${params[i]}")
	   getName $file1
	 else
	 	 args+=(${params[i]})
	 	 getName ${params[i]}
	 fi
	 outputPaths+=("${outputDir}/${filename}")
	done
	#execute the commands
 for ((i = 0; i < $MULTISAMPLE; i++ )); do
  #echo $filename
  echo mkdir -p ${outputPaths[i]}
  mkdir -p ${outputPaths[i]}
		echo kallisto quant ${cmdSlice[@]} -o ${outputPaths[i]} ${args[i]}
		kallisto quant ${cmdSlice[@]} -o ${outputPaths[i]} ${args[i]}
 done	
fi
