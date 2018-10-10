#!/bin/bash
#Ling-Hong Hung 2018 lhhunghimself@gmail.com
#Use as you want - attribute if you can

#cli shell wrapper that generates an R script for sleuth based on the options shown
#separate sleuth prep from the analyses functions which require models

#parameters are
# descriptionFile csv dataframe -d
# full_model text -f
#read_bootstrap_tpm binary -r 
#extra_bootstrap_summary    -e 
#transformation_function default log(x+0.5) -t 
#num_cores default -n 1

#sleuth_prep(sample_to_covariates, full_model = NULL,
#  filter_fun = basic_filter, target_mapping = NULL, max_bootstrap = NULL,
#  norm_fun_counts = norm_factors, norm_fun_tpm = norm_factors,
#  aggregation_column = NULL, read_bootstrap_tpm = FALSE,
#  extra_bootstrap_summary = FALSE, transformation_function = log_transform,
#  num_cores = max(1L, parallel::detectCores() - 1L), ...)

# first create the script

declare -A prepParams

function createRscript {

	#find file s2c file
	for ((i = 0; i < ${#params[@]}; i++ )); do
	 if [ "${params[$i]}" == '-d' ]; then
	  j=$((i + 1))
	  filePath=${params[j]}
	  params=(${params[@]:0:$i} ${params[@]:$((j + 1))})
	  break
	 fi
	done
	
	echo '#!/usr/bin/env Rscript'
	echo 'library(sleuth)'
	echo 's2c <- read.table('"'""${filePath}""'"', header = TRUE, stringsAsFactors = FALSE)'
	
	
	newParams=()
	
	for ((i = 0; i < ${#params[@]}; i++ )); do
	 arg=${params[i]}
	 #string args
	 case $arg in
	  "--full_model" | "--filter_fun" | "--norm_fun_counts" | "--norm_fun_tpm" | "--transformation_function" | "--aggregation_column")
	   #strip the --
	   arg="${arg:2:${#arg}-2}"
	   i=$((i + 1))
	   nextArg=${params[i]}
	   nextArg=$(sed -e 's/^"//' -e 's/"$//' <<<"$nextArg")
	   prepParams[$arg]="'""$nextArg""'"
	  ;;
	 #integer args
	  "--num_cores")
	   arg="${arg:2:${#arg}-2}"
	   i=$((i + 1))
	   nextArg=${params[i]}
	   prepParams[$arg]=$nextArg
	  ;;
	 #boolean parameters
	  "--read_bootstrap_tpm" | "--extra_bootstrap_summary")
	   arg="${arg:2:${#arg}-2}"
	   prepParams[$arg]='TRUE'
	  ;;
	 #column argument
	  "--column" )
	   i=$((i + 1))
	   nextArg=${params[i]}
	   nextArg=$(sed -e 's/^"//' -e 's/"$//' <<<"$nextArg")
	   columnArg="~$nextArg"
	  ;;
	  #use lrt argument
	  "--wald" )
	   i=$((i + 1))
	   nextArg=${params[i]}
	   nextArg=$(sed -e 's/^"//' -e 's/"$//' <<<"$nextArg")
	   which_beta=$nextArg
	  ;;
	  #top N genes 
	  "--nGenes" )
	   i=$((i + 1))
	   nextArg=${params[i]}
	   nGenes=$nextArg
	  ;;
	  #output_file 
	  "--output_file" )
	   i=$((i + 1))
	   nextArg=${params[i]}
	   nextArg=$(sed -e 's/^"//' -e 's/"$//' <<<"$nextArg")
	   outputFile=$nextArg
	  ;;
	  #gene_names
	  "--geneNamesFile" )
	   i=$((i + 1))
	   geneNamesFile=${params[i]}
	  ;;
	  "--qvalue" )
	   i=$((i + 1))
	   nextArg=${params[i]}
	   qval=$nextArg
	  ;;    
	 esac
	done
	if [ -n "$geneNamesFile" ]; then
	 echo 't2g<- read.table('"'""${geneNamesFile}""'"',header = TRUE, stringsAsFactors = FALSE)'
	 if [ -n "prepParams[$arg]" ]; then
	  prepCmd="so <- sleuth_prep(s2c, ${columnArg} , target_mapping = t2g "
	 else
	  prepCmd="so <- sleuth_prep(s2c, ${columnArg} , target_mapping = t2g, aggregation_column =""'"'ext_gene'"'"
	 fi
	else
	 prepCmd="so <- sleuth_prep(s2c, ${columnArg} "
	fi
	for i in "${!prepParams[@]}"
	do
	  prepCmd+=",$i=${prepParams[$i]}"
	done
	prepCmd+=")"
	echo $prepCmd
	
	fitCmdFull="so <- sleuth_fit(so, ${columnArg},""'"'full'"'"')'
	fitCmdReduced='so <- sleuth_fit(so, ~1,'"'"'reduced'"'"')'
	if [ -z "$which_beta" ]; then
	 #lrt model
	 echo "so <- sleuth_fit(so, ${columnArg},""'"full"'"')'
	 echo 'so <- sleuth_fit(so, ~1,'"'"'reduced'"'"')'
	 echo 'so <- sleuth_lrt(so, '"'"'reduced'"'"', '"'"'full'"')"
	 echo 'full_results <- sleuth_results(so, '"'"'reduced:full'"'"','"'"'lrt'"'"', show_all = FALSE)'
	else
	 #wald model
	 echo 'so <- sleuth_wt(so, '"'"$which_beta"'"')'
	 echo 'sleuth_significant <- sleuth_results(so, '"'"$which_beta"'"')'
	fi 
	#now get the table
	echo "sleuth_significant <- dplyr::filter(full_results, qval <= $qval)"
	if [ -n "$nGenes" ]; then
	 echo 'write.table(sleuth_significant[1:'"${nGenes}"',],file = '"'""$outputFile""'"', sep="\t",row.names=FALSE)'
	else
	 echo 'write.table(sleuth_significant,file = '"'""$outputFile""'"', sep="\t",row.names=FALSE)'
	fi
}

params=("$@")
createRscript > ./runSleuth.R
chmod +x ./runSleuth.R
bash -c './runSleuth.R'



