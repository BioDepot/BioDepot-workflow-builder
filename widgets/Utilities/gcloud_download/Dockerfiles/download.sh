#!/bin/bash

confirmationFile=$1
bucket=$2
outputDir=$3

project_id=$(jq -r '.project_id' $confirmationFile)
gcloud config set project $project_id
echo "gcloud auth activate-service-account --key-file=$confirmationFile"
gcloud auth activate-service-account --key-file=$confirmationFile

gsutil -q ls gs://${bucket}/ &>/dev/null
 
if [ $? != 0 ]; then
	echo "cannot find cloud directory ${bucket}"
	exit -1
fi

mkdir -p $outputDir

if [ -z $DIRS ]; then
	echo "no directories to download"
elif [ "$DIRS" == "[]" ]; then
	if [[ $(gsutil ls gs://${bucket}) ]]; then
		echo "downloading the entire bucket $bucket"
		echo "gsutil -m cp -r gs://${bucket}/* $outputDir/."
		gsutil -m cp -r gs://${bucket}/* $outputDir/.
	fi
else
	darray=( $(echo $DIRS | jq -r '.[]') )
	if [ -z $darray ]; then
		echo "cannot parse $DIRS"
	else
		for dir in "${darray[@]}"; do
			gsutil -q ls gs://${bucket}/${dir} &>/dev/null
			if [ $? == 0 ]; then
				if [[ $(gsutil ls gs://${bucket}/${dir}) ]]; then
					echo "downloading ${bucket}/${dir}"
					echo "gsutil -m cp -r gs://${bucket}/${dir}/* $outputDir/."
					gsutil -m cp -r gs://${bucket}/${dir}/* $outputDir/.
				fi
			else
				echo "cannot find or cannot access ${bucket}/${dir}"
			fi
		done
	fi
fi
if [ -z $FILES ]; then
	echo "no files to download"
else
	echo $FILES
	farray=( $(echo $FILES | jq -r '.[]') )
	for f in "${farray[@]}"; do
		gsutil -q ls gs://${bucket}/${f} &>/dev/null
		if [ $? == 0 ]; then
			echo "downloading ${bucket}/${f}"
			echo "gsutil -m cp gs://${bucket}/${f} $outputDir/."
			gsutil -m cp gs://${bucket}/${f} $outputDir/.
		else
			echo "cannot find or cannot access ${bucket}/${f}"
		fi
	done
fi
