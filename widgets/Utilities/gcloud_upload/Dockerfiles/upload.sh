#!/bin/bash

confirmationFile=$1
bucket=$2
localDir=$3
cloudDir=$4

project_id=$(jq -r '.project_id' $confirmationFile)
gcloud config set project $project_id
echo "gcloud auth activate-service-account --key-file=$confirmationFile"
gcloud auth activate-service-account --key-file=$confirmationFile

echo "gsutil ls gs://${bucket}/ &> /dev/null "
gsutil -q ls gs://${bucket} &> /dev/null
if [ $? == 0 ]; then
	if [[ -d $localDir ]]; then
		echo "gsutil -m cp -r $localDir/* gs://${bucket}/${cloudDir}/"
		gsutil -m cp -r $localDir/* gs://${bucket}/${cloudDir}/
	else
		echo "cannot find local directory $localDir"
		exit -3
	fi
	exit 0
else
	echo "cannot access bucket $bucket "
	exit -2
fi
