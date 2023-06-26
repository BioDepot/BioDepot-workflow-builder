#!/bin/bash

#things to do - can obtain manifest from endpoint for any guid to get file and md5 hash so we can verify before proceeding

gdcapiurl='https://api.gdc.cancer.gov'
exitCode=0

#create a temporary directory for API download and mv the downloads when complete an no error
#remove the temporary directory if there is a interrupt or error
#if there is an error save the uid for download using gen3
#need to also check for gen3 errors

function cleanup(){
	#use global tempDir as this will depend on the GUID being downloaded or the manifest being downloaded
	#make sure that basename tempdir starts with temp so that we don't remove something disastrous with a typo
	basetemp=$(basename $tempDir)
	[[ -n "$basetemp" && "${basetemp:0:4}" = "temp" ]] && rm -rf $tempDir
}

#call cleanup if we exit not so nicely
trap "cleanup  -1 " SIGINT INT TERM

function getFilename(){
	local tempDir="$(mktemp -d)"
	#make a temporary directory without write permissions to force curl to quit after obtaining filename
	chmod -w $tempDir
	pushd $tempDir > /dev/null
	local filename=$(su user -c "curl -OJ --header \"X-Auth-Token: $token\" $gdcapiurl/data/$1 \
		|& grep Warning | sed 's/.* \(.*\):.*/\1/' | grep -v Warning")
	popd > /dev/null
	rmdir $tempDir
	echo "$filename"
}

function decompressFile(){
	local filename="$1"
	case $filename in
	*.tar.xz | *.txz)
		tar -xJf "$filename"
		;;
	*.tar.bz2 | *.tbz2)
		tar -xjf "$filename"
		;;
	*.tar.gz | *.tgz)
		tar -xzf "$filename"
		;;
	*.tar)
		tar -xf "$filename"
		;;
	*.gz)
		gunzip "$filename"
		;;
	*.bz2)
		bunzip2 "$filename"
		;;
	*.zip)
		unzip "$filename"
		rm "$filename"
		;;
	*)
		echo "Decompress not supported for $filename"
		return 1
		;;
	esac
	return $?
}

#there is a bug with the multi-download api in that the files are not properly tarballed - so use only singleDownload endpoint
# returns 1 if successful otherwise 0
function singleDownloadGET(){
	local myid=$1
	local downloaded=0
	local filename=$(getFilename $myid)
	if [ -n "$filename" ]; then
		echo "Filename is $filename"
	else
		echo "ERROR: failed to retrieve filename"
		exit 1
	fi
	# check if we can skip download for decompressed files
	if [[ -n "$decompress" && -n "$noClobber" && -s "$downloadDir/$myid.log" ]]; then
		local skipDownload=true
		while read f; do
			if [ ! -e "$downloadDir/$f" ]; then
				# if we are here then one of the extracted objects does not exist
				skipDownload=false
				break
			fi
		done < $downloadDir/$myid.log
		if $skipDownload; then
			echo "Skipping download for $myid"
			return 1
		fi
	elif [[ -z "$decompress" && -n "$noClobber" && -s "$downloadDir/$filename" ]]; then
		echo "Skipping download for $myid"
		return 1
	fi
	tempDir=$(mktemp -d -p $downloadDir -t tempXXXXXX)
	if [ -z "$tempDir" ]; then
		echo "ERROR: failed to create temp directory to store download"
		return $downloaded
	fi
	pushd $tempDir > /dev/null
	# if file already exists just decompress
	if [[ -n "$decompress" && -n "$noClobber" && -s "$downloadDir/$filename" ]]; then
		downloaded=1
		decompressFile "$downloadDir/$filename"
		# store a log file to prevent script from downloading again
		find -not -name . > $downloadDir/$myid.log
		# set dot glob to move hidden files and directories
		shopt -s dotglob
		# move everything from the temp folder to intended download directory
		mv * $downloadDir
		# unset dot glob to avoid trouble
		shopt -u dotglob
	elif ! curl -OJ --header "X-Auth-Token: $token" $gdcapiurl/data/$myid; then
		echo "curl error $?"
	elif [ -f "$myid" ]; then
		cat $myid
	elif [ -f "data" ]; then
		cat data
	elif [ -n "$(ls -A)" ]; then
		# download successful
		downloaded=1
		if [ -n "$decompress" ]; then
			decompressFile "$(ls -A)"
			# store a log file to prevent script from downloading again
			find -not -name . > $downloadDir/$myid.log
		fi
		# set dot glob to move hidden files and directories
		shopt -s dotglob
		# move everything from the temp folder to intended download directory
		mv * $downloadDir
		# unset dot glob to avoid trouble
		shopt -u dotglob
	fi
	popd > /dev/null
	cleanup
	return $downloaded
}

# manifest should be either .txt or .json format
function convertManifest(){
	if basename "$manifest" | grep -q '\.json$'; then
		cmd="grep 'object_id' $manifest | grep -Eo '[0-9a-f]{8}-([0-9a-f]{4}-){3}[0-9a-f]{12}'"
	elif basename "$manifest" | grep -q '\.txt$'; then
		cmd="tail -n +2 $manifest | cut -f 1"
	else
		echo "Manifest must be either .json or .txt"
		return 1
	fi
	guidsArray=($(bash -c "$cmd"))
}

function downloadWithToken(){
	if [ -n "$manifest" ]; then
		convertManifest || return 1
	else
		guidsArray=($(convertJsonToArrayNoQuotes $guids))
	fi
	local myArray=()
	for guid in "${guidsArray[@]}"; do
		singleDownloadGET $guid && myArray+=($guid)
	done
	guidsArray=(${myArray[@]})
	return ${#myArray[@]}
}

function convertJsonToArrayNoQuotes(){
	#echos out a string that can be converted to a bash array to get around not being able to return string or array
	local string=$1
	if [ ${string:0:1} = '[' ];then
		#remove square brackets at beginning and end with substring
		#put spaces with , replacement
		#remove all " in this case
		string=$(echo "${string:1:${#string}-2}" | sed 's/\,/ /g' | sed 's/\"/ /g')
	fi
	echo "$string"
}

function downloadErrorCheck(){
	local cmd="$1"
	local log="$2"
	echo $cmd
	$cmd 2> >(tee -a $log >&2)
	local rtn=$?
	local errors=$(grep -i 'error' $log)
	if [ $rtn != 0 ]; then
		echo "Exiting with error code: $rtn"
		exitCode=$rtn
		return $rtn
	elif [ -n "$errors" ]; then
		echo "Exiting with error message:"
		echo "$errors"
		exitCode=1
		return $exitCode
	fi
}

function singleDownload(){
	[ -z $guidsArray ] && guidsArray=($(convertJsonToArrayNoQuotes $guids))
	for guid in "${guidsArray[@]}"; do
		local cmd="gen3-client download-single --profile=$profile --no-prompt --guid=$guid ${flags[@]}"
		local log="/tmp/log$guid"
		downloadErrorCheck "$cmd" "$log" || return
	done
}

function multiDownload(){
	echo "Downloading using manifest"
	local cmd="gen3-client download-multiple --profile=$profile --no-prompt ${flags[@]}"
	local log="/tmp/logManifest"
	downloadErrorCheck "$cmd" "$log"
}

#check if both guid and manifest given
#For now force user to use one or the other - otherwise this may cause difficulties with auto-multithread with the manifest being downloaded each time
#they can use two instances of the widget if they really want to do this and not merge the manifest
#if we pass a flag to indicate multi-thread execution (which we may) then we might modify this
if [[ -n $guids && -n $manifest ]]; then
	echo "Please choose either a GUID or a manifest file not both"
	echo "You can merge the GUID into the manifest or use two instances of the widget to download both"
	exit 1
fi

#First try with the old api
if [[ -n "$gdctoken" && -f "$gdctoken" ]]; then
	token=$(cat $gdctoken)
	#if there are no files left in guidsArray then all files have been successfully downloaded - exit
	downloadWithToken && exit 0
fi

echo "no gdc token given - use gen3 fence to download"
echo "Attempting to authenticate using gen3 fence service"
#Authenticate the container by copying or creating a config
#the commons might change so we do not hardcode the url
[ -z $datacommons_url ] && datacommons_url="https://nci-crdc.datacommons.io/"

#assume config file basename is config when given by user
#otherwise assume that the file given by user is a credentials file
credBasename=$(basename $cred)
if [ $credBasename == "config" ]; then
	mkdir -p /root/.gen3
	cp $cred /root/.gen3/config
elif ! gen3-client configure --profile=$profile --cred=$cred --apiendpoint=$datacommons_url; then
	echo "was not successful in creating new profile"
	exit 1
fi

#check if we have a configuration defined now
if [ -f "/root/.gen3/config" ]; then
	#now we can begin download
	flags=( "$@" )
	gen3-client auth --profile=$profile
	if [ -z $manifest ]; then
		singleDownload
	else
		multiDownload
	fi
else
	echo "must provide a valid config or credentials file"
	exit 1
fi

exit $exitCode
