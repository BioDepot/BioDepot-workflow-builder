#!/bin/bash

awsdir=$1
bucket=$2
outputDir=$3

error=0
mkdir -p outputDir
cp -r $awsdir/* /root/.aws || exit 1

copy_wildcard(){
 echo "parsing wildcards in $1"
 local my_str=$1
 local max_attempts=4
 #if there is no / then we search from the base bucket
 local my_glob=""
 local wildcard=$my_str
 if [[ $glob == */* ]]; then
	#split into a glob (directory) and wildcard string
	no_wc="${my_str%%['!'@#\$%^\&*()+]*}"
	my_glob="${no_wc%/*}"/
	wildcard="${my_str#$my_glob}"
 fi
 local command=(nice aws s3 cp --exclude "*" --include="$wildcard" --recursive s3://$bucket/$my_glob $outputDir)
 local attempts
 for attempts in {1..$max_attempts}; do
	echo "${command[@]}"
 	if "${command[@]}" ; then
   	return
 	fi
 done
 echo "error in ${command[@]}"
 exit 1

}

copy_directory(){
 echo "copying directory object $1"
 local attempts
 local max_attempts=2
 local command=(nice aws s3 cp --recursive s3://$bucket/$1 $outputDir/$dest) 
 for attempts in {1..$max_attempts}; do
	echo "${command[@]}"
 	if "${command[@]}" ; then
   	return
 	fi
 done
 echo "error in ${command[@]}"
 exit 1

}

copy_file(){
 echo "copying file object $1"
 destination=basename $1
 local attempts
 local max_attempts=2
 local command=(nice aws s3 cp s3://$bucket/$1 $outputDir/$dest) 
 for attempts in {1..$max_attempts}; do
	echo "${command[@]}"
 	if "${command[@]}" ; then
   	return
 	fi
 done
 echo "error in ${command[@]}"
 exit 1
}

copy(){
   local my_glob=$1
   echo "$my_glob"
   if [[ $my_glob == *['!'@#\$%^\&*()+]* ]]; then
     copy_wildcard $my_glob || error=1
   elif [ "${my_glob: -1}" == "/" ]; then
     copy_directory $my_glob || error=1
   else
    copy_file $my_glob || error=1
   fi	
}

multiCopy(){
 lasti=$((${#globs[@]} - 1))
 for i in $(seq 0 ${lasti}); do
  if ( mkdir $lockDir/lock$i 2> /dev/null ); then
   glob=${globs[i]}
   echo "thread $1 copying $glob"
   copy $glob
  fi
 done
}
if [ -z $DIRS ] || [ "$DIRS" == "[]" ]; then
    echo "no bucket object given to download"
	exit 1
fi
globs=( $(echo $DIRS | jq -r '.[]') )

if [ -z $nThreads ] || (( $nThreads == 1 )) || (( $nThreads == 0 )); then
	#use single thread
	echo "Using single thread"
	for glob in "${globs[@]}"; do
		copy $glob
	done
else
	lockDir=/tmp/locks.$$
	mkdir -p $lockDir
	for i in $(seq 2 $nThreads); do
	  multiCopy $i &
	done
	multiCopy 1 &
	wait
	rm -rf $lockDir
fi
exit $error
