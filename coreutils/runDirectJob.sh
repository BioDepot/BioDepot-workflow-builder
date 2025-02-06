#!/bin/bash

#Arguments are passed are the #unique_id $base_directory_for_logs commands
#pass the parent directories lockDir so that the cleanup can be saf

#echo "$@"d
#tempDir is a directory in /tmp
logPrint(){
 echo "$@" >> ${logDir}/log0
}
cleanup(){
    echo "cleaning up $1"
    rm -rf $bwbDataDir
    find $1 -type f -name 'pid.*' 2> /dev/null | while read file; do
        echo ${file}
        cid=$(cat $file)
        echo docker stop ${cid} 2> /dev/null 
        docker stop ${cid} 2> /dev/null 
    done
    cd /tmp/$tempDir && rm locks -rf
}
gatherData(){
 #array declaration is by default local
 declare -a allData=()
 #gather all the threads
 for ((i=0; i<${#myjobs[@]}; ++i)); do
    threadData=()
	keyfiles=($(ls $bwbDataDir/output$i))
    logPrint "$(ls $bwbDataDir/outputs$i) keyfiles ${#keyfiles[@]} ${keyfiles[@]}"
	if (( ${#keyfiles[@]} )); then
		allData[$i]=$(for keyfile in ${keyfiles[@]}; do
			printf '%s\0%s\0' "$keyfile" "$(cat $bwbDataDir/output$i/$keyfile)"
		done | jq -Rs ' split("\u0000") | . as $a | reduce range(0; length/2) as $i ({}; . + {($a[2*$i]): ($a[2*$i + 1]|fromjson? // .)})')
	else
	   allData[$i]={}
	fi
 done
 dataString=$(printf '%s\n' "${allData[@]}"  | jq -s .)
 logPrint "datastring $dataString"
}

exitstatus=0
outputFile=$1
shift
tempDir=$1
shift
logBaseDir=$1
shift
logDir="${logBaseDir}/${tempDir}/logs"
mkdir -p $logDir


logPrint "tempDir $tempDir"
logPrint "logBaseDir $logBaseDir"
#echo "logDir is $logDir"

myjobs=( "$@" )
logPrint "$@"
lockDir=/tmp/${tempDir}/locks
echo "mkdir -p $lockDir"
mkdir -p $lockDir
errorDir=/tmp/${tempDir}/errors
#must pass through the mountpoint - write to mountpoint/.bwb

hostDataDir=$BWBHOSTSHARE/$tempDir
bwbDataDir=$BWBSHARE/$tempDir
#wait until all the variables i.e. lockDir have been defined

trap "cleanup ${lockDir} -1 " SIGINT INT TERM

runJob(){
    for ((i=0; i<${#myjobs[@]}; ++i)); do
        #make a lock directory will fail if it exists
        #can replace this with another signaling/messaging method - but need to know when a job is taken
        if (mkdir $lockDir/lock$i 2> /dev/null ); then
            cmd="${myjobs[i]}"
            #write the pid of the process in the name of a file in the lock directory
            #this will also contain the cid of the docker process so that it can be aborted
            mkdir -p $bwbDataDir/output$i
            cmdStr="docker run -i --rm  --init --cidfile=$lockDir/lock$i/pid.$BASHPID -v $hostDataDir/output$i:/tmp/output " 
            cmdStr="$cmdStr $cmd "
            echo "$cmdStr"
            eval $cmdStr
            if [ "$?" -eq "0" ]; then
                echo "job$i exited successfully"
            else
                exitstatus=1
                mkdir -p $errorDir/job$i.$?
            fi
            rm $lockDir/lock$i/pid.$BASHPID
        fi
    done
    exit $exitstatus
}
if [ -z $NWORKERS ]; then
    NWORKERS=1
fi
for i in $(seq 1 $NWORKERS); do
    echo "starting job with thread $i"
    runJob $i 2>&1 | tee -a ${logDir}/log$i >> ${logDir}/log0 &
done	

#catch sigint and term and cleanup anyway
wait
gatherData
logPrint "output is $dataString"
echo "$dataString" > $outputFile
cleanup ${lockDir}
exit $exitstatus


