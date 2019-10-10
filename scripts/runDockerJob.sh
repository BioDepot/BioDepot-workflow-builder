#!/bin/bash

#Arguments are passed are the #unique_id $base_directory_for_logs commands
#pass the parent directories lockDir so that the cleanup can be saf

#echo "$@"d
#tempDir is a directory in /tmp
tempDir=$1
shift
logBaseDir=$1
shift
echo "tempDir $tempDir"
echo "logBaseDir $logBaseDir"
#echo "logDir is $logDir"

myjobs=( "$@" )
echo "$@"
lockDir=/tmp/${tempDir}/locks
echo "mkdir -p $lockDir"
mkdir -p $lockDir

errorDir=/tmp/${tempDir}/errors
echo "running job with $NWORKERS threads"

#wait until all the variables i.e. lockDir have been defined
for ((i=1; i<${#myjobs[@]}; ++i)); do
    cmd="${myjobs[i]}"
    echo "job $i is docker run -i --rm  --init --cidfile=$lockDir/lock$i/pid.$BASHPID $cmd"
done
trap "cleanup ${lockDir} -1 " SIGINT INT TERM

cleanup(){
    echo "cleaning up $1"
    find $1 -type f -name 'pid.*' 2> /dev/null | while read file; do
        echo ${file}
        cid=$(cat $file)
        echo docker stop ${cid} 2> /dev/null 
        docker stop ${cid} 2> /dev/null 
    done
    cd /tmp/$tempDir && rm locks -rf
}

runJob(){
    for ((i=0; i<${#myjobs[@]}; ++i)); do
        #make a lock directory will fail if it exists
        #can replace this with another signaling/messaging method - but need to know when a job is taken
        if (mkdir $lockDir/lock$i 2> /dev/null ); then
            cmd="${myjobs[i]}"
            #write the pid of the process in the name of a file in the lock directory
            #this will also contain the cid of the docker process so that it can be aborted
            cmdStr="docker run -i --rm  --init --cidfile=$lockDir/lock$i/pid.$BASHPID" 
            cmdStr="$cmdStr $cmd"
            echo "$cmdStr"
            eval $cmdStr
            if [ "$?" -eq "0" ]; then
                echo "job$i exited successfully"
            else
                mkdir -p $errorDir/job$i.$?
            fi
            rm $lockDir/lock$i/pid.$BASHPID
        fi
    done
    exit 0
}
logDir="${logBaseDir}/${tempDir}/logs"
mkdir -p $logDir
for i in $(seq 1 $NWORKERS); do
    echo "starting job with thread $i"
    runJob $i 2>&1 | tee -a ${logDir}/log$i >> ${logDir}/log0 &
done	

#catch sigint and term and cleanup anyway
wait
cleanup ${lockDir}
exit 0


