#!/bin/bash

#Arguments are passed are the commands
#A lock directory will be made by the starting process with the pid in the name
#if interrupted - still cleanup

#echo "$@"
logDir=$1
shift
#echo "logDir is $logDir"

myjobs=( "$@" )
echo "$@"
lockDir=/tmp/locks.$$
echo "mkdir -p $lockDir"
mkdir -p $lockDir
echo nThreads $NWORKERS

#wait until all the variables i.e. lockDir have been defined
for ((i=1; i<${#myjobs[@]}; ++i)); do
    cmd="${myjobs[i]}"
    echo "job $i is docker run -i --rm  --init --cidfile=$lockDir/lock$i/pid.$BASHPID $cmd"
done
trap "cleanup ${lockDir} " SIGINT INT TERM

cleanup(){
    echo "cleaning up $1"
    find $1 -type f -name 'pid.*' 2> /dev/null | while read file; do
        echo ${file}
        cid=$(cat $file)
        echo docker stop ${cid} 2> /dev/null 
        docker stop ${cid} 2> /dev/null 
    done
    cd /tmp && rm locks.$$ -rf
    #cd /data/.bwb && rm $logDir
    exit
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
        rm $lockDir/lock$i/pid.$BASHPID
        fi
    done
    exit
}
if [ -z "$logDir" ]; then
	for i in $(seq 1 $NWORKERS); do
		echo "starting job with thread $i"
		runJob $i &
	done
else
	for i in $(seq 1 $NWORKERS); do
		echo "starting job with thread $i"
		runJob $i 2>&1 | tee -a /data/.bwb/${logDir}/log$i >> /data/.bwb/${logDir}/log0 &
	done	
fi
#catch sigint and term and cleanup anyway

wait
cleanup ${lockDir}
exit

