#!/bin/bash

#First argument is the number of threads - the other arguments are passed are the commands
#A lock directory will be made by the starting process with the pid in the name
#It spawns nthreads runJob processes and then 

#if interrupted - still cleanup
myjobs=( "$@" )

nThreads=$1
lockDir=/tmp/locks.$$
mkdir -p $lockDir
echo nThreads $nThreads
#wait until all the variables i.e. lockDir have been defined
trap "cleanup ${lockDir} " SIGINT INT TERM
cleanup(){
    echo "cleaning up $1"
    echo " find $1 -type f -name 'pid.* "
    find $1 -type f -name 'pid.*' 2> /dev/null | while read file; do
        echo ${file}
        cid=$(cat $file)
        echo docker stop ${cid} 2> /dev/null 
        docker stop ${cid} 2> /dev/null 
    done
    rm ${lockDir} -rf
    exit
}
runJob(){
    for ((i=1; i<${#myjobs[@]}; ++i)); do
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

for i in $(seq 1 $nThreads); do
      runJob $i &
done
#catch sigint and term and cleanup anyway

wait
cleanup ${lockDir}
exit

