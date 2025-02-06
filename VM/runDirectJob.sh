#!/bin/bash

#Arguments are passed are the #unique_id $base_directory_for_logs commands
#pass the parent directories lockDir so that the cleanup can be saf

#extract the key pair values from $1
logPrint(){
 echo "$@" >> "${logDir}"/log0
}
cleanup(){
    echo "cleaning up $1"
    rm -rf "$bwbDataDir"
    find $1 -type f -name 'pid.*' 2> /dev/null | while read -r file; do
        echo "${file}"
        cid=$(cat "$file")
        echo docker stop "${cid}" 2> /dev/null 
        docker stop "${cid}" 2> /dev/null 
    done
    cd /tmp/"$tempDir" && rm locks -rf
}
gatherData(){
    #array declaration is by default local
    declare -a allData=()
    #gather all the threads
    for ((i=0; i<${#myjobs[@]}; ++i)); do
        readarray -t keyfiles < <(ls "${bwbDataDir}/output${i}")
        logPrint "$(ls "${bwbDataDir}/output${i}")" "keyfiles" "${#keyfiles[@]}" "${keyfiles[@]}"
        if (( ${#keyfiles[@]} )); then
            allData[$i]=$(for keyfile in "${keyfiles[@]}"; do
                printf '%s\0%s\0' "$keyfile" "$(cat "${bwbDataDir}/output${i}/${keyfile}")"
            done | jq -Rs 'split("\u0000") | . as $a | reduce range(0; length/2) as $i ({}; . + {($a[2*$i]): ($a[2*$i + 1] | fromjson? // .)})')
        else
            allData[$i]={}
        fi
    done
    dataString=$(printf '%s\n' "${allData[@]}"  | jq -s .)
    echo "dataString $dataString"
    logPrint "datastring $dataString"
}
exitstatus=0
runJob(){
    local i=$1
    
    for ((i=0; i<${#myjobs[@]}; ++i)); do
        #make a lock directory will fail if it exists
        #can replace this with another signaling/messaging method - but need to know when a job is taken
        if (mkdir "$lockDir/lock$i" 2> /dev/null ); then
            cmd="${myjobs[i]}"
            #write the pid of the process in the name of a file in the lock directory
            #this will also contain the cid of the docker process so that it can be aborted
            mkdir -p "$bwbDataDir/output$i"
            cmdStr="$cmdStr $cmd "
            echo "$cmdStr"
            eval "$cmdStr"
            if [ "$?" -eq "0" ]; then
                echo "job$i exited successfully"
            else
                exitstatus=1
                mkdir -p "$errorDir/job$i.$?"
            fi
            rm "$lockDir/lock$i/pid.$BASHPID"
        fi
    done
    exit $exitstatus
}



if [ -z "$@" ]; then
  echo "Usage: $0 '<json_string>'"
  exit 1
fi

json_input="$@"
echo "json_input $json_input"
#convert json_input into a json object
json_object=$(echo "$json_input" | jq -c .)
echo "json_object $json_object"
#extract the console_json object for json_object
console_json=$(echo "$json_object" | jq -r '.console_json')
logBaseDir=$(echo "$console_json" | jq -r '.baseLogDir')
tempDir=$(echo "$console_json" | jq -r '.processDir')
outputFile=$(echo "$console_json" | jq -r '.outputfile')
NWORKERS=$(echo "$console_json" | jq -r '.NWORKERS')
logDir="${logBaseDir}/${tempDir}/logs"
mkdir -p $logDir
logPrint "tempDir $tempDir"
logPrint "logBaseDir $logBaseDir"

lockDir=/tmp/${tempDir}/locks
echo "mkdir -p $lockDir"
mkdir -p "$lockDir"
errorDir=/tmp/${tempDir}/errors
#must pass through the mountpoint - write to mountpoint/.bwb

bwbDataDir=$BWBSHARE/$tempDir
#wait until all the variables i.e. lockDir have been defined

trap 'cleanup "${lockDir}" -1' SIGINT INT TERM

#extract the direct_json object for json_object
direct_json=$(echo "$json_object" | jq -r '.direct_json')
echo "direct_json $direct_json"
#find the number of members in the direct_json object
num_members=$(echo "$direct_json" | jq -r 'length')
echo "num_members $num_members"
#loop over the members of the direct_json object
myjobs=()
for ((i=0; i<num_members; i++)); do
    echo "member $i"
    #extract the member
    member=$(echo "$direct_json" | jq -r ".[$i]")
    echo "member $member"
    #this is a list of 1 element - extract the element
    member=$(echo "$member" | jq -r '.[0]')
    echo "member $member"
    #extract the command from the member
    cmd=$(echo "$member" | jq -r '.command')
    echo "cmd $cmd"
    #extract the env from the member
    env=$(echo "$member" | jq -r '.env')
    echo "env $env"
    env_assignments=$(echo "$env" | jq -r 'to_entries | map("\(.key)=\(.value|@sh)") | join(" ")')
    command_string="${env_assignments} ${cmd}"
    echo "command_with_env: ${command_string}"
    myjobs+=("$command_string")
done

echo "myjobs $myjobs"
if [ -z "$NWORKERS" ]; then
    NWORKERS=1
fi
for i in $(seq 1 "$NWORKERS"); do
    echo "starting job with thread $i"
    runJob "$i" 2>&1 | tee -a "${logDir}/log$i" >> "${logDir}/log0" &
done	

#catch sigint and term and cleanup anyway
wait
gatherData
logPrint "output is $dataString"
echo "$dataString" > "$outputFile"
cleanup "${lockDir}"
exit "$exitstatus"


