#!/bin/bash
# Arguments: output JSON file, unique proc.ID, log directory, Docker commands.
# A failed worker must not be hidden by tee or by a background subshell.
set -o pipefail

logPrint() {
    echo "$@" >> "$logDir/log0"
}

cleanup() {
    local file cid
    # Stop only containers belonging to this invocation, before removing data.
    for file in "$lockDir"/lock*/pid.*; do
        [ -f "$file" ] || continue
        cid=$(cat "$file")
        [ -z "$cid" ] || docker stop "$cid" >/dev/null 2>&1
    done
    if [ "$ownsDataDir" = true ]; then
        rm -rf -- "$bwbDataDir"
    fi
    rm -rf -- "$lockDir"
}

fail() {
    echo "Bwb job failed: $*" >&2
    if [ -d "${logDir:-}" ]; then
        logPrint "Bwb job failed: $*"
    fi
    mkdir -p "$errorDir/setup" 2>/dev/null
    exit 1
}

gatherData() {
    local i keyfile
    local -a allData=() keyfiles=()
    for ((i=0; i<${#myjobs[@]}; ++i)); do
        keyfiles=("$bwbDataDir/output$i/"*)
        if [ ${#keyfiles[@]} -gt 0 ]; then
            allData[$i]=$(for keyfile in "${keyfiles[@]}"; do
                printf '%s\0%s\0' "${keyfile##*/}" "$(cat "$keyfile")"
            done | jq -Rs 'split("\u0000") | . as $a | reduce range(0; (length/2|floor)) as $i ({}; . + {($a[2*$i]): ($a[2*$i + 1]|fromjson? // .)})') || return 1
        else
            allData[$i]='{}'
        fi
    done
    dataString=$(printf '%s\n' "${allData[@]}" | jq -s .) || return 1
    logPrint "output is $dataString"
    printf '%s\n' "$dataString" > "$outputFile"
}

if [ "$#" -lt 4 ]; then
    echo 'Usage: runDockerJob.sh OUTPUT proc.ID LOGDIR DOCKER_COMMAND ...' >&2
    exit 1
fi
outputFile=$1
tempDir=$2
logBaseDir=$3
shift 3
myjobs=("$@")
# These paths are later removed. Reject traversal and empty identifiers.
if [[ ! $tempDir =~ ^proc\.[a-zA-Z0-9_-]+$ ]]; then
    echo 'Invalid Bwb process directory identifier' >&2
    exit 1
fi
lockDir="/tmp/$tempDir/locks"
errorDir="/tmp/$tempDir/errors"
ownsDataDir=false
shopt -s nullglob

if ! mkdir -p "$logBaseDir/$tempDir/logs"; then
    logBaseDir=/tmp/.bwb
fi
logDir="$logBaseDir/$tempDir/logs"
mkdir -p "$logDir" || fail "Cannot create log directory $logDir"
logPrint "tempDir $tempDir"
logPrint "logBaseDir $logBaseDir"

if [[ ${BWBSHARE:-} != /* || ${BWBHOSTSHARE:-} != /* ]]; then
    fail 'BWBSHARE and BWBHOSTSHARE must be absolute, nonempty paths'
fi
if [[ ! ${BWBSHARE//\//} || ! ${BWBHOSTSHARE//\//} ]]; then
    fail 'The filesystem root cannot be used as the shared workspace'
fi
bwbDataDir="${BWBSHARE%/}/$tempDir"
hostDataDir="${BWBHOSTSHARE%/}/$tempDir"
NWORKERS=${NWORKERS:-1}
[[ $NWORKERS =~ ^[1-9][0-9]*$ ]] || fail 'NWORKERS must be a positive integer'

# Create rather than merely testing the parent's write bits. Existing writable
# workspaces under an unwritable mount root are deliberately supported.
mkdir -p "$BWBSHARE" || fail "Cannot create shared workspace $BWBSHARE"
mkdir "$bwbDataDir" || fail "Cannot create new job directory $bwbDataDir"
ownsDataDir=true
trap cleanup EXIT
trap 'exit 130' INT
trap 'exit 143' TERM
mkdir -p "$lockDir" || fail "Cannot create lock directory $lockDir"

# No Docker job is launched unless every metadata output directory is ready.
for ((i=0; i<${#myjobs[@]}; ++i)); do
    mkdir "$bwbDataDir/output$i" || fail "Cannot create $bwbDataDir/output$i"
done

runJob() {
    local i cmdStr rc workerStatus=0
    for ((i=$1-1; i<${#myjobs[@]}; ++i)); do
        if mkdir "$lockDir/lock$i" 2>/dev/null; then
            # Quote our generated arguments even though widget commands retain
            # their existing shell-command interface.
            printf -v cmdStr '%q ' docker run -i --rm --init \
                "--cidfile=$lockDir/lock$i/pid.$BASHPID" \
                -v "$hostDataDir/output$i:/tmp/output"
            cmdStr+="${myjobs[i]}"
            echo "$cmdStr"
            eval "$cmdStr"
            rc=$?
            if [ "$rc" -eq 0 ]; then
                echo "job$i exited successfully"
            else
                echo "job$i failed with exit code $rc" >&2
                workerStatus=1
                mkdir -p "$errorDir/job$i.$rc"
            fi
            rm -f -- "$lockDir/lock$i/pid.$BASHPID"
        fi
    done
    return "$workerStatus"
}

workerPids=()
for ((i=1; i<=NWORKERS; ++i)); do
    runJob "$i" 2>&1 | tee -a "$logDir/log$i" >> "$logDir/log0" &
    workerPids+=("$!")
done
exitstatus=0
for workerPid in "${workerPids[@]}"; do
    wait "$workerPid" || exitstatus=1
done
if [ "$exitstatus" -ne 0 ]; then
    fail 'One or more Docker workers failed; see the job log'
fi
gatherData || fail 'Cannot collect or write Docker job outputs'
