#/bin/bash
if [ $# -eq 0 ]; then
    cmd=orange-canvas
else
    cmd="$*"
fi
if [ -z $STARTING_WORKFLOW ]; then
workflow=""
else
    workflow="__init $STARTING_WORKFLOW"
fi
echo $workflow
echo $cmd
while true
do
pgrep -x "orange-canvas" 1>/dev/null && exit 0
echo $cmd $workflow
$cmd $workflow &
workflow="" 
pid="$!"
mkdir -p "/tmp/pid.$pid"
wait $pid
if [ -d "/tmp/pid.$pid" ]; then
    if [ -f "/tmp/pid.$pid/workflow" ]; then
        workflow=`cat /tmp/pid.$pid/workflow`
        rm -rf "/tmp/pid.$pid"
    else
        rm -rf "/tmp/pid.$pid"
        break
    fi
fi
done

