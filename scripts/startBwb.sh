#/bin/bash
if [ $# -eq 0 ]; then
    cmd=orange-canvas
else
    cmd="$*"
fi
while true
do
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

