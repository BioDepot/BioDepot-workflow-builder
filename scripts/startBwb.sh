#/bin/bash
if [ $# -eq 0 ]; then
    cmd=orange-canvas
else
    cmd="$*"
fi

while true
do
$cmd &
pid="$!"
mkdir -p "/tmp/pid.$pid"
wait $pid
if [ -d "/tmp/pid.$pid" ]; then
    if [ -f "/tmp/pid.$pid/cmd" ]; then
        cmd= `cat /tmp/pid.$pid/cmd`
        rm -rf "/tmp/pid.$pid"
    else
        rm -rf "/tmp/pid.$pid"
        break
    fi
fi
done

