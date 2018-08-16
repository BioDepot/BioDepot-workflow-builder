#/bin/bash

while true
do
orange-canvas &
pid="$!"
echo $pid
touch /tmp/pid.$pid
wait $pid
if [ -f "/tmp/pid.$pid" ]; then
	rm "/tmp/pid.$pid"
	break
fi
done

