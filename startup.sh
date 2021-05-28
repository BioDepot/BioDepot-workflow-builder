#!/bin/bash

# Determine an available socket
setdisplay() {
	i=1
	while [ -S "/tmp/.X11-unix/X$i" ]; do
		let "i+=1"
	done
	sed -i "s#\(command.*:\)[0-9]*#\1$i#" /etc/supervisor/conf.d/supervisord.conf
}

setdisplay

mkdir -p /var/run/sshd
cd /web && ./run.py > /var/log/web.log 2>&1 &
nginx -c /etc/nginx/nginx.conf

cp /root/.fluxbox/bwb.svg /orange3/Orange/canvas/icons/orange-canvas.svg
echo `hostname` > /etc/dockerid
#startScheduler.sh &> /tmp/schedulerLog &
