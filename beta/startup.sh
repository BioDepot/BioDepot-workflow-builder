#!/bin/bash

# Determine an available socket
((display=$(ls -rv /tmp/.X11-unix/ | grep -om1 '[0-9]*') + 1))
sed -i "/^command/ s/ :[0-9]\+/ :${display}/" /etc/supervisor/conf.d/supervisord.conf

mkdir -p /var/run/sshd
cd /web && ./run.py > /var/log/web.log 2>&1 &
nginx -c /etc/nginx/nginx.conf

cp /root/.fluxbox/bwb.svg /orange3/Orange/canvas/icons/orange-canvas.svg
echo `hostname` > /etc/dockerid
#startScheduler.sh &> /tmp/schedulerLog &
