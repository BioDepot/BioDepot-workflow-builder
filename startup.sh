#!/bin/bash

mkdir -p /var/run/sshd
cd /web && ./run.py > /var/log/web.log 2>&1 &
nginx -c /etc/nginx/nginx.conf

cp /root/.fluxbox/bwb.svg /orange3/Orange/canvas/icons/orange-canvas.svg
echo `hostname` > /etc/dockerid
