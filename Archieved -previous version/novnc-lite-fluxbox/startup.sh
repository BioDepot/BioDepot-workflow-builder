#!/bin/bash

mkdir -p /var/run/sshd

# create an ubuntu user
# PASS=`pwgen -c -n -1 10`

#for alpine need to create sudo
groupadd sudo
echo '%sudo ALL=(ALL) NOPASSWD: ALL' >> /etc/sudoers
PASS=ubuntu

id -u ubuntu &>/dev/null || useradd --create-home --shell /bin/bash --user-group --groups adm,sudo ubuntu
echo "ubuntu:$PASS" | chpasswd

cd /web && ./run.py > /var/log/web.log 2>&1 &
#cd /web && ./run.py  &
nginx -c /etc/nginx/nginx.conf
exec /usr/bin/supervisord -n -c /etc/supervisor/supervisord.conf &
