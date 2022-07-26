#!/bin/bash
set -Eeo pipefail
user=$1
pass=$2
useradd  "$user" || echo "User $user already exists"
$(echo "$user:$pass" | chpasswd)  || { echo 'password change failed' ; exit 1; }
usermod -a -G ftpaccess $user
service ssh restart