#!/bin/bash

#find container id
hostid=$(hostname)
[ -z "$hostid" ] && hostid=$(head -1 /proc/self/cgroup | cut -d/ -f3)
docker kill $hostid