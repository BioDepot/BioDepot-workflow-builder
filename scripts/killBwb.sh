#!/bin/bash

#find hostname
hostid=$(head -1 /proc/self/cgroup | cut -d/ -f3)
docker kill $hostid