#!/bin/bash
jsonFile=$1
cpuCount=$2
memory=$3

cmd="curl -XPOST \"http://127.0.0.1:8080/1.0.0/schedule-job?cpu_count=${cpuCount}&memory=${memory}\" -H \"Content-Type: application/json\" --data \"@${jsonFile}\""
echo $cmd
eval $cmd
