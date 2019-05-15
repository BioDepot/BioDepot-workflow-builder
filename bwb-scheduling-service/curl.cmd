#!/bin/bash

namespace=$1
jsonFile=$2
cpuCount=$3
memory=$4

cmd="curl -XPOST \"http://127.0.0.1:8080/1.0.0/schedule-job?cpu_count=${cpuCount}&memory=${memory}&namespace=${namespace}\" -H \"Content-Type: application/json\" --data \"@${jsonFile}\""
echo $cmd
eval $cmd
