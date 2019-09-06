#!/bin/bash
container=$1
logDir=$2
worker=$3

pwd 2>&1 | tee -a ${logDir}/log${worker} >> ${logDir}/log0

docker start -a -i $container 2>&1 | tee -a ${logDir}/log${worker} >> ${logDir}/log0
