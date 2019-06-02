#!/bin/bash
container=$1
logDir=$2
worker=$3

function cleanup {
	docker stop $container 2> /dev/null 
}
trap cleanup SIGINT INT TERM

echo $worker $container 2>&1 | tee -a ${logDir}/log${worker} >> ${logDir}/log0
docker start -a -i $container 2>&1 | tee -a ${logDir}/log${worker} >> ${logDir}/log0
