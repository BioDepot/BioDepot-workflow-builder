#!/bin/bash
namespace=$1
jobDone=""
monitor(){
	retvalue=`curl -s localhost:8080/1.0.0/status?namespace=${namespace}`
	echo $retvalue | grep 'does not exist' &> /dev/null
	if [ $? == 0 ]; then
		exit
	fi
		
	echo $retvalue | jq  'values[]' | grep 1 &> /dev/null
	if [ $? == 0 ]; then
		sleep 1
	else
		jobDone=1
	fi
	
}

while [ -z "$jobDone" ]; do
	monitor
done
