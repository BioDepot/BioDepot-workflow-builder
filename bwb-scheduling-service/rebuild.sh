#!/bin/bash
MY_PATH="`dirname \"$0\"`"              # relative
MY_PATH="`( cd \"$MY_PATH\" && pwd )`"  # absolutized and normalized
if [ -z "$MY_PATH" ] ; then
  exit 1  
fi
cd $MY_PATH && docker build -t biodepot/bwb-scheduling-service . && docker build -t biodepot/local-agent -f Dockerfile-agent .
