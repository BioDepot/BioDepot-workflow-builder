#!/bin/bash
sudo docker build -t kallisto-build -f Dockerfile-build .
VERSION=`sudo docker run --rm -i -v ${PWD}:/data kallisto-build bash -c 'kallisto version' | awk '{print $3}'`
DATE=$(date +%D | sed 's/\///g')
TAG="biodepot/kallisto:${VERSION}__ubuntu-16.04__${DATE}"
sudo docker run --rm -i -v ${PWD}:/data kallisto-build bash -c 'cp /usr/local/bin/kallisto /data/.'
sudo docker build -t ${TAG} .
