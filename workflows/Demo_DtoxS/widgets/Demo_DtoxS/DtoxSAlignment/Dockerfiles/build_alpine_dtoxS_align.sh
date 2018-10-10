#!/bin/sh
 sudo docker build -t "alpine-bwa-builder" build
 sudo docker run --rm -v ${PWD}:/local alpine-bwa-builder cp bwa /local/bwa
 sudo docker rmi alpine-bwa-builder
 sudo docker build -t "biodepot/dtoxs_alignment:1.0__alpine-3.7__bwa-0.715-r1140__python-2.7.14__072818" ./
 sudo rm bwa
