#!/bin/bash
docker build  -t debian-star-build -f Dockerfile-build .
docker run --rm  -v ${PWD}:/data debian-star-build cp STAR /data/STARbin
docker build -t biodepot/star:2.6.0c__debian-8.11-slim__072918 .
