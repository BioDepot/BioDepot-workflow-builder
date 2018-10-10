#!/bin/bash
docker build  -t debian-star-build -f Dockerfile-build .
docker run --rm  -v ${PWD}:/data debian-star-build cp STAR /data/STARbin
docker build -t biodepot/debian-star:8.11-slim-2.6.0c .
