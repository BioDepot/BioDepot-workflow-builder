FROM ubuntu:16.04
MAINTAINER Ling-Hong Hung
RUN apt-get update && apt-get install -y libhdf5-10 zlib1g \
    && apt-get autoclean -y \
    && apt-get autoremove -y \
    && rm -rf /var/lib/apt/lists/*
COPY ./kallisto /bin/kallisto
COPY runKallisto.sh /usr/local/bin/

