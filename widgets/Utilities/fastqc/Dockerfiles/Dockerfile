FROM ubuntu:18.04
RUN apt-get update && apt-get -y install fastqc default-jre  libfindbin-libs-perl \
    && apt-get autoclean -y \
    && apt-get autoremove -y \
    && rm -rf /var/lib/apt/lists/* 
 
