FROM ubuntu:22.04
RUN apt-get update && apt-get -y install fastqc default-jre  libfindbin-libs-perl libhtsjdk-java libngs-java \
    && apt-get autoclean -y \
    && apt-get autoremove -y \
    && rm -rf /var/lib/apt/lists/* 
 
