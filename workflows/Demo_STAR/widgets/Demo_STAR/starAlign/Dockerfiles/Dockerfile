#note that version 2.60 does not compile with 18.04
From debian:8.11-slim
MAINTAINER "Ling-Hong Hung" lhhunghimself@gmail.com 
RUN apt-get update && apt-get install -y gzip bzip2 \
    && apt-get autoclean -y \
    && apt-get autoremove -y \
    && rm -rf /var/lib/apt/lists/*
COPY STARbin /bin/STAR
COPY runstar.sh /bin/runstar.sh

