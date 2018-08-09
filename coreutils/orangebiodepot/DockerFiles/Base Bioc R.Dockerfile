FROM ubuntu:18.04

MAINTAINER Ling-Hong Hung lhhunghimself@gmail.com

# Prepare R environment
ENV RHOME_DIR /usr/local/rhome
ENV PATH $RHOME_DIR/bin:$PATH
RUN mkdir -p $RHOME_DIR

# R pre-requisites
#To get R's blas and lapack must compile from source NOT from deb

RUN apt-get update && \
    apt-get install -y --no-install-recommends apt-utils fonts-dejavu \
    build-essential xorg-dev gcc gcc-multilib gobjc++ gfortran libblas-dev libcairo2-dev liblzma-dev gobjc++ libreadline-dev aptitude \
    libbz2-dev libpcre3-dev libcurl4-openssl-dev libssl-dev libxml2-dev \
    software-properties-common wget texinfo texlive texlive-fonts-extra 

RUN cd /tmp && wget https://cran.r-project.org/src/base/R-3.5.0.tar.gz && \
    tar -xzvf R-latest.tar.gz && \
    cd /tmp/R-* && ./configure && \
    cd /tmp/R-* && make -j 8 && \
    cd /tmp/R-* && make install rhome=$RHOME_DIR && rm -rf /tmp/R-*

RUN Rscript -e "source('https://bioconductor.org/biocLite.R')
CMD ["/bin/bash"]


