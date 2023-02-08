FROM ubuntu:22.04

MAINTAINER Ling-Hong Hung lhhunghimself@gmail.com

#needed so that there is no query for time zone
ARG DEBIAN_FRONTEND=noninteractive

#define different versions
ARG RBASE_VERSION=4
ARG RBASE_SUBVERSION=4.2.2
ARG BIOC_VERSION=3.16
#older versions of bioconductor us BiocInstaller
ARG BIOC_INSTALLER=BiocManager

ENV TZ Etc/UTC
# Prepare R environment
ENV RHOME_DIR /usr/local/rhome
ENV PATH $RHOME_DIR/bin:$PATH
RUN mkdir -p $RHOME_DIR

# R pre-requisites
#To get R's blas and lapack must compile from source NOT from deb

RUN apt-get update && \
    apt-get install -y --no-install-recommends apt-utils fonts-dejavu \
    build-essential xorg-dev gcc gobjc++ gfortran libblas-dev libcairo2-dev liblzma-dev libreadline-dev aptitude \
    libbz2-dev libpcre3-dev libcurl4-openssl-dev libssl-dev libxml2-dev \
    software-properties-common wget texinfo texlive texlive-fonts-extra default-jdk && \
	cd /tmp && wget https://cran.r-project.org/src/base/R-$RBASE_VERSION/R-$RBASE_SUBVERSION.tar.gz && \
    tar -xzvf R-$RBASE_SUBVERSION.tar.gz && \
    cd /tmp/R-* && ./configure && \
    cd /tmp/R-* && make -j 8 && \
    cd /tmp/R-* && make install rhome=$RHOME_DIR \
    && rm -rf /tmp/R-* \
    && apt-get remove -y build-essential wget\
    && apt-get autoclean -y \
    && apt-get autoremove -y \
    && rm -rf /var/lib/apt/lists/*

RUN echo R_UNZIPCMD="'"internal"'"  > /root/.Renviron


#Versions 3.5 and newer use BiocManager
RUN if [ "$BIOC_INSTALLER" = BiocManager ] ; then \
    Rscript -e 'install.packages("BiocManager",repos="http://cran.us.r-project.org",lib="/usr/local/rhome/lib/")' ; \
  else \
    RUN Rscript -e 'install.packages("BiocInstaller",repos="http://bioconductor.org/packages/$BIOC_VERSION/bioc",lib="/usr/local/rhome/lib/")' ; \
  fi
