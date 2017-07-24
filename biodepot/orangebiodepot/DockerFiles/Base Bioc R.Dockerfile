FROM ubuntu:16.04

MAINTAINER Jiaming Hu <huj22@uw.edu>

# R pre-requisites
RUN apt-get update && \
    apt-get install -y --no-install-recommends apt-utils fonts-dejavu \
    build-essential xorg-dev gcc gcc-multilib gobjc++ gfortran libblas-dev libcairo2-dev liblzma-dev gobjc++ libreadline-dev aptitude \
    libbz2-dev libpcre3-dev libcurl4-openssl-dev libssl-dev libxml2-dev \
    software-properties-common wget texinfo texlive texlive-fonts-extra 

# Prepare R environment
ENV RHOME_DIR /usr/local/rhome

ENV PATH $RHOME_DIR/bin:$PATH

RUN mkdir -p $RHOME_DIR

#To get R's blas and lapack must compile from source NOT from deb
RUN cd /tmp && wget https://cran.r-project.org/src/base/R-latest.tar.gz && \
    tar -xzvf R-latest.tar.gz && \
    cd /tmp/R-* && ./configure && \
    cd /tmp/R-* && make -j 8 && \
    cd /tmp/R-* && make install rhome=$RHOME_DIR && rm -rf /tmp/R-*


#install components of bioconductor for networkBMA
RUN Rscript -e "source('https://bioconductor.org/biocLite.R');biocLite(c('BMA','Rcpp','RcppArmadillo','RcppEigen','BH','leaps','XML', 'xml2'),ask=FALSE)"

RUN Rscript -e "source('https://bioconductor.org/biocLite.R');biocLite(c('GenomicFeatures','GenomicAlignments','BiocParallel','DESeq2'),ask=FALSE)"

CMD ["/bin/bash"]


