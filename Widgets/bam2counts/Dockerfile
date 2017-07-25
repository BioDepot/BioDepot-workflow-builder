FROM ubuntu:16.04

MAINTAINER Jiaming Hu <huj22@uw.edu>

# R pre-requisites
RUN echo 'deb http://cran.rstudio.com/bin/linux/ubuntu xenial/' | tee -a /etc/apt/sources.list

RUN gpg --keyserver keyserver.ubuntu.com --recv-key E084DAB9 && \
    gpg -a --export E084DAB9 | apt-key add -

RUN apt-get update && \
    apt-get install -y --no-install-recommends apt-utils \
    build-essential gcc aptitude \
    libcurl4-openssl-dev libssl-dev libxml2-dev \
    wget r-base r-base-dev && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*


#install components of bioconductor for networkBMA
RUN Rscript -e "source('https://bioconductor.org/biocLite.R');biocLite(c('BMA','Rcpp','RcppArmadillo','RcppEigen','BH','leaps','XML', 'xml2'),ask=FALSE)" && \
    Rscript -e "source('https://bioconductor.org/biocLite.R');biocLite(c('GenomicFeatures','GenomicAlignments','BiocParallel','DESeq2'),ask=FALSE)"

ADD bam2counts.R /home/root/
ADD counts2DESeq.R /home/root/

CMD ["Rscript /home/root/bam2counts.R"]


