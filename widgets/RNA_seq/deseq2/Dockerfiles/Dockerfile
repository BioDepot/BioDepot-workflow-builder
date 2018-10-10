FROM biodepot/bioconductor:3.7-ubuntu-16.04-r-3.5.1
MAINTAINER Ling-Hong Hung
RUN echo 'source("https://bioconductor.org/biocLite.R")' >> /tmp/script.R && \
    echo 'biocLite("DESeq2", dependencies=TRUE)' >> /tmp/script.R && \
    Rscript /tmp/script.R && \
    rm /tmp/script.R
COPY runDESeq2.sh /bin/
