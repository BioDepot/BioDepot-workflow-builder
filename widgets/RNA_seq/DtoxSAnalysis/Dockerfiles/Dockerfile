FROM biodepot/ubuntu-bioc-r:16.04-3.6-3.4.3
#Mostly from Trevor's Dockerfile except for versioning and clean up
MAINTAINER Ling-Hong Hung "lhhung@uw.edu"

# Install required R Packages
RUN echo 'source("https://bioconductor.org/biocLite.R")' >> /tmp/script.R && \
    echo 'biocLite(c("locfit", "edgeR", "matrixStats", "gplots", "compare"), dependencies=TRUE)' >> /tmp/script.R && \
    Rscript /tmp/script.R && \
    rm /tmp/script.R
    
ENV HOME /home/user
WORKDIR $HOME

# Add input directories/files
ADD Programs/ $HOME/Programs/
ADD Scripts/ $HOME/Scripts/
