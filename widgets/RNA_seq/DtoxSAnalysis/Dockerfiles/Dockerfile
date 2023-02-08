FROM biodepot/bioc-r:3.6__r-3.4.3__ubuntu-16.04__818d4296__c361c682__d890a38f
#Mostly from Trevor's Dockerfile except for versioning and clean up
MAINTAINER Ling-Hong Hung "lhhung@uw.edu"

# Install required R Packages
RUN echo 'library("BiocInstaller",lib="/usr/local/rhome/lib/")' >> /tmp/script.R && \
    echo 'biocLite(c("locfit", "edgeR", "matrixStats", "gplots", "compare"), dependencies=TRUE)' >> /tmp/script.R && \
    Rscript /tmp/script.R && \
    rm /tmp/script.R
    
ENV HOME /home/user
WORKDIR $HOME

# Add input directories/files
ADD Programs/ $HOME/Programs/
ADD Scripts/ $HOME/Scripts/
