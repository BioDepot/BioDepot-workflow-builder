FROM biodepot/bioc-r:3.16-ubuntu-22.04-r-4.2.2__e3b0c442__3863bd4b__394dc8fc
RUN apt-get update && apt-get -y install firefox python3-all sudo python3-pip libfontconfig1-dev libharfbuzz-dev libfribidi-dev git libtiff-dev libzmq3-dev build-essential libgit2-dev
RUN echo 'library("BiocManager",lib="/usr/local/rhome/lib/")' >> /tmp/script.R && \
    echo 'BiocManager::install("DESeq2", dependencies=TRUE, update=TRUE, ask=FALSE)' >> /tmp/script.R &&  \
    Rscript /tmp/script.R && \
    rm /tmp/script.R
COPY runDESeq2.sh /bin/
