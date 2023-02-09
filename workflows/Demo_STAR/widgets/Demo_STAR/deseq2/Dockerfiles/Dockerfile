FROM biodepot/bioc-r:3.16-ubuntu-22.04-r-4.2.2__e3b0c442__3863bd4b__394dc8fc
RUN apt-get update && apt-get -y install libfontconfig1-dev libharfbuzz-dev libfribidi-dev git libtiff-dev libzmq3-dev build-essential libgit2-dev \
    && echo 'library("BiocManager",lib="/usr/local/rhome/lib/")' >> /tmp/script.R \
    && echo 'BiocManager::install("DESeq2", dependencies=TRUE, update=TRUE, ask=FALSE)' >> /tmp/script.R \
    && Rscript /tmp/script.R \
    && apt-get -y remove build-essential \
    && rm /tmp/script.R \
    && apt-get autoclean -y \
    && apt-get autoremove -y \
    && rm -rf /var/lib/apt/lists/*
COPY runDESeq2.sh /bin/
