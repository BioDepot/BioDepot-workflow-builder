FROM biodepot/bioc-r:3.18-r-4.3.2__bookworm-slim
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
