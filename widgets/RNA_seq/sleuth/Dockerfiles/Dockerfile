FROM biodepot/bioc-r:3.16-ubuntu-22.04-r-4.2.2__e3b0c442__3863bd4b__394dc8fc

ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get update && apt-get install -y libssl-dev libgit2-dev libcurl4-gnutls-dev libxml2-dev libharfbuzz-dev r-bioc-rhdf5 libfribidi-dev libtiff-dev curl\
    && apt-get autoclean -y \
    && apt-get autoremove -y \
    && rm -rf /var/lib/apt/lists/*
  
COPY install.R /home/root/install.R
RUN Rscript /home/root/install.R
COPY runSleuth.sh /

