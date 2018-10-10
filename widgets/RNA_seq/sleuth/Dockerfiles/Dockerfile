FROM biodepot/rbase:3.5.1__ubuntu-18.04

#For sleuth

#install libgit2 separate as it has a (false) dependency on libcurl4-gnutils-dev 
#which conflicts with the required libcur4-openssl-dev
#install libcurl4 here and curl as libcurl4-openssl-dev seems buggy

RUN apt-get update && apt-get install -y libgit2-dev libcurl4 curl \    
    && apt-get autoclean -y \
    && apt-get autoremove -y \
    && rm -rf /var/lib/apt/lists/*

RUN apt-get update && apt-get install -y libssl-dev libcurl4-openssl-dev libxml2-dev \
    && apt-get autoclean -y \
    && apt-get autoremove -y \
    && rm -rf /var/lib/apt/lists/*
COPY install_18.04.R /home/root/install.R

RUN Rscript /home/root/install.R
COPY runSleuth.sh /

