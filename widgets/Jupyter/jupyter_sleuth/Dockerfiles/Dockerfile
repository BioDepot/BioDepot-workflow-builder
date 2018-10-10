FROM biodepot/rbase:3.4.4__ubuntu-16.04

#For sleuth

#install libgit2 separate as it has a (false) dependency on libcurl4-gnutils-dev 
#which conflicts with the required libcur4-openssl-dev

RUN apt-get update && apt-get install -y libgit2-dev \    
    && apt-get autoclean -y \
    && apt-get autoremove -y \
    && rm -rf /var/lib/apt/lists/*

RUN apt-get update && apt-get install -y libssl-dev libcurl4-openssl-dev libxml2-dev \
    && apt-get autoclean -y \
    && apt-get autoremove -y \
    && rm -rf /var/lib/apt/lists/*
COPY install.R /home/root/install.R
RUN Rscript /home/root/install.R

RUN apt-get update && apt-get -y install build-essential python3-all python3-pip \
    libncurses5-dev libncursesw5-dev libzmq3-dev zlib1g-dev \
    && pip3 install --upgrade pip \
    && pip install jupyter \
    && apt-get autoclean -y \
    && apt-get autoremove -y \
    && rm -rf /var/lib/apt/lists/*

COPY installKernel.R  /home/root/installKernel.R
RUN Rscript /home/root/installKernel.R

RUN apt-get update && apt-get -y install firefox \
    && apt-get autoclean -y \
    && apt-get autoremove -y \
    && rm -rf /var/lib/apt/lists/*
COPY .mozilla /root/.mozilla

