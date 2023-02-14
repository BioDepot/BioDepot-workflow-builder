FROM debian:bookworm-slim

#needed so that there is no query for time zone
ARG DEBIAN_FRONTEND=noninteractive

<<<<<<< HEAD
#define different versions
ARG RBASE_VERSION=4
ARG RBASE_SUBVERSION=4.2.2

#older versions of bioconductor us BiocInstaller
ARG BIOC_INSTALLER=BiocManager

ENV TZ Etc/UTC
# Prepare R environment
ENV RHOME_DIR /usr/local/rhome
ENV PATH $RHOME_DIR/bin:$PATH
RUN mkdir -p $RHOME_DIR

# R pre-requisites
#To get R's blas and lapack must compile from source NOT from deb
RUN apt-get -y update
RUN apt-get install -y epiphany epiphany-browser epiphany-browser-data xdg-desktop-portal-gnome xdg-desktop-portal-gtk 
RUN apt-get remove -y dbus-user-session
RUN apt-get install dbus-x11

RUN apt-get update && \
    apt-get install -y --no-install-recommends apt-utils fonts-dejavu \
    build-essential xorg-dev gcc gobjc++ gfortran libblas-dev libcairo2-dev liblzma-dev libreadline-dev aptitude \
    libbz2-dev libpcre3-dev libcurl4-openssl-dev libssl-dev libxml2-dev \
    software-properties-common wget texinfo texlive texlive-fonts-extra default-jdk && \
	cd /tmp && wget https://cran.r-project.org/src/base/R-$RBASE_VERSION/R-$RBASE_SUBVERSION.tar.gz && \
    tar -xzvf R-$RBASE_SUBVERSION.tar.gz && \
    cd /tmp/R-* && ./configure && \
    cd /tmp/R-* && make -j 8 && \
    cd /tmp/R-* && make install rhome=$RHOME_DIR \
    && rm -rf /tmp/R-* \
    && apt-get remove -y build-essential wget\
    && apt-get autoclean -y \
    && apt-get autoremove -y \
    && rm -rf /var/lib/apt/lists/*

RUN echo R_UNZIPCMD="'"internal"'"  > /root/.Renviron

RUN apt-get update && apt-get -y install firefox-esr xdg-utils python3-all sudo python3-pip libfontconfig1-dev libharfbuzz-dev libfribidi-dev git libtiff-dev libzmq3-dev build-essential libgit2-dev \
 && python3 -m pip install pip==22.3.1 \
 && pip install jupyter \
 && R -e "install.packages(c('IRdisplay', 'devtools' ),repos = 'http://cran.us.r-project.org'); devtools::install_github('IRkernel/IRkernel'); IRkernel::installspec()" \
 && apt-get remove -y build-essential \
 && apt-get autoclean -y \
 && apt-get autoremove -y \
 && rm -rf /var/lib/apt/lists/*
#This gives an error signal even though it is successful we put it in a combo so that the error is ignored
RUN xdg-settings set default-web-browser firefox.desktop || echo 1
COPY .mozilla ~/.mozilla

#To prevent helvetica from being used 
RUN echo 'options(bitmapType="cairo")' >> ~/.Rprofile
=======
sudo apt update && sudo apt install midori
>>>>>>> 5b6f3b2871e51f9991a7496d22445d0bae789a86
