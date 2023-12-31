FROM biodepot/jupyter-bioc-r:6.5.4__bioc-r-3.18-r-4.3.2__bookworm-slim

ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get update && apt-get install -y libssl-dev libgit2-dev libcurl4-gnutls-dev libxml2-dev libharfbuzz-dev r-bioc-rhdf5 libfribidi-dev libtiff-dev curl\
    && apt-get autoclean -y \
    && apt-get autoremove -y \
    && rm -rf /var/lib/apt/lists/*
  
COPY install.R /home/root/install.R
RUN Rscript /home/root/install.R
COPY installDeps.R /home/root/installDeps.R
RUN Rscript /home/root/installDeps.R


