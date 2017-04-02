FROM ubuntu:latest
MAINTAINER Daniel Kristiyanto

ENV DEBIAN_FRONTEND noninteractive
EXPOSE 6080
WORKDIR /root

## REQUIRED PACKAGES
RUN apt-get update -y
RUN apt-get install -y git x11vnc wget unzip xvfb openbox geany menu \
    build-essential python3 python3-dev python3-pip virtualenv libssl-dev \
    net-tools rox-filer feh python3-pyqt5 libqt5webkit5-dev python3-pyqt5.qtsvg \
    python3-pyqt5.qtwebkit

## NOVNC
RUN git clone https://github.com/kanaka/noVNC.git && \
    cd noVNC/utils && git clone https://github.com/kanaka/websockify websockify

## ORANGE3
RUN rm /bin/sh && ln -s /bin/bash /bin/sh
RUN virtualenv --python=python3 --system-site-packages orange3venv
RUN source orange3venv/bin/activate
RUN git clone https://github.com/biolab/orange3.git 
RUN pip3 install --upgrade pip
RUN pip3 install -r orange3/requirements-core.txt
RUN pip3 install -r orange3/requirements-gui.txt
RUN pip3 install docker numpy pysam
RUN pip3 install -e orange3

## BIODEPOT
ADD biodepot biodepot
RUN pip3 install -e biodepot

## CLEAN UP
RUN apt-get autoclean && apt-get autoremove && rm -rf /var/lib/apt/lists/*

## DESKTOP SETTINGS
ADD Desktop/menu.xml  /root/.config/openbox/menu.xml 
ADD Desktop/bg.png /root/.config/openbox/bg.png
RUN echo "feh /root/.config/openbox/bg.png & rox-filer /data & orange-canvas" \ 
    >> /root/.config/openbox/autostart
ADD Desktop/rc.xml /root/.config/openbox/rc.xml

## CMD
ADD Desktop/novnc.sh /root/novnc.sh
RUN chmod 0755 /root/novnc.sh



## START
WORKDIR /data
CMD /root/novnc.sh

