FROM ubuntu:18.04
MAINTAINER lhhung<lhhung@uw.edu>

ENV DEBIAN_FRONTEND noninteractive
ENV HOME /root
ENV PIP_DISABLE_PIP_VERSION_CHECK 1

#base files/utils to be used inside container
RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        curl \
        feh \
        fluxbox \
        fonts-wqy-microhei \
        gtk2-engines-murrine \
        jq \
        language-pack-gnome-zh-hant \
        language-pack-zh-hant \
        libbz2-dev \
        libgl1-mesa-dri \
        liblzma-dev \
        libqt5webkit5-dev \
        libssl1.0 \
        mesa-utils \
        nano \
        net-tools \
        nginx \
        pwgen \
        python3-pyqt5 \
        python3-pyqt5.qtsvg \
        python3-pyqt5.qtwebkit \
        rsync \
        supervisor \
        ttf-ubuntu-font-family \
        virtualenv \
        wget \
        x11vnc \
        xterm \
        xvfb \
        zlib1g-dev \
    && rm -rf /var/lib/apt/lists/*

#files for vnc framebuffer
ADD noVNC /noVNC/
RUN dpkg -i /noVNC/x11vnc*.deb \
    && apt-get autoremove -y --purge

#files for web interface noVNC
ADD web /web/
RUN apt-get update \
    && apt-get install -y \
        build-essential \
        python-dev \
        python-pip \
    && pip install -r /web/requirements.txt \
    && apt-get remove -y --purge --auto-remove build-essential \
    && rm -rf /var/lib/apt/lists/*

#files for orange and biodepot
RUN ln -fs /bin/bash /bin/sh
RUN virtualenv --python=python3 --system-site-packages orange3venv
RUN source orange3venv/bin/activate
COPY orange3 orange3
RUN apt-get update \
    && apt-get install -y \
        build-essential \
        python3-dev \
        python3-pip \
    && python3 -m pip install --upgrade \
        pip==20.0.1 \
        setuptools \
        wheel \
    && pip3 install -r orange3/requirements-core.txt \
    && pip3 install -r orange3/requirements-gui.txt \
    && pip3 install \
        beautifulsoup4 \
        docker \
        pysam \
    && pip3 install -e orange3 \
    && apt-get remove -y --purge --auto-remove build-essential \
    && rm -rf /var/lib/apt/lists/*

#install Docker-ce
RUN apt-get update \
    && apt-get install -y \
        apt-transport-https \
        gnupg2 \
        software-properties-common \
    && curl -fsSL https://download.docker.com/linux/ubuntu/gpg | APT_KEY_DONT_WARN_ON_DANGEROUS_USAGE=1 apt-key add - \
    && add-apt-repository -y \
        "deb [arch=amd64] https://download.docker.com/linux/ubuntu bionic stable" \
    && apt-get update \
    && apt-get install -y \
        containerd.io \
        docker-ce \
        docker-ce-cli \
    && apt-get remove -y --purge --auto-remove \
        apt-transport-https \
        gnupg2 \
        software-properties-common \
    && rm -rf /var/lib/apt/lists/*

#nginx and supervisor setup
ADD supervisord.conf /etc/supervisor/conf.d/
ADD nginx.conf /etc/nginx/sites-enabled/default

#jsonpickle
RUN pip3 install --user jsonpickle

#put biodepot here and keep pip for rapid updates
ADD biodepot biodepot
RUN pip3 install -e biodepot

ADD startup.sh /
EXPOSE 6080
WORKDIR /data

#install docker-compose
RUN curl -L "https://github.com/docker/compose/releases/download/1.23.2/docker-compose-$(uname -s)-$(uname -m)" -o /usr/local/bin/docker-compose \
    && chmod +x /usr/local/bin/docker-compose

#Change app name to Bwb
RUN sed -i 's/"Orange Canvas"/"Bwb"/' /orange3/Orange/canvas/config.py

#set up some config files
COPY fluxbox_config/ /root/.fluxbox/
COPY user_config/ /root/

#patch orange3
COPY orangePatches/schemeedit.py /orange3/Orange/canvas/document/schemeedit.py
COPY orangePatches/canvasmain.py /orange3/Orange/canvas/application/canvasmain.py
COPY orangePatches/widgetsscheme.py /orange3/Orange/canvas/scheme/widgetsscheme.py
COPY orangePatches/signalmanager.py /orange3/Orange/canvas/scheme/signalmanager.py
COPY orangePatches/link.py /orange3/Orange/canvas/scheme/link.py
COPY orangePatches/signals.py /orange3/Orange/widgets/utils/signals.py
COPY orangePatches/linkitem.py /orange3/Orange/canvas/canvas/items/linkitem.py
COPY orangePatches/__main__.py /orange3/Orange/canvas/__main__.py

#add bwb start scripts
COPY scripts/startBwb.sh /usr/local/bin/startBwb.sh
COPY scripts/runDockerJob.sh /usr/local/bin/runDockerJob.sh
COPY scripts/startScheduler.sh /usr/local/bin/startScheduler.sh
COPY scripts/build_workflow_containers.sh /usr/local/bin/build_workflow_containers.sh
COPY scripts/whiteListToolDock.py /usr/local/bin/whiteListToolDock.py
COPY scripts/addWorkflowsToToolDock.py /usr/local/bin/addWorkflowsToToolDock.py
COPY scripts/addWidgetToToolDock.sh /usr/local/bin/addWidgetToToolDock.sh
COPY scripts/removeWidgetFromToolDock.sh /usr/local/bin/removeWidgetFromToolDock.sh

#add widgets and workflows
ADD widgets /widgets/
ADD workflows /workflows/
ADD notebooks /notebooks/
ADD templates /templates/
ADD coreutils /coreutils/
ADD icons /icons/
ADD tutorialFiles /tutorialFiles
ADD serverSettings.json /biodepot

#start it up
CMD /startup.sh \
    && /usr/bin/supervisord -n -c /etc/supervisor/supervisord.conf
