FROM ubuntu:18.04
MAINTAINER lhhung<lhhung@uw.edu>
#Dockerfile for widget development container
#comment to force rebuild
ENV DEBIAN_FRONTEND noninteractive
ENV HOME /root
ENV PIP_DISABLE_PIP_VERSION_CHECK 1
#base files/utils to be used inside container
RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        feh \
        fluxbox \
        fonts-wqy-microhei \
        gtk2-engines-murrine \
        language-pack-gnome-zh-hant \
        language-pack-zh-hant \
        libgl1-mesa-dri \
        libqt5webkit5-dev \
        libssl-dev \
        mesa-utils \
        nano \
        net-tools \
        nginx \
        pwgen \
        python3-pyqt5 \
        python3-pyqt5.qtsvg \
        python3-pyqt5.qtwebkit \
        sudo \
        supervisor \
        ttf-ubuntu-font-family \
        virtualenv \
        x11vnc \
        xterm \
        xvfb \
    && apt-get clean \
    && apt-get autoremove -y --purge \
    && rm -rf /var/lib/apt/lists/*

#files for vnc framebuffer
ADD noVNC /noVNC/
RUN apt-get update \
    && apt-get install -y \
        libssl1.0 \
        wget \
    && dpkg -i /noVNC/x11vnc*.deb \
    && apt-get clean \
    && apt-get autoremove -y --purge \
    && rm -rf /var/lib/apt/lists/*

#files for web interface noVNC
ADD web /web/
RUN apt-get update \
    && apt-get install -y \
        build-essential \
        python-dev \
        python-pip \
        python3-pip \
    && pip install --upgrade wsgiref \
    && python3 -m pip install --upgrade \
        pip \
        setuptools \
        wheel \
    && pip install -r /web/requirements.txt \
    && pip3 install docker \
    && apt-get remove -y --purge \
        build-essential \
        python-dev \
        python-pip \
        python3-pip \
    && apt-get clean \
    && apt-get autoremove -y --purge \
    && rm -rf /var/lib/apt/lists/*

#files for orange and biodepot
RUN ln -fs /bin/bash /bin/sh
RUN virtualenv --python=python3 --system-site-packages orange3venv
RUN source orange3venv/bin/activate
COPY orange3 orange3
RUN apt-get update \
    && apt-get install -y \
        build-essential \
        libbz2-dev \
        liblzma-dev \
        python-dev \
        python-pip \
        python3-dev \
        python3-pip \
        zlib1g-dev \
    && python3 -m pip install --upgrade \
        pip \
        setuptools \
        wheel \
    && pip3 install -r orange3/requirements-core.txt \
    && pip3 install -r orange3/requirements-gui.txt \
    && pip3 install \
        beautifulsoup4 \
        docker \
        pysam \
    && pip3 install -e orange3 \
    && apt-get remove -y --purge build-essential \
    && apt-get clean \
    && apt-get autoremove -y --purge \
    && rm -rf /var/lib/apt/lists/*

#install Docker-ce
RUN apt-get update \
    && apt-get install -y \
        apt-transport-https \
        curl \
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
    && apt-get remove -y --purge \
        apt-transport-https \
        curl \
        gnupg2 \
        software-properties-common \
    && apt-get clean \
    && apt-get autoremove -y --purge \
    && rm -rf /var/lib/apt/lists/*

#nginx and supervisor setup
ADD supervisord.conf /etc/supervisor/conf.d/
ADD nginx.conf /etc/nginx/sites-enabled/default

#jsonpickle
RUN pip3 install --user jsonpickle

#put biodepot here and keep pip for rapid updates
ADD biodepot biodepot
RUN pip3 install pip==20.0.1
RUN pip3 install -e biodepot

ADD startup.sh /
EXPOSE 6080
WORKDIR /data

#install rsync curl docker-compose and jq
RUN apt-get update \
    && apt-get install -y \
        curl \
        jq \
        rsync \
    && curl -L "https://github.com/docker/compose/releases/download/1.23.2/docker-compose-$(uname -s)-$(uname -m)" -o /usr/local/bin/docker-compose \
    && chmod +x /usr/local/bin/docker-compose \
    && apt-get clean \
    && apt-get autoremove -y --purge \
    && rm -rf /var/lib/apt/lists/*

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
