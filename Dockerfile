FROM ubuntu:18.04
ARG TARGETARCH
# Setup demo environment variables
ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        curl \
        dbus-x11 \
        feh \
        firefox \
        fluxbox \
        fonts-wqy-microhei \
        gtk2-engines-murrine \
        gvfs-backends \
        jq \
        language-pack-gnome-zh-hant \
        language-pack-zh-hant \
        libbz2-dev \
        libgl1-mesa-dri \
        liblzma-dev \
        libqt5webkit5-dev \
        libssl1.0 \
        libwebkit2gtk-4.0 \
        mesa-utils \
        nano \
        nemo \
        net-tools \
        nginx \
        openssh-server \
        novnc \
        pwgen \
        python3-pyqt5 \
        python3-pyqt5.qtsvg \
        python3-pyqt5.qtwebkit \
        rsync \
        supervisor \
        ttf-ubuntu-font-family \
        wget \
        x11vnc \
        xdg-utils \
        xterm \
        xvfb \
        zlib1g-dev \
    && rm -rf /var/lib/apt/lists/*

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
ARG DOCKER_GPG=/etc/apt/keyrings/docker.gpg
RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        gnupg \
    && mkdir -p /etc/apt/keyrings \
    && curl -fsSL https://download.docker.com/linux/ubuntu/gpg | gpg --dearmor -o $DOCKER_GPG \
    && arch=$(dpkg --print-architecture) \
    && codename=$(awk -F= '/CODENAME/ {print $2;exit}' /etc/os-release) \
    && echo "deb [arch=$arch signed-by=$DOCKER_GPG] https://download.docker.com/linux/ubuntu $codename stable" > \
        /etc/apt/sources.list.d/docker.list \
    && apt-get update \
    && apt-get install -y --no-install-recommends \
        containerd.io \
        docker-ce \
        docker-ce-cli \
    && apt-get remove -y --purge --auto-remove \
        gnupg \
    && rm -rf /var/lib/apt/lists/*

#jsonpickle
RUN pip3 install --user jsonpickle

#put biodepot here and keep pip for rapid updates
ADD widgets widgets
ADD biodepot biodepot

#This script is necessary for customization
COPY scripts/generate_setup.sh /usr/local/bin/generate_setup.sh

RUN pip3 install -e biodepot

#install docker-compose
RUN curl -L "https://github.com/docker/compose/releases/download/1.23.2/docker-compose-$(uname -s)-$(uname -m)" -o /usr/local/bin/docker-compose \
    && chmod +x /usr/local/bin/docker-compose

#Change app name to Bwb
RUN sed -i 's/"Orange Canvas"/"Bwb"/' /orange3/Orange/canvas/config.py

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
COPY orangePatches/discovery.py /orange3/Orange/canvas/registry/discovery.py

#add bwb start scripts
COPY scripts/startBwb.sh /usr/local/bin/startBwb.sh
COPY scripts/startSingleBwb.sh /usr/local/bin/startSingleBwb.sh
COPY scripts/runDockerJob.sh /usr/local/bin/runDockerJob.sh
COPY scripts/startScheduler.sh /usr/local/bin/startScheduler.sh
COPY scripts/build_workflow_containers.sh /usr/local/bin/build_workflow_containers.sh
COPY scripts/whiteListToolDock.py /usr/local/bin/whiteListToolDock.py
COPY scripts/addWorkflowsToToolDock.py /usr/local/bin/addWorkflowsToToolDock.py
COPY scripts/addWidgetToToolDock.sh /usr/local/bin/addWidgetToToolDock.sh
COPY scripts/removeWidgetFromToolDock.sh /usr/local/bin/removeWidgetFromToolDock.sh
COPY scripts/generate_setup.sh usr/local/bin/generate_setup.sh
COPY executables /usr/local/bin/executables
#add widgets and workflows

RUN groupadd ftpaccess
COPY sshd_config /etc/ssh/sshd_config
COPY startSftp.sh /usr/local/bin/startSftp.sh
COPY scripts/findResolution.sh /usr/local/bin/findResolution.sh

ADD workflows /workflows/
ADD notebooks /notebooks/
ADD templates /templates/
ADD coreutils /coreutils/
ADD icons /icons/
ADD tutorialFiles /tutorialFiles
ADD serverSettings.json /biodepot

ADD websockify /websockify
ADD noVNC /noVNC
ADD startup.sh /
ADD nginx.conf /etc/nginx/sites-enabled/default
ADD supervisord.conf /etc/supervisor/conf.d/

#setup arch dependent executables

COPY scripts/setExecutablesArch.sh /usr/local/bin/setExecutablesArch.sh
RUN setExecutablesArch.sh /usr/local/bin/executables $TARGETARCH

WORKDIR /data
CMD /startup.sh && /usr/bin/supervisord -n -c /etc/supervisor/supervisord.conf
EXPOSE 6080
EXPOSE 22

