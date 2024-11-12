#!/bin/bash
DOCKER_GPG=/etc/apt/keyrings/docker.gpg
 apt-get install -y --no-install-recommends \
        curl \
        dbus-x11 \
        feh \
        firefox \
        fluxbox \
        fonts-wqy-microhei \
        geany \
        gtk2-engines-murrine \
        gvfs-backends \
        jq \
        language-pack-gnome-zh-hant \
        language-pack-zh-hant \
        libbz2-dev \
        libgl1-mesa-dri \
        liblzma-dev \
        libssl1.0 \
        libwebkit2gtk-4.0 \
        mesa-utils \
        nano \
        nemo \
        net-tools \
        openssh-server \
        pwgen \
        python3-pyqt5 \
        python3-pyqt5.qtsvg \
        python3-pyqt5.qtwebkit \
        rsync \
        ttf-ubuntu-font-family \
        wget \
        xdg-utils \
        lxterminal \
        xvfb \
        zlib1g-dev \
        zenity
  git clone https://github.com/BioDepot/BioDepot-workflow-builder.git
  cp -r Biodepot-workflow-builder/orange3 /orange3
  cp -r Biodepot-workflow-builder/orangePatches /orangePatches
apt-get install -y \
        build-essential \
        python3-dev \
        python3-pip \
    && python3 -m pip install --upgrade \
        pip==20.0.1 \
        setuptools \
        wheel 
 cd / && pip3 install -r orange3/requirements-gui.txt 
         beautifulsoup4 \
         docker \
         pysam 
 pip3 install -e orange3

apt-get install -y --no-install-recommends gnupg 
mkdir -p /etc/apt/keyrings && curl -fsSL https://download.docker.com/linux/ubuntu/gpg | gpg --dearmor -o $DOCKER_GPG 
arch=$(dpkg --print-architecture) \
     && codename=$(awk -F= '/CODENAME/ {print $2;exit}' /etc/os-release) \
     && echo "deb [arch=$arch signed-by=$DOCKER_GPG] https://download.docker.com/linux/ubuntu $codename stable" > \
         /etc/apt/sources.list.d/docker.list \
     && apt-get update \
     && apt-get install -y --no-install-recommends \
         containerd.io \
         docker-ce \
         docker-ce-cli
 pip3 install --user jsonpickle
 cp -r Biodepot-workflow-builder/widgets /widgets
 cp -r Biodepot-workflow-builder/biodepot /biodepot
 cp -r Biodepot-workflow-builder/coreutils /coreutils

 apt-get install -y xorg 
 cp -r Biodepot-workflow-builder/scripts/generate_setup.sh /usr/local/bin/generate_setup.sh
 chmod +x /usr/local/bin/generate_setup.sh
 pip3 install -e /biodepot
 curl -L "https://github.com/docker/compose/releases/download/1.23.2/docker-compose-$(uname -s)-$(uname -m)" -o /usr/local/bin/docker-compose \
     && chmod +x /usr/local/bin/docker-compose
 sed -i 's/"Orange Canvas"/"Bwb"/' /orange3/Orange/canvas/config.py
 rm -rf ~/.fluxbox
 cp -r Biodepot-workflow-builder/fluxbox_config ~/.fluxbox
 cp -r Biodepot-workflow-builder/user_config/* ~/. 
 cp /orangePatches/schemeedit.py /orange3/Orange/canvas/document/schemeedit.py
 cp /orangePatches/canvasmain.py /orange3/Orange/canvas/application/canvasmain.py
 cp /orangePatches/widgetsscheme.py /orange3/Orange/canvas/scheme/widgetsscheme.py
 cp /orangePatches/signalmanager.py /orange3/Orange/canvas/scheme/signalmanager.py
 cp /orangePatches/link.py /orange3/Orange/canvas/scheme/link.py
 cp /orangePatches/signals.py /orange3/Orange/widgets/utils/signals.py
 cp /orangePatches/linkitem.py /orange3/Orange/canvas/canvas/items/linkitem.py
 cp /orangePatches/__main__.py /orange3/Orange/canvas/__main__.py
 cp /orangePatches/discovery.py /orange3/Orange/canvas/registry/discovery.py

cp Biodepot-workflow-builder/scripts/startBwb.sh /usr/local/bin/startBwb.sh
cp Biodepot-workflow-builder/scripts/startSingleBwb.sh /usr/local/bin/startSingleBwb.sh
cp Biodepot-workflow-builder/scripts/runDockerJob.sh /usr/local/bin/runDockerJob.sh
cp Biodepot-workflow-builder/scripts/startScheduler.sh /usr/local/bin/startScheduler.sh
cp Biodepot-workflow-builder/scripts/build_workflow_containers.sh /usr/local/bin/build_workflow_containers.sh
cp Biodepot-workflow-builder/scripts/whiteListToolDock.py /usr/local/bin/whiteListToolDock.py
cp Biodepot-workflow-builder/scripts/addWorkflowsToToolDock.py /usr/local/bin/addWorkflowsToToolDock.py
cp Biodepot-workflow-builder/scripts/addWidgetToToolDock.sh /usr/local/bin/addWidgetToToolDock.sh
cp Biodepot-workflow-builder/scripts/removeWidgetFromToolDock.sh /usr/local/bin/removeWidgetFromToolDock.sh
cp Biodepot-workflow-builder/scripts/generate_setup.sh /usr/local/bin/generate_setup.sh
cp -r Biodepot-workflow-builder/executables /usr/local/bin/executables
groupadd ftpaccess
cp Biodepot-workflow-builder/sshd_config /etc/ssh/sshd_config
cp Biodepot-workflow-builder/startSftp.sh /usr/local/bin/startSftp.sh
cp Biodepot-workflow-builder/scripts/findResolution.sh /usr/local/bin/findResolution.sh

cp -r Biodepot-workflow-builder/workflows /workflows
cp -r Biodepot-workflow-builder/notebooks /notebooks
cp -r Biodepot-workflow-builder/templates /templates
cp -r Biodepot-workflow-builder/icons /icons
cp -r Biodepot-workflow-builder/tutorialFiles /tutorialFiles
cp Biodepot-workflow-builder/serverSettings.json /biodepot

cp Biodepot-workflow-builder/startup.sh /

#cp Biodepot-workflow-builder/ config files for dev tools
cp -r Biodepot-workflow-builder/dev-files/geany/ /root/.config/
cp -r Biodepot-workflow-builder/VM/xorg.conf /etc/X11/xorg.conf
cp -r Biodepot-workflow-builder/VM/*.sh /usr/local/bin/
cp Biodepot-workflow-builder/VM/menu /root/.fluxbox/menu
cp /root/.fluxbox/bwb.svg /orange3/Orange/canvas/icons/.