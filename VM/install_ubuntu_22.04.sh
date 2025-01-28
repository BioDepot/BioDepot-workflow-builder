#!/bin/bash
#run this as root from /
set -e
USERNAME=ubuntu
DOCKER_GPG=/etc/apt/keyrings/docker.gpg
BIODEPOT=/home/$USERNAME/BioDepot-workflow-builder
export DEBIAN_FRONTEND=noninteractive
sudo apt-get -y update
sudo apt-get install -y --no-install-recommends \
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
        supervisor \
        wget \
        xdg-utils \
        lxterminal \
        x11vnc \
        xvfb \
        zlib1g-dev \
        zenity

  #cd / &&  git clone https://github.com/BioDepot/$BIODEPOT.git 
  #cp -r $BIODEPOT/orange3 /orange3
  #cp -r $BIODEPOT/orangePatches /orangePatches
sudo apt-get install -y \
        build-essential \
        python3-dev \
        python3-pip \
    && python3 -m pip install --upgrade \
        pip==20.0.1 \
        setuptools \
        wheel 
sudo  rsync -a  $BIODEPOT/VM/orange3/ /orange3/ && sudo rsync -a $BIODEPOT/orangePatches/ /orangePatches/
sudo pip3 install -r /orange3/requirements-gui.txt \
         beautifulsoup4 \
         docker \
         pysam
sudo pip3 install setuptools==65.5.1
sudo pip3 install -e /orange3

sudo apt-get install -y --no-install-recommends gnupg 
sudo mkdir -p /etc/apt/keyrings && sudo curl -fsSL https://download.docker.com/linux/ubuntu/gpg | sudo gpg --dearmor -o $DOCKER_GPG 
arch=$(dpkg --print-architecture) \
     && codename=$(awk -F= '/CODENAME/ {print $2;exit}' /etc/os-release) \
     && sudo echo "deb [arch=$arch signed-by=$DOCKER_GPG] https://download.docker.com/linux/ubuntu $codename stable" > \
         /etc/apt/sources.list.d/docker.list \
     && sudo apt-get update \
     && sudo apt-get install -y --no-install-recommends \
         containerd.io \
         docker-ce \
         docker-ce-cli
 sudo pip3 install --user jsonpickle
 sudo cp -r $BIODEPOT/widgets /widgets \
  && sudo cp -r $BIODEPOT/biodepot /biodepot \
  && sudo cp -r $BIODEPOT/VM/coreutils /coreutils

 sudo apt-get install -y xorg 
 sudo cp -r $BIODEPOT/scripts/generate_setup.sh /usr/local/bin/generate_setup.sh
 sudo chmod +x /usr/local/bin/generate_setup.sh
 sudo pip3 install -e /biodepot
 sudo curl -L "https://github.com/docker/compose/releases/download/1.23.2/docker-compose-$(uname -s)-$(uname -m)" -o /usr/local/bin/docker-compose \
     && chmod +x /usr/local/bin/docker-compose
 sudo sed -i 's/"Orange Canvas"/"Bwb"/' /orange3/Orange/canvas/config.py
 sudo rm -rf ~/.fluxbox
 sudo rm -rf ~/.config/biolab.si
 sudo cp -r $BIODEPOT/fluxbox_config ~/.fluxbox \
 && cp -r $BIODEPOT/user_config/* ~/. \
 && cp /orangePatches/schemeedit.py /orange3/Orange/canvas/document/schemeedit.py \
 && cp /orangePatches/canvasmain.py /orange3/Orange/canvas/application/canvasmain.py \
 && cp /orangePatches/widgetsscheme.py /orange3/Orange/canvas/scheme/widgetsscheme.py \
 && cp /orangePatches/signalmanager.py /orange3/Orange/canvas/scheme/signalmanager.py \
 && cp /orangePatches/link.py /orange3/Orange/canvas/scheme/link.py \
 && cp /orangePatches/signals.py /orange3/Orange/widgets/utils/signals.py \
 && cp /orangePatches/linkitem.py /orange3/Orange/canvas/canvas/items/linkitem.py \
 && cp /orangePatches/__main__.py /orange3/Orange/canvas/__main__.py \
 && cp /orangePatches/discovery.py /orange3/Orange/canvas/registry/discovery.py 

sudo cp $BIODEPOT/scripts/startBwb.sh /usr/local/bin/startBwb.sh
sudo cp $BIODEPOT/scripts/startSingleBwb.sh /usr/local/bin/startSingleBwb.sh
sudo cp $BIODEPOT/scripts/runDockerJob.sh /usr/local/bin/runDockerJob.sh
sudo cp $BIODEPOT/scripts/startScheduler.sh /usr/local/bin/startScheduler.sh
sudo cp $BIODEPOT/scripts/build_workflow_containers.sh /usr/local/bin/build_workflow_containers.sh
sudo cp $BIODEPOT/scripts/whiteListToolDock.py /usr/local/bin/whiteListToolDock.py
sudo cp $BIODEPOT/scripts/addWorkflowsToToolDock.py /usr/local/bin/addWorkflowsToToolDock.py
sudo cp $BIODEPOT/scripts/addWidgetToToolDock.sh /usr/local/bin/addWidgetToToolDock.sh
sudo cp $BIODEPOT/scripts/removeWidgetFromToolDock.sh /usr/local/bin/removeWidgetFromToolDock.sh
sudo cp $BIODEPOT/scripts/generate_setup.sh /usr/local/bin/generate_setup.sh
sudo cp -r $BIODEPOT/executables /usr/local/bin/executables
sudo groupadd ftpaccess
sudo cp $BIODEPOT/sshd_config /etc/ssh/sshd_config
sudo cp $BIODEPOT/startSftp.sh /usr/local/bin/startSftp.sh
sudo cp $BIODEPOT/scripts/findResolution.sh /usr/local/bin/findResolution.sh

sudo cp -r $BIODEPOT/workflows /workflows
sudo cp -r $BIODEPOT/notebooks /notebooks
sudo cp -r $BIODEPOT/templates /templates
sudo cp -r $BIODEPOT/icons /icons
sudo cp -r $BIODEPOT/tutorialFiles /tutorialFiles
sudo cp $BIODEPOT/serverSettings.json /biodepot

sudo cp $BIODEPOT/startup.sh /

#sudo cp $BIODEPOT/ config files for dev tools
sudo cp -r $BIODEPOT/dev-files/geany/ /root/.config/
sudo cp -r $BIODEPOT/VM/xorg.conf /etc/X11/xorg.conf
sudo cp -r $BIODEPOT/VM/*.sh /usr/local/bin/
sudo cp $BIODEPOT/VM/menu /root/.fluxbox/menu
sudo cp /root/.fluxbox/bwb.svg /orange3/Orange/canvas/icons/orange-canvas.svg
sudo cp -r ~/biolab.si ~/.config/
#for the user
sudo usermod -aG docker $USERNAME
sudo chown -R $USERNAME:$USERNAME /biodepot /$BIODEPOT /orange3 /widgets /workflows /coreutils /orangePatches /icons 
sudo rsync -av /$BIODEPOT/VM/user_config/ /home/$USERNAME/ && chown -R $USERNAME:$USERNAME /home/$USERNAME


#log out and login as $USERNAME 
#supervisorctl -c ~/supervisor/supervisord.conf reread
#supervisorctl -c ~/supervisor/supervisord.conf update
#supervisorctl -c ~/supervisor/supervisord.conf start all
