#!/bin/bash
#cd /home/$USERNAME
#git clone https://github.com/BioDepot/BioDepot-workflow-builder.git
#Then run this script using 
# /home/$USERNAME/BioDepot-workflow-builder/VM/install_ubuntu_22.04.sh
set -e
USERNAME=ubuntu
DOCKER_GPG=/etc/apt/keyrings/docker.gpg


BIODEPOT=/home/$USERNAME/BioDepot-workflow-builder
export DEBIAN_FRONTEND=noninteractive
echo "package_name package_name/option select value" |  debconf-set-selections
 apt-get -y update
 apt-get install -y --no-install-recommends \
        curl \
        build-essential \
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
        python3.12-dev\
         python3.12-venv\
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
python3.12 -m venv /venv
source /venv/bin/activate
python3.12 -m ensurepip --upgrade
python3.12 -m pip install --upgrade pip setuptools 
  rsync -av  $BIODEPOT/VM/orange3/ /orange3/
 pip3 install -r /orange3/requirements-gui.txt \
         beautifulsoup4 \
         docker \
         pysam
 pip3 install -e /orange3

 apt-get install -y --no-install-recommends gnupg 
 mkdir -p /etc/apt/keyrings &&  curl -fsSL https://download.docker.com/linux/ubuntu/gpg |  gpg --dearmor -o $DOCKER_GPG 
arch=$(dpkg --print-architecture) \
     && codename=$(awk -F= '/CODENAME/ {print $2;exit}' /etc/os-release) \
     &&  echo "deb [arch=$arch signed-by=$DOCKER_GPG] https://download.docker.com/linux/ubuntu $codename stable" > \
         /etc/apt/sources.list.d/docker.list \
     &&  apt-get update \
     &&  apt-get install -y --no-install-recommends \
         containerd.io \
         docker-ce \
         docker-ce-cli
  rsync -av $BIODEPOT/widgets/ /widgets/ \
  &&  rsync -av $BIODEPOT/biodepot/ /biodepot/ \
  &&  rsync -av $BIODEPOT/VM/coreutils/ /coreutils/

  apt-get install -y xorg 
  cp  $BIODEPOT/scripts/generate_setup.sh /usr/local/bin/generate_setup.sh
  chmod +x /usr/local/bin/generate_setup.sh
  pip3 install -e /biodepot
  curl -L "https://github.com/docker/compose/releases/download/1.23.2/docker-compose-$(uname -s)-$(uname -m)" -o /usr/local/bin/docker-compose \
     && chmod +x /usr/local/bin/docker-compose
  sed -i 's/"Orange Canvas"/"Bwb"/' /orange3/Orange/canvas/config.py
  rm -rf ~/.fluxbox
  rm -rf ~/.config/biolab.si
  rsync -av $BIODEPOT/VM/user_config/ ~/ 

 cp $BIODEPOT/scripts/startSingleBwb.sh /usr/local/bin/startSingleBwb.sh
 cp $BIODEPOT/scripts/runDockerJob.sh /usr/local/bin/runDockerJob.sh
 cp $BIODEPOT/scripts/startScheduler.sh /usr/local/bin/startScheduler.sh
 cp $BIODEPOT/scripts/build_workflow_containers.sh /usr/local/bin/build_workflow_containers.sh
 cp $BIODEPOT/scripts/whiteListToolDock.py /usr/local/bin/whiteListToolDock.py
 cp $BIODEPOT/scripts/addWorkflowsToToolDock.py /usr/local/bin/addWorkflowsToToolDock.py
 cp $BIODEPOT/scripts/addWidgetToToolDock.sh /usr/local/bin/addWidgetToToolDock.sh
 cp $BIODEPOT/scripts/removeWidgetFromToolDock.sh /usr/local/bin/removeWidgetFromToolDock.sh
 cp $BIODEPOT/scripts/generate_setup.sh /usr/local/bin/generate_setup.sh
 cp -r $BIODEPOT/executables /usr/local/bin/executables
 groupadd ftpaccess
 cp $BIODEPOT/sshd_config /etc/ssh/sshd_config
 cp $BIODEPOT/startSftp.sh /usr/local/bin/startSftp.sh
 cp $BIODEPOT/scripts/findResolution.sh /usr/local/bin/findResolution.sh

 cp -r $BIODEPOT/workflows /workflows
 cp -r $BIODEPOT/notebooks /notebooks
 cp -r $BIODEPOT/templates /templates
 cp -r $BIODEPOT/icons /icons
 cp -r $BIODEPOT/tutorialFiles /tutorialFiles
 cp $BIODEPOT/serverSettings.json /biodepot

 cp $BIODEPOT/startup.sh /

# cp $BIODEPOT/ config files for dev tools
 cp -r $BIODEPOT/dev-files/geany/ /root/.config/
 cp -r $BIODEPOT/VM/xorg.conf /etc/X11/xorg.conf
 cp -r $BIODEPOT/VM/*.sh /usr/local/bin/
 cp $BIODEPOT/VM/menu /root/.fluxbox/menu
 cp /root/.fluxbox/bwb.svg /orange3/Orange/canvas/icons/orange-canvas.svg
# cp -r ~/biolab.si ~/.config/
#for the user
# groupadd $USERNAME
 usermod -aG docker $USERNAME
 chown -R $USERNAME:$USERNAME /biodepot $BIODEPOT /orange3 /widgets /workflows /coreutils  /icons 
 rsync -av $BIODEPOT/VM/user_config/ /home/$USERNAME/ && chown -R $USERNAME:$USERNAME /home/$USERNAME
