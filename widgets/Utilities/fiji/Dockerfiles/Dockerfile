FROM openjdk:11.0.13-slim
RUN apt-get -y update && apt-get -y install default-jdk tar wget curl && rm -rf /var/lib/apt/lists/*
RUN wget -c https://downloads.imagej.net/fiji/archive/20201104-1356/fiji-linux64.tar.gz \
    -O - | tar -xz -C /usr/local/bin/
ENV PATH $PATH:/usr/local/bin/Fiji.app
ADD startfiji.sh /usr/local/bin/startfiji.sh
RUN ImageJ-linux64 --update update
