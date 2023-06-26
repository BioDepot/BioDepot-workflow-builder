FROM ubuntu:20.04
LABEL maintainer="Bob Schmitz<schmitzr1984@gmail.com>"
ARG IGV_VERSION="2.9.2"
# This flag is needed to render Java in Bwb
ENV JAVA_TOOL_OPTIONS="-Dsun.java2d.xrender=false"
RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        bash \
        openjdk-11-jre \
        unzip \
        wget \
    && wget -nv https://data.broadinstitute.org/igv/projects/downloads/$(echo ${IGV_VERSION} | sed 's|.[0-9]\+$||')/IGV_Linux_${IGV_VERSION}_WithJava.zip \
    && unzip -q *.zip -x "*/jdk-11/*" \
    && rm *.zip \
# Remove some benign warnings
    && sed -i 's|--add-exports java.desktop/com.sun.java.swing.plaf.windows=jide.common ||; \
               s|--add-exports java.desktop/sun.awt.windows=jide.common ||' IGV*/igv.args \
    && mv IGV*/* /usr/local/bin \
    && rmdir IGV* \
    && apt-get purge --auto-remove -y \
        unzip \
        wget \
    && rm -rf /var/lib/apt/lists/*
COPY start-igv.sh /usr/local/bin
