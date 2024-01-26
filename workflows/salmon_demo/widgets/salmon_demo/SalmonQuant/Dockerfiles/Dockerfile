ARG VERSION="1.10.1"
From debian:bookworm-slim as builder
ENV DEBIAN_FRONTEND=noninteractive
ARG VERSION
RUN apt-get update && apt-get install -y build-essential g++ libbz2-dev libz-dev curl libboost-all-dev liblzma-dev ca-certificates zlib1g-dev libcurl4-openssl-dev curl unzip autoconf apt-transport-https ca-certificates gnupg software-properties-common wget cmake
RUN curl -L https://github.com/COMBINE-lab/salmon/archive/refs/tags/v${VERSION}.tar.gz |  tar -zvxf -
RUN  cd salmon-${VERSION} && \
    mkdir build && \
    cd build && \
    cmake .. -DCMAKE_INSTALL_PREFIX=/usr/local/salmon && make && make install
From debian:bookworm-slim
ARG VERSION
ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get update && apt-get install -y libhwloc15 \
    && apt-get autoclean -y \
    && apt-get autoremove -y \
    && rm -rf /var/lib/apt/lists/*
COPY --from=builder /usr/local/salmon/ /usr/local/
COPY runSalmon.sh /usr/local/bin/runSalmon.sh
