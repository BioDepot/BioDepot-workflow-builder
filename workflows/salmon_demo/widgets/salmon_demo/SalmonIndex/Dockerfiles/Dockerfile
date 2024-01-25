ARG VERSION=0.48.0

FROM debian:buster-slim as builder
ARG VERSION

RUN apt-get update && apt-get install -y build-essential cmake autoconf automake curl \
 zlib1g-dev libhdf5-dev libbz2-dev liblzma-dev

RUN curl -L https://github.com/pachterlab/kallisto/archive/refs/tags/v${VERSION}.tar.gz | tar -zvxf - 
RUN cd kallisto*/ext/htslib && autoheader && autoconf
RUN cd kallisto* && mkdir -p build && cd build && cmake .. -DUSE_HDF5=ON && make && make install

FROM debian:buster-slim
MAINTAINER Ling-Hong Hung
RUN apt-get update && apt-get install -y libhdf5-dev zlib1g-dev libbz2-dev liblzma-dev \
    && apt-get autoclean -y \
    && apt-get autoremove -y \
    && rm -rf /var/lib/apt/lists/*
COPY --from=builder /usr/local/bin/kallisto /usr/local/bin/kallisto
COPY runKallisto.sh /usr/local/bin/

LABEL "tag.1"="$VERSION"
LABEL "latest"="True"

#Note that bookworm uses autoconf 2.71 which is buggy and causes problems with the htslib compile