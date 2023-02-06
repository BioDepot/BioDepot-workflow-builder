ARG VERSION=0.48.0

FROM debian:bookworm-slim as builder
ARG VERSION

RUN apt-get update && apt-get install -y build-essential cmake autoconf automake curl \
 zlib1g-dev libhdf5-dev libbz2-dev liblzma-dev

RUN curl -L https://github.com/pachterlab/kallisto/archive/refs/tags/v${VERSION}.tar.gz | tar -zvxf - 
RUN cd kallisto*/ext/htslib && autoheader && autoconf
RUN cd kallisto* && mkdir -p build && cd build && cmake ..  && make && make install

FROM debian:bookworm-slim
MAINTAINER Ling-Hong Hung
RUN apt-get update && apt-get install -y libhdf5-10 zlib1g \
    && apt-get autoclean -y \
    && apt-get autoremove -y \
    && rm -rf /var/lib/apt/lists/*
COPY --from=builder ./kallisto /bin/kallisto
COPY runKallisto.sh /usr/local/bin/

LABEL "tag.1"="$VERSION"
LABEL "latest"="True"