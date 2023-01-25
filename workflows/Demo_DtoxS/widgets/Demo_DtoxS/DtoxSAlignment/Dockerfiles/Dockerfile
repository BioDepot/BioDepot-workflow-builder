ARG VERSION=0.7.15

FROM debian:buster-slim as builder
ARG VERSION

RUN apt-get update && apt-get install -y curl make build-essential zlib1g-dev 
RUN curl -L https://github.com/lh3/bwa/archive/refs/tags/v${VERSION}.tar.gz | tar -zxvf - 
RUN cd bwa* && make -j 8

FROM debian:buster-slim
ARG VERSION
ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get update && apt-get -y install python python-numpy \
 && apt-get autoclean -y \
 && apt-get autoremove -y \
 && rm -rf /var/lib/apt/lists/*
COPY --from=builder /bwa-${VERSION}/bwa /usr/local/bwa
COPY Programs /root/Programs

ENV PATH="$PATH:/usr/local/bin"

LABEL "tag.0"="$VERSION"
LABEL "latest"="True"


