ARG VERSION="2.6.0c"
From debian:bookworm-slim as builder
ENV DEBIAN_FRONTEND=noninteractive
ARG VERSION
RUN apt-get update && apt-get install -y build-essential g++ libbz2-dev libz-dev curl
RUN curl -L https://github.com/alexdobin/STAR/archive/refs/tags/${VERSION}.tar.gz |  tar -zvxf -
RUN cd STAR*/source && make clean && make -j 8 STARstatic
#new STAR binary in source directory
From debian:bookworm-slim
ARG VERSION
ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get update && apt-get install -y gzip bzip2 \
    && apt-get autoclean -y \
    && apt-get autoremove -y \
    && rm -rf /var/lib/apt/lists/*
COPY --from=builder /STAR*/source/STAR /usr/local/bin/STAR
COPY runstar.sh /usr/local/bin/runstar.sh
ENV PATH="$PATH:/usr/local/bin" 
LABEL "tag1"="$VERSION"
