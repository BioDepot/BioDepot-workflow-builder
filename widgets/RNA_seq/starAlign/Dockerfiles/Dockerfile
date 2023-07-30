ARG VERSION="2.6.0c"
From debian:bookworm-slim as builder
ENV DEBIAN_FRONTEND=noninteractive
ARG VERSION
ARG TARGETARCH
RUN apt-get update && apt-get install -y build-essential g++ libbz2-dev libz-dev curl xxd
RUN curl -L https://github.com/alexdobin/STAR/archive/refs/tags/${VERSION}.tar.gz |  tar -zvxf -
RUN if [ "$TARGETARCH" = "arm64" ] ;\
    then cd STAR*/source && make clean && make -j 8  CXXFLAGS_SIMD='-march=armv8.5-a+simd' STARstatic; \
    else  cd STAR*/source && make clean && make -j 8 STARstatic; \
    fi 
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
