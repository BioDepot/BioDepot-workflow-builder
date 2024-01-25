FROM debian:buster-slim as builder
RUN apt-get update && apt-get install -y gcc libc-dev libxext-dev
COPY docker_xnoshm.c .
RUN gcc docker_xnoshm.c -shared -o docker_xnoshm.so
FROM debian:buster-slim
RUN apt-get -y update && apt-get -y install gnumeric && apt-get autoclean -y \
    && apt-get autoremove -y \
    && rm -rf /var/lib/apt/lists/*
COPY --from=builder docker_xnoshm.so /usr/local/lib/docker_xnoshm.so
ENV LD_PRELOAD=/usr/local/lib/docker_xnoshm.so
