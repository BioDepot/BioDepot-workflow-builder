FROM debian:buster-slim as builder
RUN apt-get update && apt-get -y install git tar gzip
RUN git clone https://github.com/smikkelsendk/bwa-for-arm.git
RUN tar -xzvf  bwa-for-arm/bin/bwa0717-aarch64.tar.gz
RUN tar -xzvf  bwa-for-arm/bin/bwa0717-x64.tar.gz
RUN mv bwa0717-aarch64 /usr/local/bin/bwa-arm64
RUN mv bwa0717-x64 /usr/local/bin/bwa-amd64

FROM debian:buster-slim
RUN apt-get update && apt-get -y install python python-numpy
COPY --from=builder /usr/local/bin /usr/local/bin
RUN ln  -s /usr/local/bin/bwa-$(dpkg --print-architecture) /usr/local/bin/bwa
COPY Programs /root/Programs

