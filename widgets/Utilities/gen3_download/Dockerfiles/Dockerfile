FROM alpine:3.12 as builder
ARG GEN3_CLIENT_VERSION=2021.03

RUN apk add --update --no-cache binutils \
    && wget https://github.com/uc-cdis/cdis-data-client/releases/download/${GEN3_CLIENT_VERSION}/dataclient_linux.zip \
    && unzip *zip \
    && strip gen3-client

# Store only the resulting binary in the final image
# Resulting in significantly smaller docker image size
FROM alpine:3.12
MAINTAINER Hong Hung
RUN adduser --disabled-password -s /bin/bash user
RUN apk add --update --no-cache bash ca-certificates curl
# Create a symlink to make the precompiled utility work with Alpine
RUN ln -s lib lib64 \
    && cd lib \
    && ln -s ld-musl-x86_64.so.1 ld-linux-x86-64.so.2
COPY --from=builder /gen3-client /usr/local/bin/
COPY runG3Download.sh /usr/local/bin/runG3Download.sh
