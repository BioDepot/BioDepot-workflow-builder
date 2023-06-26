FROM python:3.8.0-alpine3.10
ARG CLI_VERSION=1.16.272
RUN apk upgrade && apk add --no-cache nano bash zip
RUN pip install --upgrade pip && pip install --no-cache-dir awscli==$CLI_VERSION && \
    mkdir /root/.aws
COPY upload.sh /usr/local/bin/upload.sh
