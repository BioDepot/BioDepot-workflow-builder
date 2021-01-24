FROM alpine:3.12.0
RUN adduser --disabled-password -s /bin/bash  user
RUN apk update upgrade --no-cache && apk add bash curl wget unzip tar bzip2\
 && rm -rf /var/cache/apk*
COPY download.sh /root/.
