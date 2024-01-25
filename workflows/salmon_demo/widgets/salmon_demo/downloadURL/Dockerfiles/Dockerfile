FROM alpine:3.15
RUN adduser --disabled-password -s /bin/bash user
RUN apk update upgrade --no-cache && apk add bash curl wget unzip tar bzip2 grep \
 && rm -rf /var/cache/apk*
COPY download.sh /root/.
