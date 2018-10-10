FROM alpine:3.7
RUN apk update upgrade --no-cache && apk add bash curl wget unzip tar bzip2\
 && rm -rf /var/cache/apk*
COPY download.sh /root
