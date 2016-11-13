FROM  alpine:edge
MAINTAINER Ling-Hong Hung <lhhunghimself@gmail.com>

ENV DEBIAN_FRONTEND noninteractive
ENV HOME /root
RUN apk update upgrade --no-cache && apk add --no-cache ca-certificates fluxbox bash nano musl python xvfb net-tools openssh nginx xterm shadow sudo 

COPY x11vnc/ /
RUN apk add --allow-untrusted x11vnc-0.9.13-r0.apk && rm *.apk APKINDEX.tar.gz 
ADD noVNC /noVNC
RUN sed -i -- "s/ps -p/ps -o pid | grep/g" /noVNC/utils/launch.sh

ADD web /web/
RUN apk add --no-cache gcc python2-dev musl-dev py2-pip && pip install --upgrade pip && pip install --upgrade supervisor && pip install -r /web/requirements.txt \
     && apk del gcc musl-dev py2-pip python2-dev
RUN apk add --no-cache wget && wget https://bootstrap.pypa.io/ez_setup.py -O - | python && apk del wget

ADD alpine_nginx.conf /etc/nginx/nginx.conf
ADD nginx.conf /etc/nginx/sites-enabled/default
ADD startup.sh /
ADD supervisord_base.conf /etc/supervisor/supervisord.conf
ADD supervisord2.conf /etc/supervisor/conf.d/supervisord.conf

RUN mkdir /var/log/supervisor

EXPOSE 6080
WORKDIR /root

RUN apk del *-dev
#ENTRYPOINT ["/startup.sh"]

