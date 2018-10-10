From debian:8.11-slim
MAINTAINER "Ling-Hong Hung" lhhunghimself@gmail.com 
RUN apt-get update && apt-get install -y build-essential g++ libbz2-dev libz-dev
COPY STAR-2.6.0c/ /src/
WORKDIR /src/source
RUN make clean
RUN make STARstatic
