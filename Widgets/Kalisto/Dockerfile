FROM ubuntu:14.04

MAINTAINER xlianguw@uw.edu
RUN echo 'debconf debconf/frontend select Noninteractive' | debconf-set-selections

RUN apt-get update -y
RUN apt-get install build-essential -y
RUN apt-get install git -y
RUN apt-get install cmake -y
RUN apt-get install libc-dev -y
RUN apt-get install zlib1g-dev -y
RUN apt-get install libhdf5-dev -y
RUN git clone https://github.com/pachterlab/kallisto.git
WORKDIR /kallisto
RUN mkdir build
WORKDIR build
RUN cmake ..
RUN make
RUN make install

ADD test.sh /kallisto/test

CMD ["/bin/bash"]
