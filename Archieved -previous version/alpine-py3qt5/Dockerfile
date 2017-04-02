FROM alpine:edge
#minimal needed to compile the entire PyQt5 package 
#PyQt5 is in /build - should be copied to the appropriate place i.e. /xx/python3.5/sites-available

RUN apk update upgrade --no-cache && apk add  build-base python3 python3-dev python2-dev py-dbus-dev qt5-qtbase-dev py3-sip py3-sip-dev libx11-dev qt5-qtsvg-dev qt5-qtwebkit-dev qt5-qttools-dev g++
RUN ln -s /usr/lib/libpython3.5m.so.1.0 /usr/lib/libpython3.5.so

COPY PyQt5_gpl-5.7.tar.gz PyQt5_gpl-5.7.tar.gz
RUN tar -xzvf  PyQt5_gpl-5.7.tar.gz && mkdir /build
WORKDIR  PyQt5_gpl-5.7
RUN python configure.py --verbose --qmake /usr/lib/qt5/bin/qmake -d  /build --confirm-license --qsci-api --target-py-version=3.5.2 --sip-incdir=/usr/include/python3.5m
RUN  make -j 8
RUN  make install
# for debugging
RUN apk add bash nano 
