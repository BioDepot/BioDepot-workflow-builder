FROM debian:jessie-slim
RUN apt-get update && apt-get -y install gnumeric \
    && apt-get autoclean -y \
    && apt-get autoremove -y \
    && rm -rf /var/lib/apt/lists/*
