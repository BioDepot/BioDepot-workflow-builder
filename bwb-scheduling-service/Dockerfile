FROM python:3.7.2-slim-stretch

RUN mkdir -p /usr/src/app
WORKDIR /usr/src/app

COPY requirements.txt /usr/src/app/

RUN apt-get update && apt-get install -y build-essential \
    && pip3 install --no-cache-dir -r requirements.txt \
    && apt-get remove -y build-essential \
    && apt-get autoclean -y \
    && apt-get autoremove -y \
    && rm -rf /var/lib/apt/lists/*

COPY . /usr/src/app

EXPOSE 8080

ENTRYPOINT ["/usr/src/app/docker_init.sh"]
