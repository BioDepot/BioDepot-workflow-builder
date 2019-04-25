import argparse
import concurrent.futures
import json
import logging
import os
import socket
from multiprocessing import cpu_count as get_cpu_count
from threading import Thread
from time import sleep
from urllib import parse

import redis
import requests
from psutil import virtual_memory

from agent_lib import FlaskAppWrapper


def app_thread(app):
    app.run()


class DockerLocalAgent:
    API_VERSION = "1.0.0"

    def __init__(
        self, argv, api_url, hostport, hostname=None, cpu_count=None, mem_limit=None
    ):
        parser = argparse.ArgumentParser(description="local-agent extra arguments")
        parser.add_argument("--redis_host", type=str, required=True)
        parser.add_argument("--redis_port", type=int, required=True)
        parsed_args = parser.parse_args(argv)
        self.__redis_host__ = parsed_args.redis_host
        self.__redis_port__ = parsed_args.redis_port

        app = FlaskAppWrapper(
            hostname=hostname, hostport=hostport, run_command_fnc=self.run_command
        )
        self.__t__ = Thread(target=app_thread, kwargs={"app": app})

        self.__api_url__ = api_url
        self.__hostname__ = hostname or socket.gethostbyname(socket.gethostname())
        self.__hostport__ = hostport
        self.__cpu_count__ = cpu_count or get_cpu_count()
        self.__mem_limit__ = int(mem_limit or virtual_memory().total / 1024 / 1024)

    @staticmethod
    def __execute_command__(command):
        image = command["image"]
        args = " ".join(command.get("args", []))

        volume_params = []
        volumes = command.get("volumes", [])
        for volume in volumes:
            host_dir = volume["host_dir"]
            mount_dir = volume["mount_dir"]
            mode = volume["mode"]
            volume_params.append("-v %s:%s:%s" % (host_dir, mount_dir, mode))
        volume_param = " ".join(volume_params)

        env_params = []
        env_vars = command.get("env", [])
        for env_var in env_vars:
            key = env_var["key"]
            val = env_var["val"]
            env_params.append("-e %s=%s" % (key, val))
        env_param = " ".join(env_params)
        docker_command = "docker run -i --rm %s %s %s %s" % (
            env_param,
            volume_param,
            image,
            args,
        )
        output = os.popen(docker_command).read()
        return output

    @staticmethod
    def run_docker(redis_host, redis_port, queue_id):
        # must connect in the function and not before
        r = redis.Redis(host=redis_host, port=redis_port)
        while True:
            queue_object = r.rpop(queue_id)
            if queue_object is None:
                break

            try:
                command = json.loads(queue_object)
                DockerLocalAgent.__execute_command__(command['command'])
            except Exception as e:
                logging.exception(e)

    def run_command(self, redis_host, redis_port, queue_id, max_workers=None):
        if not max_workers:
            max_workers = self.__cpu_count__
        with concurrent.futures.ProcessPoolExecutor(
            max_workers=max_workers
        ) as executor:
            for _ in range(max_workers):
                executor.submit(
                    DockerLocalAgent.run_docker, redis_host, redis_port, queue_id
                )
        return "Commands fed"

    def register_agent(self):
        endpoint = "/%s/%s" % (self.API_VERSION, "register-host")
        uri = parse.urljoin(self.__api_url__, endpoint)
        r = requests.post(
            uri,
            params={
                "host_name": self.__hostname__,
                "host_port": self.__hostport__,
                "core_count": self.__cpu_count__,
                "memory": self.__mem_limit__,
                "redis_host": self.__redis_host__,
                "redis_port": self.__redis_port__,
            },
        )
        if r.status_code != 200:
            raise RuntimeError("Failed to register agent", r.content)

    def start_agent(self):
        self.__t__.start()

        uri = "http://%s:%s/ping" % (self.__hostname__, self.__hostport__)
        while True:
            sleep(1)
            try:
                r = requests.get(uri)
            except Exception as e:
                logging.exception("Failed to validate agent: %s" % e)
                continue
            if r.status_code == 200:
                break

    def join(self):
        self.__t__.join()
