import argparse
import json
import logging
import pickle
import shlex
import socket
import subprocess
from multiprocessing import cpu_count as get_cpu_count
from subprocess import Popen, PIPE, run
from threading import Thread
from time import sleep
from urllib import parse
from uuid import uuid4, uuid5



import os, sys
import redis
import requests
from psutil import virtual_memory

from agent_lib import FlaskAppWrapper


def app_thread(app):
    app.run()

class DockerLocalAgent:
    API_VERSION = "1.0.0"

    def __init__(self, argv, api_url, hostport, hostname=None, cpu_count=None, mem_limit=None):
        parser = argparse.ArgumentParser(description="local-agent extra arguments")
        parser.add_argument("--redis_host", type=str, required=True)
        parser.add_argument("--redis_port", type=int, required=True)
        parsed_args = parser.parse_args(argv)
        self.__redis_host__ = parsed_args.redis_host
        self.__redis_port__ = parsed_args.redis_port

        app = FlaskAppWrapper(hostname=hostname, hostport=hostport, run_command_fnc=self.run_command,
                              status_fnc=self.status)
        self.__t__ = Thread(target=app_thread, kwargs={'app': app})

        self.__api_url__ = api_url
        self.__hostname__ = hostname or socket.gethostbyname(socket.gethostname())
        self.__hostport__ = hostport
        self.__cpu_count__ = cpu_count or get_cpu_count()
        self.__mem_limit__ = int(mem_limit or virtual_memory().total / 1024 / 1024)

    @staticmethod
    def __create_container__(job_id,worker_id,command):
        image = command['image']
        args = " ".join(command.get('args', []))
        volume_params = []
        volumes = command.get('volumes', [])
        for volume in volumes:
            host_dir = volume['host_dir']
            mount_dir = volume['mount_dir']
            mode = volume['mode']
            volume_params.append("-v %s:%s:%s" % (host_dir, mount_dir, mode))
        volume_param = " ".join(volume_params)

        env_params = []
        env_vars = command.get('env', [])
        for env_var in env_vars:
            key = env_var['key']
            val = env_var['val']
            env_params.append("-e %s=%s" % (key, val))
        env_param = " ".join(env_params)
       
        #raise RuntimeError("log_dir is {} job_id is {} worker_id is {}\n".format(log_dir,job_id,worker_id))
        #docker_command = "docker create --rm --name {}.{} alpine echo 'hello world {}' ".format(job_id,worker_id,worker_id)
        docker_command = "docker create --rm --name {}.{} {} {} {} {}".format(job_id,worker_id,env_param, volume_param, image, args)
        sys.stderr.write(docker_command)
        docker_create = run(shlex.split(docker_command), stdout=PIPE, stderr=PIPE)
        if docker_create.returncode:
            raise RuntimeError("Failed to create container\nReturn code {}\nstderr is {}\ndocker command is {}\n".format(docker_create.returncode,docker_create.stderr,docker_command))
            return 0
        sys.stderr.write('created container {}\n'.format(docker_create.stdout.decode('utf-8').rstrip()))
        return docker_create.stdout.decode('utf-8').rstrip()
    @staticmethod
    def __execute_container__(container):
        docker_command='docker start -a -i {}'.format(container)
        sys.stderr.write('{}\n'.format(docker_command))
        os.system(docker_command)
        return 0

    def fetch_and_run (self, redis_obj,command_queue_id,error_queue_id,job_id,worker_id):
        while True:
            command = redis_obj.rpop(command_queue_id)
            if command is None:
                return 0
            try:
                command_cli = json.loads(command)
                container = DockerLocalAgent.__create_container__(job_id, worker_id,command_cli)
                try:
                    ret_code=DockerLocalAgent.__execute_container__(container)
                    if ret_code:
                        redis_obj.lpush(error_queue_id,command)
                except Exception as e:
                    redis_obj.lpush(error_queue_id,command)
                    return 1
            except Exception as e:
                logging.exception(e)
                if error_queue_id:
                    redis_obj.lpush(error_queue_id,command)
                return 1

    def run_docker(self, redis_host, redis_port, command_queue_id, job_id, worker_id):
        # must connect in the function and not before
        r = redis.Redis(host=redis_host, port=redis_port)
        error_queue_id='errors.{}.{}'.format(job_id,worker_id)     
        retValue=self.fetch_and_run(r,command_queue_id,error_queue_id,job_id,worker_id)
        return retValue
        
    def run_command(self, redis_host, redis_port, command_queue_id, job_id, worker_id):
        return_code = self.run_docker(redis_host=redis_host, redis_port=redis_port, command_queue_id=command_queue_id, job_id=job_id, worker_id=worker_id)
        return json.dumps({'worker_id' : worker_id,'status' : return_code})
   
    def status(self, container_name):
        docker_command = "docker inspect %s --format='{{.State.ExitCode}}'" % container_name
        process = Popen(shlex.split(docker_command), stdin = PIPE, stdout = PIPE, stderr = PIPE, shell = False)
        out, err = process.communicate()
        if process.returncode != 0:
            logging.exception(err)
            raise RuntimeError("Failed to get docker status code")
        exit_code = out[:-1].decode()
        return exit_code

    def register_agent(self):
        endpoint = "/%s/%s" % (self.API_VERSION, "register-host")
        uri = parse.urljoin(self.__api_url__, endpoint)
        r = requests.post(uri, params={'host_name': self.__hostname__, 'host_port': self.__hostport__,
                                       'core_count': self.__cpu_count__, 'memory': self.__mem_limit__,
                                       'redis_host': self.__redis_host__, 'redis_port': self.__redis_port__})
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

