# -*- coding: utf-8 -*-
from docker import Client
import requests, json
from PyQt4.QtCore import QThread, SIGNAL


class DockerClient:
    def __init__(self, url, name):
        self.url = url
        self.name = name
        self.cli = Client(base_url=url)

    def getName(self):
        return self.name

    def getUrl(self):
        return self.url

    def images(self):
        return self.cli.images(all=True)

    def has_image(self, name, version):
        repo_tag = name + ':' + version
        for image in self.cli.images():
            if repo_tag in image['RepoTags']:
                return True
        return False

    def remove_image(self, id, force=False):
        self.cli.remove_image(id, force=force)

    def containers(self, all=True):
        return self.cli.containers(all=all)

    """
    volumes is a dict mapping host directory to container directory
    {
        "/Users/host/directory": "path/to/container/directory"
    }
    commands is a list of bash commands to run on container
        ["pwd", "touch newfile.txt"]
    """
    def create_container(self, name, volumes=None, commands=None):
        # TODO should we use Image ID instead of Image Name?
        host_config = None
        if type(volumes) is dict:
            binds = []
            for host_dir, container_dir  in volumes.items():
                binds.append(host_dir + ":" + container_dir)
            host_config = self.cli.create_host_config(binds=binds)
            volumes = list(volumes.values())
        if type(commands) is list:
            commands = "bash -c \"" + ' && '.join(commands) + "\""
        return self.cli.create_container(image=name,
                                         volumes=volumes,
                                         command=commands,
                                         stdin_open=True,
                                         host_config=host_config)

    def start_container(self, id):
        return self.cli.start(id)

    def container_running(self, id):
        for container in self.containers(all=False):
            if container['Id'] == id:
                return True
        return False

    def remove_container(self, id, force=False):
        self.cli.remove_container(id, force=force)

    def stop_container(self, id):
        self.cli.stop(id)

    def pause_container(self, id):
        self.cli.pause(id)

    def unpause_container(self, id):
        self.cli.unpause(id)

    def version(self):
        return self.cli.version()

    def info(self):
        return self.cli.info()

    def volumes(self):
        return self.cli.volumes()['Volumes']

    def remove_volume(self, name):
        self.cli.remove_volume(name)


class PullImageThread(QThread):
    def __init__(self, cli, name, version):
        QThread.__init__(self)
        self.cli = cli
        self.name = name
        self.version = version

    def __del__(self):
        self.wait()

    def run(self):
        repo_tag = self.name + ':' + self.version
        try:
            # Docker splits downloads into multiple parts
            # We create a dict mapping id to % finished for each part (progress)
            # The total progress is the mean of individual progresses
            progs = dict()
            for line in self.cli.pull(repo_tag, stream=True):
                for status in line.decode('utf-8').split('\r\n')[:-1]:
                    line = json.loads(status)
                    statusStr = line['status']
                    if statusStr == 'Pulling fs layer':
                        progs[line['id']] = 0
                    # First 50% progress is Downloading
                    elif statusStr == 'Downloading':
                        progDetail = line['progressDetail']
                        if len(progDetail) > 1:
                            progs[line['id']] = progDetail['current'] / progDetail['total'] * 50
                    # Last 50% progress is Extracting
                    elif statusStr == 'Extracting':
                        progDetail = line['progressDetail']
                        if len(progDetail) > 1:
                            progs[line['id']] = 50 + (progDetail['current'] / progDetail['total'] * 50)
                    if (len(progs) > 0):
                        self.emit(SIGNAL('pull_progress'), sum(progs.values()) / len(progs))

        except requests.exceptions.RequestException:
            print('Connection Exception!')