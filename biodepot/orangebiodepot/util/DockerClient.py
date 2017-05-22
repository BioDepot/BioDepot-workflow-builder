# -*- coding: utf-8 -*-
from docker import APIClient
import requests, json
from PyQt5.QtCore import QThread, pyqtSignal
import os


class DockerClient:
    def __init__(self, url, name):
        self.url = url
        self.name = name
        self.cli = APIClient(base_url=url)

    def getClient(self):
        return self.cli

    def getName(self):
        return self.name

    def getUrl(self):
        return self.url

    def images(self):
        return self.cli.images(all=True)

    def has_image(self, name, version):
        repo_tag = name + ':' + version
        for image in self.cli.images():
            if not image['RepoTags']: continue #DK fix NoneType is not iterable
            elif repo_tag in image['RepoTags']: return True
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
    def create_container(self, name, volumes=None, commands=None, environment=None):
        # TODO should we use Image ID instead of Image Name?
        host_config = None
        if type(volumes) is dict:
            binds = []
            for host_dir, container_dir  in volumes.items():
                binds.append(self.to_host_directory(host_dir) + ":" + container_dir)
            host_config = self.cli.create_host_config(binds=binds)
            volumes = list(volumes.values())
        if type(commands) is list:
            commands = "bash -c \"" + ' && '.join(commands) + "\""
        return self.cli.create_container(image=name,
                                         volumes=volumes,
                                         command=commands,
                                         environment=environment,
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

    '''
        Convert path of container back to host path
    '''
    def to_host_directory(self, path):
        source = ''
        destination = ''

        # locate BwB container
        for c in self.cli.containers():
            if c['Image'] == 'biodepot/bwb':
                # found BwB container, locate source and destination
                for m in c['Mounts']:
                    if 'docker.sock' in m['Source']:
                        continue
                    source = m['Source']
                    destination = m['Destination']

        if source is '' or destination is '':
            return path

        destination = os.path.join(destination, '')

        # if the path is not mapping from host, nothing will be done
        if destination not in path :
            return path

        abspath = os.path.join(source, path[path.find(destination)+len(destination) : ])
        return abspath


class DockerThread_BuildImage(QThread):
    build_process = pyqtSignal(['QString'])
    build_complete = pyqtSignal(int)
    def __init__(self, cli, name, dockerfile):
        QThread.__init__(self)
        self.docker = cli
        self.dockerfile = dockerfile
        self.name = name

    def __del__(self):
        self.wait()

    def run(self):
        imagename = self.name
        try:
            import io
            f = io.BytesIO(self.dockerfile.encode('utf-8'))
            for rawline in self.docker.getClient().build(fileobj=f, tag=imagename, stream=True):
                for jsonstr in rawline.decode('utf-8').split('\r\n')[:-1]:
                    log = jsonstr
                    try:
                        line = json.loads(jsonstr)
                        log = line['stream']
                    except ValueError as e:
                        print (e)
                    except TypeError as e:
                        print (e)
                    self.build_process.emit(log)

        except requests.exceptions.RequestException as e:
            # TODO emit error
            print(e)

class PullImageThread(QThread):
    pull_progress = pyqtSignal(int)
    def __init__(self, cli, name, version):
        QThread.__init__(self)
        self.docker = cli
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
            for line in self.docker.getClient().pull(repo_tag, stream=True):
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
                        self.current_progress = sum(progs.values()) / len(progs)
                        self.pull_progress.emit(self.current_progress)


        except requests.exceptions.RequestException:
            # TODO emit error
            print('Connection Exception!')