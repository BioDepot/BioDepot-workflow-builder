import os, sys
import json
import requests
import subprocess
from docker import APIClient
from PyQt5.QtCore import QThread, pyqtSignal,QProcess, Qt
from PyQt5 import QtWidgets, QtGui, QtCore
import socket
import datetime

class ConsoleProcess():
    #subclass that attaches a process and pipes the output to textedit widget console widget
    def __init__(self, console=None, errorHandler=None, finishHandler=None):
        self.process=QProcess()
        self.console=console
        self.cidFile=""
        self.state='stopped'
        if console:
            self.process.readyReadStandardOutput.connect(lambda: self.writeConsole(self.process, console,self.process.readAllStandardOutput,Qt.white))
            self.process.readyReadStandardError.connect(lambda: self.writeConsole(self.process, console,self.process.readAllStandardError,Qt.red))
        if finishHandler:
            self.process.finished.connect(finishHandler)
            
    def writeConsole(self,process,console,read,color):
        console.setTextColor(color)
        console.append(read().data().decode('utf-8',errors="ignore"))
        #console.append('</span>')
    
    def writeMessage(self,message,color=Qt.green):
        #for bwb messages
        self.console.setTextColor(color)
        self.console.append(message)
        
    def stop(self,message=None):
        self.state='stopped'
        #we need to write our own code to
        try: 
            with open (self.cidFile,'r') as f:
                cid=f.read()
        except Exception as e:
            sys.stderr.write("unable to read cidFile\n")
            return
        stopCmd='docker stop {} '.format(cid)
        sys.stderr.write('Stop command: {}'.format(stopCmd))
        os.system(stopCmd)
        self.process.kill()
        self.cleanup()
        if message:
            self.writeMessage(message)
            
    def cleanup(self):
        os.unlink(self.cidFile)

        
class DockerClient:
    def __init__(self, url, name):
        self.url = url
        self.name = name
        self.cli = APIClient(base_url=url)
        self.bwb_instance_id = socket.gethostname()
        self.bwbMounts={}
        self.findVolumeMappings()
        self.logFile=None

    def getClient(self):
        return self.cli

    def getName(self):
        return self.name

    def getUrl(self):
        return self.url

    def images(self):
        return self.cli.images(all=True)

    def has_image(self, name, version="latest"):
        if not name:
            return False
        repoTag = name + ':' + version
        conId=subprocess.check_output(['docker', 'images', '-q', repoTag])
        if conId:
            return True
        return False

    def remove_image(self, id, force=False):
        self.cli.remove_image(id, force=force)
        
    def pull_image(self, id):
        self.cli.pull(id)
        
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
    def create_container_cli(self, name, volumes=None, commands=None, environment=None, hostVolumes=None, consoleProc=None, exportGraphics=False, portMappings=None,testMode=False,logFile=None):
        #reset logFile when it is not None - can be "" though - this allows an active reset
        if logFile is not None:
            self.logFile = logFile
            
        #skips DockerPy and creates the command line equivalent
        volumeMappings=''
        for container_dir, host_dir in hostVolumes.items():
            volumeMappings=volumeMappings+"-v {}:{} ".format(self.to_best_host_directory(host_dir),container_dir)
        envs=''
        for env, var in environment.items():
            #strip whitespace
            env.strip()
            #strip quotes if present
            if env[0] == env[-1] and env.startswith(("'",'"')):
                env=env[1:-1]
            envs=envs+ "-e {}={} ".format(env,var)
        #create container
        consoleProc.cidFile='/tmp/'+ str(datetime.datetime.now().date()) + '_' + str(datetime.datetime.now().time()).replace(':', '.')
        dockerBaseCmd='docker run -i --rm '
        if exportGraphics:
            dockerBaseCmd+='-e DISPLAY=:1 -v /tmp/.X11-unix:/tmp/.X11-unix '
        if portMappings:
            for mapping in portMappings:
                dockerBaseCmd+='-p {} '.format(mapping)
        dockerCmd=dockerBaseCmd + ' --init --cidfile={} {} {} {} {}'.format(consoleProc.cidFile,volumeMappings,envs,name,commands)
        sys.stderr.write('Docker command is\n{}\n'.format(dockerCmd))
        consoleProc.state='running'
        if testMode:
            dockerCmd=dockerBaseCmd + ' --init {} {} {} {}'.format(volumeMappings,envs,name,commands)
            #Do not test for logFile - this may be None if it is not the first widget in testMode
            if self.logFile:
                with open (self.logFile,'a') as f:
                    f.write('    {}\n'.format(dockerCmd))
            consoleProc.process.start('echo',[dockerCmd])
        else:   
            consoleProc.process.start('/bin/bash',['-c',dockerCmd])

        
    def create_container(self, name, volumes=None, commands=None, environment=None, hostVolumes=None):
        #hostVolues is a dict with keys being the container volumes
        # TODO should we use Image ID instead of Image Name?
        host_config = None
  
        if not (hostVolumes is None):
            binds = []
            for container_dir, host_dir in hostVolumes.items():
                binds.append(self.to_best_host_directory(host_dir) + ":" + container_dir)
            host_config = self.cli.create_host_config(binds=binds)
            volumes = list(hostVolumes.keys())           
        elif type(volumes) is dict:
        # this is backwards - it is possible to have the same host directory mapped to multiple containers but not the other way
        # keep this so as not to break early widgets
            binds = []
            for host_dir, container_dir in volumes.items():
                binds.append(self.to_best_host_directory(host_dir) + ":" + container_dir)
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
    
    def findVolumeMappings(self):
        for c in self.cli.containers():
            container_id = c['Id']
            if len(container_id) < 12: 
                continue
            if container_id[:12] == self.bwb_instance_id:
                for m in c['Mounts']:
                    if not ('/var/run' in m['Source']):
                        self.bwbMounts[m['Source']]=m['Destination']

    def to_best_host_directory(self, path, returnNone=False):
        if self.bwbMounts == {}:
            return path
        bestPath = None
        for source, dest in self.bwbMounts.items():
            absPath=self.to_host_directory(path, source, dest)
            if absPath is not None:
                if bestPath is None:
                    bestPath=absPath
                elif len(absPath) < len(bestPath):
                    bestPath=absPath
        if bestPath is None:
            if returnNone:
                return None
            return path
        return bestPath
        
    def to_host_directory(self, path, source,dest):
        cleanDestination = os.path.normpath(dest)
        cleanPath= os.path.normpath(path)
        cleanSource=os.path.normpath(source)        
        #check if it is already relative to host path
        if cleanSource in cleanPath:
            return path

        # if the path is not mapping from host, will return path
        if cleanDestination not in cleanPath:
            return None
        abspath = os.path.normpath(str.join(os.sep,(cleanSource, path[path.find(cleanDestination) + len(cleanDestination):])))
        return abspath

class DockerThread_BuildImage(QThread):
    build_process = pyqtSignal(str)
    build_complete = pyqtSignal(int)

    def __init__(self, cli, name, path, dockerfile):
        QThread.__init__(self)
        self.docker = cli
        self.dockerfile = dockerfile
        self.name = name
        self.buildpath = path

    def __del__(self):
        self.wait()

    def run(self):
        print('Start building image {0} ---->>> \n docker file: {1} \n use path: {2}'.format(self.name, self.dockerfile, self.buildpath))
        try:
            # f = io.BytesIO(self.dockerfile.encode('utf-8'))
            for rawline in self.docker.getClient().build(path=self.buildpath, tag=self.name, dockerfile=self.dockerfile, rm=True):
                for jsonstr in rawline.decode('utf-8').split('\r\n')[:-1]:
                    log = jsonstr
                    try:
                        line = json.loads(jsonstr)
                        log = line['stream']
                    except ValueError as e:
                        print(e)
                    except TypeError as e:
                        print(e)
                    except KeyError as e:
                        log = ', '.join("{!s}={!r}".format(key, val) for (key, val) in line.items())
                    except:
                        log = ''
                    # print (log)
                    self.build_process.emit(log)

        except requests.exceptions.RequestException as e:
            self.build_process.emit(e.explanation)
            self.build_complete.emit(1)
        except Exception as e:
            self.build_process.emit(str(e))
            self.build_complete.emit(1)
            return


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
