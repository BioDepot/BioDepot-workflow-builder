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
        self.state='stopped'
        if console:
            self.process.readyReadStandardOutput.connect(lambda: self.writeConsole(self.process, console,self.process.readAllStandardOutput,Qt.white))
            self.process.readyReadStandardError.connect(lambda: self.writeConsole(self.process, console,self.process.readAllStandardError,Qt.red))
        if finishHandler:
            self.process.finished.connect(finishHandler)
    def addIterateSettings(self,settings):
        env=QtCore.QProcessEnvironment.systemEnvironment()
        attrs=[]
        groupSizes=[]
        ramSizes=[]
        maxThreads=1
        env.insert("NWORKERS", "{}".format(settings['nWorkers']))
        if 'iteratedAttrs' not in settings or not settings['iteratedAttrs']:
            return
        for attr in settings['iteratedAttrs']:
            attrs.append(attr)
            if attr in settings['data'] and 'threads' in settings['data'][attr] and settings['data'][attr]['threads']:
                if int(settings['data'][attr]['threads']) > maxThreads:
                    maxThreads=int(settings['data'][attr]['threads'])
            if attr in settings['data'] and 'groupSize' in settings['data'][attr] and settings['data'][attr]['groupSize']:
                groupSizes.append(settings['data'][attr]['groupSize'])
            else:
                groupSizes.append('1')
            if attr in settings['data'] and 'ram' in settings['data'][attr] and settings['data'][attr]['ram']:
                ramSizes.append(settings['data'][attr]['ram'])
            else:
                ramSizes.append('0')            
        #need to add code to pass environment arrays to process
        if attrs:
            env.insert("ITERATEDATTRS",":".join(attrs))
            env.insert("GROUPSIZES",":".join(groupSizes))
            env.insert('RAMSIZES',":".join(ramSizes))
        self.process.setProcessEnvironment(env)

    def addServerSettings(self,settings):
        pass
    
    def writeConsole(self,process,console,read,color):
        console.setTextColor(color)
        console.append(read().data().decode('utf-8',errors="ignore"))
    
    def writeMessage(self,message,color=Qt.green):
        #for bwb messages
        self.console.setTextColor(color)
        self.console.append(message)
        
    def stop(self,message=None):
        self.state='stopped'
        #the runDockerJob.sh cleans itself up when interrupted
        self.process.terminate()
        if message:
            self.writeMessage(message)

        
class CmdJson:
    def __init__(self):
        self.jsonObj={}
        self.jsonObj['args']=[]
        self.jsonObj['image']=name
        self.jsonObj['deps']=""
        self.jsonObj['envs']=[]
        self.jsonObj['volumes']=[]

    def addBaseArgs(self,cmd):
        self.jsonObj['args']=['-i', '--rm', '--init',cmd]

    def addVolume(self,host_dir,mount_dir,mode):
        volumeMapping['host_dir']=host_dir
        volumeMapping['mount_dir']=container_dir
        volumeMapping['mode']=mode
        self.jsonObj['volumes'].append(volumeMapping)
        
    def addEnv(self,key,val):
        key.strip()
        #strip quotes if present
        if key[0] == key[-1] and key.startswith(("'",'"')):
            key=key[1:-1]
        envDict['key']=key
        envDict['val']=val
        self.jsonObj['envs'].append(envDict)
    
    #need to fix this so that we can add parameters for maxWorkers and threads per worker    
    def addThreadsRam(self,nThreads,ram):
        self.jsonObj['nThreads']=nThreads
        self.jsonObj['ram']=ram

        
class TaskJson:
    def __init__(self,name,description,cmds):
        self.jsonObj['commands']={}
        #cmds is an list of cmdJson object
        self.jsonObj['commands']['command']=cmds 
        self.jsonObj['commands']['name']=name
        self.jsonObj['commands']['description']=description
        
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
    
    def findMaxIterateValues(settings):
        maxThreads=0
        maxRam=0
        for attr in settings['iteratedAttrs']:
            if attr in settings['data'] and 'threads' in settings['data'][attr] and settings['data'][attr]['threads']:
                if int(settings['data'][attr]['threads']) > maxThreads:
                    maxThreads=int(settings['data'][attr]['threads'])
            if attr in settings['data'] and 'ram' in settings['data'][attr] and settings['data'][attr]['ram']:
                ramSize=int(settings['data'][attr]['ram'])
                if ramSizes > maxRam:
                    maxRam =ramSizes
        return maxThreads, maxRam
        
    def create_container_external(self, name, volumes=None, cmds=None, environment=None, hostVolumes=None,  exportGraphics=False, portMappings=None,testMode=False,logFile=None,scheduleSettings=None,iterateSettings=None):
        cmdsJson=[]
        for cmd in cmds:
            cmdJson=CmdJson()
            cmdJson.addBaseArgs(cmd)
            for env, var in environment.items():
                cmdJson.addEnv(env,var)
            for container_dir, host_dir in hostVolumes.items():
                cmdJson.addVolume(host_dir,container_dir,'rw')
            if exportGraphics:
                cmdJson.addEnv('DISPLAY',':1')
                cmdJson.addVolume('/tmp/.X11-unix','/tmp/.X11-unix','rw')
            maxThreads=1
            maxRam=0
            if iterateSettings:
                maxThreads,maxRam=findMaxIterateValues(iterateSettings)
            cmdJson.addThreadsRam(maxThreads,maxRam)
        taskJson=TaskJson('testTask','test description',cmds)

    def create_container_iter(self, name, volumes=None, cmds=None, environment=None, hostVolumes=None, consoleProc=None, exportGraphics=False, portMappings=None,testMode=False,logFile=None,scheduleSettings=None,iterateSettings=None):
        #reset logFile when it is not None - can be "" though - this allows an active reset
        if logFile is not None:
            self.logFile = logFile
         
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
        #create container cmds
        #the runDockerJob.sh script takes care of the first part of the docker command and cidfile
        #docker  run -i --rm --init --cidfile=<lockfile>'
        dockerBaseFlags=""
        dockerCmds=[]
        if exportGraphics:
            dockerBaseFlags+='-e DISPLAY=:1 -v /tmp/.X11-unix:/tmp/.X11-unix '
        if portMappings:
            for mapping in portMappings:
                dockerBaseFlags+='-p {} '.format(mapping)
        for cmd in cmds:
            dockerCmds.append(dockerBaseFlags + ' {} {} {} {}'.format(volumeMappings,envs,name,cmd))
        consoleProc.state='running'
        
        #pass on iterateSettings
        if iterateSettings:
            consoleProc.addIterateSettings(iterateSettings)
        else:
            #need to have NWORKERS set
             env.insert("NWORKERS", "1")
        if testMode:
            baseCmd='docker  run -i --rm --init '
            echoStr=''
            for dockerCmd in dockerCmds:
                fullCmd=baseCmd+dockerCmd
                echoStr=echoStr+fullCmd+'\n'
            print (echoStr)
            #Do not test for logFile - this may be None if it is not the first widget in testMode
            if self.logFile:
                with open (self.logFile,'a') as f:
                    f.write(echoStr)
                   
            consoleProc.process.start('echo',[echoStr])
        else:
            consoleProc.process.start('runDockerJob.sh',dockerCmds)
                   
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
