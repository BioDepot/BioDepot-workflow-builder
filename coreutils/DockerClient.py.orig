import os, sys, re
import json
import requests
import subprocess
from docker import APIClient
from PyQt5.QtCore import QThread, pyqtSignal, QProcess, Qt
from PyQt5 import QtWidgets, QtGui, QtCore
import datetime
import uuid
import time
from pathlib import Path
#runScheduler.sh ['b2b72d36270f4200a9e', '/tmp/docker.b2b72d36270f4200a9e.json', '8', '8096']
class ConsoleProcess:
    #note that the logBaseDir is hard coded under /data/.bwb for now - this should be changed for whatever the primary volume mapping is
    # subclass that attaches a process and pipes the output to textedit widget console widget
    def __init__(self, console=None, errorHandler=None, finishHandler=None, logHandler=None):
        self.threadNumber=0
        self.processDir="proc.{}".format(uuid.uuid4().hex)
        self.startupOnly=False
        self.finishHandler=None
        self.log = QProcess()
        self.baseLogDir="/data/.bwb"
        self.logDir="{}/{}/logs".format(self.baseLogDir,self.processDir)
        self.logFile="{}/log{}".format(self.logDir,self.threadNumber)
        self.errorDir="/tmp/{}/errors".format(self.processDir)
        self.process = QProcess()
        self.console = console
        self.state = "stopped"
<<<<<<< HEAD
        self.stdoutColor=Qt.white
        self.stderrColor=Qt.red
        self.job_id=None
        self.schedule_pid=0
=======
        #call cleanup just in case there is a previous processDir with same name
        self.cleanup()
>>>>>>> f0a94431916b263daed02d8b06468b720099d07b
        if console:
            self.log.readyReadStandardOutput.connect(
                lambda: self.writeConsole(
                    self.log, console, self.log.readAllStandardOutput, self.stdoutColor
                )
            )
            self.log.readyReadStandardError.connect(
                lambda: self.writeConsole(
                    self.log, console, self.log.readAllStandardError, self.stderrColor
                )
            )        
        if finishHandler:
            self.finishHandler=finishHandler
        self.process.finished.connect(self.onFinish)

        
    def changeThreadNumber(self,newThreadNumber):
        if newThreadNumber == self.threadNumber:
            return
        self.threadNumber = newThreadNumber
        self.logFile="{}/log{}".format(self.logDir,self.threadNumber)
        self.changeConsole()
    
    def changeConsole(self):
        if self.process.state() == 0:
            self.writeFileToConsole(self.logFile)
        else:
            sys.stderr.write("terminating log\n")
            self.log.kill()
            self.log.waitForFinished(10)
            if self.process.state():
                sys.stderr.write("starting new log\n")
                self.startLog()
            else:
                self.writeFileToConsole(self.logFile)
                
    def writeFileToConsole(self,filename):
        if filename and Path(filename).is_file():
            self.log.kill()
            self.log.waitForFinished(10)            
            if (self.log.state() == 0):
                self.log.start('cat',[filename])
        
    def addIterateSettings(self, settings):
        env = QtCore.QProcessEnvironment.systemEnvironment()
        attrs = []
        groupSizes = []
        ramSizes = []
        maxThreads = 1
        sys.stderr.write("finished method init\n")
        nWorkers = 1
        if "nWorkers" in settings:
            nWorkers = settings["nWorkers"]
        env.insert("NWORKERS", "{}".format(nWorkers))

        if "iteratedAttrs" not in settings or not settings["iteratedAttrs"]:
            return
        sys.stderr.write("checking attrs in settings\n")
        for attr in settings["iteratedAttrs"]:
            sys.stderr.write("adding attr {}\n".format(attr))
            attrs.append(attr)
            if (
                attr in settings["data"]
                and "threads" in settings["data"][attr]
                and settings["data"][attr]["threads"]
            ):
                if int(settings["data"][attr]["threads"]) > maxThreads:
                    maxThreads = int(settings["data"][attr]["threads"])
            if (
                attr in settings["data"]
                and "groupSize" in settings["data"][attr]
                and settings["data"][attr]["groupSize"]
            ):
                groupSizes.append(settings["data"][attr]["groupSize"])
            else:
                groupSizes.append("1")
            if (
                attr in settings["data"]
                and "ram" in settings["data"][attr]
                and settings["data"][attr]["ram"]
            ):
                ramSizes.append(settings["data"][attr]["ram"])
            else:
                ramSizes.append("0")
        # need to add code to pass environment arrays to process
        if attrs:
            env.insert("ITERATEDATTRS", ":".join(attrs))
            env.insert("GROUPSIZES", ":".join(groupSizes))
            env.insert("RAMSIZES", ":".join(ramSizes))
        self.process.setProcessEnvironment(env)

    def addServerSettings(self, settings):
        pass

    def writeConsole(self, process, console, read, color=None):
        console.setTextColor(color)
        console.append(read().data().decode('utf-8',errors="ignore"))

    def writeMessage(self, message, color=Qt.green):
        # for bwb messages
        self.console.setTextColor(color)
        self.console.append(message)

    def stop(self, message=None):
        self.state = "stopped"
        # the runDockerJob.sh cleans itself up when interrupted
        self.process.terminate()
        self.log.terminate()
        if message:
            self.writeMessage(message)
            
    def startLog(self,schedule=False,job_id=None):
        if schedule:
            self.scheduleLog(job_id)
            return
        self.logFile='{}/log{}'.format(self.logDir,self.threadNumber)
        myPath=Path(self.logFile)
        if not myPath.is_file():
            myPath.touch()
        self.writeMessage('tail -c +0 -f {}'.format(self.logFile))
        tailParams=['-c','+0','-f',self.logFile]
        self.log.start('tail',tailParams)
    
    def startTest(self,cmdString):
        self.writeMessage(cmdString)
        if self.finishHandler:
            self.finishHandler()
            
    def startScheduleLog(self):
        
<<<<<<< HEAD
        self.log.start('getlog.sh',[self.job_id])
        
    def start(self,cmds,schedule=False,job_id=None):
=======
    def schedule(self,parms):
        job_id=parms[0]
        self.process.start('runScheduler.sh',parms)
        self.process.waitForFinished()
        self.scheduleLog(job_id)
        
    def start(self,cmds,schedule=False):
>>>>>>> f0a94431916b263daed02d8b06468b720099d07b
        self.cleanup()
        if not job_id:
            self.job_id=str(uuid.uuid4().hex)[0:19]
        self.job_id=job_id
        self.logDir='logs.{}'.format(self.job_id)
        os.makedirs('/data/.bwb/{}'.format(self.logDir))
        if schedule:
            self.startScheduleLog()
            sys.stderr.write('runScheduler.sh {}\n'.format(cmds))
            self.process.start('runScheduler.sh',cmds)
        else:
<<<<<<< HEAD
=======
            os.makedirs(self.logDir)
>>>>>>> f0a94431916b263daed02d8b06468b720099d07b
            self.startLog()
            cmds.insert(0,self.baseLogDir)
            cmds.insert(0,self.processDir)
            self.writeMessage('runDockerJob.sh {}'.format(cmds))
            self.process.start('runDockerJob.sh',cmds)
        
    def onFinish(self,code,status):
        if self.startupOnly:
            return
        if self.finishHandler:
            if not status:
                status=self.checkForErrors()
            self.finishHandler(code=code,status=status)
        self.log.terminate()

    def checkForErrors(self):
        if os.path.exists(self.errorDir):
            return -1;
        return 0
        
    def cleanup(self):
        if self.baseLogDir and self.processDir:
            path="{}/".format(self.baseLogDir,self.processDir)
            if os.path.exists(path):
                os.system("rm {} -r".format(path))
        if self.processDir:
            path="/tmp/{}".format(self.processDir)
            if os.path.exists(path):
                os.system("rm {} -r".format(path))            


class TaskJson:
    #nWorkers same as maxWorkers - both kept to not break older pipelines
    def __init__(self, imageName):
        self.jsonObj={}
        self.jsonObj["image"] = imageName
        self.jsonObj["args"] = []
        self.jsonObj["envs"] = []
        self.jsonObj["volumes"] = []

    def addBaseArgs(self, cmd):
        #self.jsonObj["command"]["args"] = ["-i", "--rm", "--init", cmd]
        self.jsonObj["args"] = [cmd]

    def addVolume(self, host_dir, container_dir, mode):
        volumeMapping = {}
        volumeMapping["host_dir"] = host_dir
        volumeMapping["mount_dir"] = container_dir
        volumeMapping["mode"] = mode
        self.jsonObj["volumes"].append(volumeMapping)

    def addEnv(self, key, val):
        key.strip()
        # strip quotes if present
        if key[0] == key[-1] and key.startswith(("'", '"')):
            key = key[1:-1]
        self.jsonObj["env"].append({"key": key, "val": val})

    # need to fix this so that we can add parameters for maxWorkers and threads per worker
    def addMaxWorkers(self, maxWorkers):
        self.jsonObj['maxWorkers']=maxWorkers
        
    def addThreadsRam(self, nThreads, ram):
        self.jsonObj["nThreads"] = nThreads
        self.jsonObj["ram"] = ram

    def addName(self, name):
        self.jsonObj["name"] = name

    def addDescription(self, description):
        self.jsonObj["description"] = description


class DockerJson:
    def __init__(self, tasksJson,job_id,name="",description="",maxWorkers=1):
        self.jsonObj = {}
        self.jsonObj["job_id"]=job_id
        self.jsonObj["name"]=name
        self.jsonObj["description"]=description
        self.jsonObj["maxWorkers"]=maxWorkers
        self.jsonObj['tasks']=[]
        for taskJson in tasksJson:
            self.jsonObj['tasks'].append(taskJson.jsonObj)

class DockerClient:
    def __init__(self, url, name):
        self.url = url
        self.name = name
        self.cli = APIClient(base_url=url)
        self.bwb_instance_id = str(
            subprocess.check_output(
                'cat /proc/self/cgroup | head -1 | cut -d "/" -f3',
                shell=True,
                universal_newlines=True,
            )
        ).splitlines()[0]
        self.bwbMounts = {}
        self.findVolumeMappings()
        self.logFile = None
        self.schedulerStarted = False

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
        repoTag = name + ":" + version
        conId = subprocess.check_output(["docker", "images", "-q", repoTag])
        if conId:
            return True
        return False

    def remove_image(self, id, force=False):
        self.cli.remove_image(id, force=force)

    def pull_image(self, id):
        self.cli.pull(id)

    def containers(self, all=True):
        return self.cli.containers(all=all)

    def findMaxIterateValues(self, settings):
        maxThreads = 0
        maxRam = 0
        if "iteratedAttrs" in settings:
            for attr in settings["iteratedAttrs"]:
                if (
                    attr in settings["data"]
                    and "threads" in settings["data"][attr]
                    and settings["data"][attr]["threads"]
                ):
                    if int(settings["data"][attr]["threads"]) > maxThreads:
                        maxThreads = int(settings["data"][attr]["threads"])
                if (
                    attr in settings["data"]
                    and "ram" in settings["data"][attr]
                    and settings["data"][attr]["ram"]
                ):
                    ramSize = int(settings["data"][attr]["ram"])
                    if ramSize > maxRam:
                        maxRam = ramSize
        return maxThreads, maxRam

    def create_container_external(
        self,
        name,
        volumes=None,
        cmds=None,
        environment=None,
        hostVolumes=None,
        consoleProc=None,
        exportGraphics=False,
        portMappings=None,
        testMode=False,
        logFile=None,
        scheduleSettings=None,
        iterateSettings=None,
        iterate=False,
    ):
        tasksJson = []
        count = 0
        cpuCount='8'
        memory='8096'
        maxWorkers=1
        if iterateSettings and 'nWorkers' in iterateSettings:
            maxWorkers=iterateSettings['nWorkers']
        for cmd in cmds:
            taskJson = TaskJson(name)
            taskJson.addBaseArgs(cmd)
            for env, var in environment.items():
                taskJson.addEnv(env, var)
            for container_dir, host_dir in hostVolumes.items():
                taskJson.addVolume(host_dir, container_dir, "rw")
            if exportGraphics:
                taskJson.addEnv("DISPLAY", ":1")
                taskJson.addVolume("/tmp/.X11-unix", "/tmp/.X11-unix", "rw")
            maxThreads = 1
            maxRam = 0
            if iterate and iterateSettings:
                maxThreads, maxRam = self.findMaxIterateValues(iterateSettings)
            taskJson.addThreadsRam(maxThreads, maxRam)
            taskJson.addName("cmdName{}".format(count))
            taskJson.addDescription("command{}".format(count))
            tasksJson.append(taskJson)
            count += 1

        job_id=str(uuid.uuid4().hex)[0:19]
        dockerJson = DockerJson(tasksJson,job_id=job_id,maxWorkers=maxWorkers)
        #jsonFile = "/data/dockerTest.json"
        jsonFile = "/tmp/docker.{}.json".format(job_id)
        with open(jsonFile, "w") as outfile:
            json.dump(dockerJson.jsonObj, outfile)
<<<<<<< HEAD
        parms=[job_id,jsonFile,cpuCount,memory]
        consoleProc.start(parms,schedule=True,job_id=job_id)

=======
        parms=[namespace,jsonFile,cpuCount,memory]
        consoleProc.start(parms,schedule=True)
    def arrayPrint(self,var):
        #remove spaces from array 
        if type(var) is list or type(var) is tuple:
            return"[{}]".format(",".join(map(repr, var)))
        else:
            return var
>>>>>>> f0a94431916b263daed02d8b06468b720099d07b
    def create_container_iter(
        self,
        name,
        volumes=None,
        cmds=None,
        environment=None,
        hostVolumes=None,
        consoleProc=None,
        exportGraphics=False,
        portMappings=None,
        testMode=False,
        logFile=None,
        scheduleSettings=None,
        iterateSettings=None,
        iterate=False,
    ):
        # reset logFile when it is not None - can be "" though - this allows an active reset
        if logFile is not None:
            self.logFile = logFile
        volumeMappings = ""
        for container_dir, host_dir in hostVolumes.items():
            volumeMappings = volumeMappings + "-v {}:{} ".format(
                self.to_best_host_directory(host_dir), container_dir
            )
        envs = ""
        for env, var in environment.items():
            # strip whitespace
            env.strip()
            # strip quotes if present
            if env[0] == env[-1] and env.startswith(("'", '"')):
                env = env[1:-1]
            envs = envs + "-e {}={} ".format(env, self.arrayPrint(var))
        # create container cmds
        # the runDockerJob.sh script takes care of the first part of the docker command and cidfile
        # docker  run -i --rm --init --cidfile=<lockfile>'
        dockerBaseFlags = ""
        dockerCmds = []
        if exportGraphics:
            dockerBaseFlags += "-e DISPLAY=:1 -v /tmp/.X11-unix:/tmp/.X11-unix "
        if portMappings:
            for mapping in portMappings:
                dockerBaseFlags += "-p {} ".format(mapping)
        for cmd in cmds:
            dockerCmds.append(
                dockerBaseFlags + " {} {} {} {}".format(volumeMappings, envs, name, cmd)
            )
        consoleProc.state = "running"
        for dcmd in dockerCmds:
            sys.stderr.write("{}\n".format(dcmd))
        # pass on iterateSettings
        if iterate and iterateSettings:
            sys.stderr.write("adding iterate settings\n")
            consoleProc.addIterateSettings(iterateSettings)
            sys.stderr.write("added iterate settings\n")
        else:
            envs + "-e NWORKERS=1 "
        if testMode:
            baseCmd = "docker  run -i --rm --init "
            echoStr = ""
            for dockerCmd in dockerCmds:
                fullCmd = baseCmd + dockerCmd
                echoStr = echoStr + fullCmd + "\n"
            print(echoStr)
            # Do not test for logFile - this may be None if it is not the first widget in testMode
            if self.logFile:
                with open(self.logFile, "a") as f:
                    f.write(echoStr)

            consoleProc.startTest(echoStr)
        else:
            sys.stderr.write("starting runDockerJob.sh\n")
            consoleProc.start(dockerCmds)

    def findVolumeMappings(self):
        for c in self.cli.containers():
            container_id = c["Id"]
            if container_id == self.bwb_instance_id:
                for m in c["Mounts"]:
                    sys.stderr.write("Container mount points include {}\n".format(m))
                    if not (
                        "/var/run" in m["Source"] or "/tmp/.X11-unix" in m["Source"]
                    ):
                        self.bwbMounts[m["Source"]] = m["Destination"]

    def to_best_host_directory(self, path, returnNone=False):
        sys.stderr.write("bwbMounts are {}\n".format(self.bwbMounts))
        if self.bwbMounts == {}:
            self.findVolumeMappings()
            sys.stderr.write(
                "bwbMounts after findVolume are {}\n".format(self.bwbMounts)
            )
        bestPath = None
        for source, dest in self.bwbMounts.items():
            absPath = self.to_host_directory(path, source, dest)
            if absPath is not None:
                if bestPath is None:
                    bestPath = absPath
                elif len(absPath) < len(bestPath):
                    bestPath = absPath
        if bestPath is None:
            if returnNone:
                return None
            return path
        return bestPath

    def to_host_directory(self, path, source, dest):
        cleanDestination = os.path.normpath(dest)
        cleanPath = os.path.normpath(path)
        cleanSource = os.path.normpath(source)
        # check if it is already relative to host path
        if cleanSource in cleanPath:
            return path

        # if the path is not mapping from host, will return path
        if cleanDestination not in cleanPath:
            return None
        abspath = os.path.normpath(
            str.join(
                os.sep,
                (
                    cleanSource,
                    path[path.find(cleanDestination) + len(cleanDestination) :],
                ),
            )
        )
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
        print(
            "Start building image {0} ---->>> \n docker file: {1} \n use path: {2}".format(
                self.name, self.dockerfile, self.buildpath
            )
        )
        try:
            # f = io.BytesIO(self.dockerfile.encode('utf-8'))
            for rawline in self.docker.getClient().build(
                path=self.buildpath, tag=self.name, dockerfile=self.dockerfile, rm=True
            ):
                for jsonstr in rawline.decode("utf-8").split("\r\n")[:-1]:
                    log = jsonstr
                    try:
                        line = json.loads(jsonstr)
                        log = line["stream"]
                    except ValueError as e:
                        print(e)
                    except TypeError as e:
                        print(e)
                    except KeyError as e:
                        log = ", ".join(
                            "{!s}={!r}".format(key, val) for (key, val) in line.items()
                        )
                    except:
                        log = ""
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
        repo_tag = self.name + ":" + self.version
        try:
            # Docker splits downloads into multiple parts
            # We create a dict mapping id to % finished for each part (progress)
            # The total progress is the mean of individual progresses
            progs = dict()
            for line in self.docker.getClient().pull(repo_tag, stream=True):
                for status in line.decode("utf-8").split("\r\n")[:-1]:
                    line = json.loads(status)
                    statusStr = line["status"]
                    if statusStr == "Pulling fs layer":
                        progs[line["id"]] = 0
                    # First 50% progress is Downloading
                    elif statusStr == "Downloading":
                        progDetail = line["progressDetail"]
                        if len(progDetail) > 1:
                            progs[line["id"]] = (
                                progDetail["current"] / progDetail["total"] * 50
                            )
                    # Last 50% progress is Extracting
                    elif statusStr == "Extracting":
                        progDetail = line["progressDetail"]
                        if len(progDetail) > 1:
                            progs[line["id"]] = 50 + (
                                progDetail["current"] / progDetail["total"] * 50
                            )
                    if len(progs) > 0:
                        self.current_progress = sum(progs.values()) / len(progs)
                        self.pull_progress.emit(self.current_progress)
        except requests.exceptions.RequestException:
            # TODO emit error
            print("Connection Exception!")
