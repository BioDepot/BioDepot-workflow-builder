import os, sys, re
import json
import requests
import subprocess
from docker import APIClient
from PyQt5.QtCore import QThread, pyqtSignal, QProcess, Qt
from PyQt5 import QtWidgets, QtGui, QtCore
from PyQt5.QtWidgets import *
from PyQt5.QtGui import *
import datetime
import uuid
import time
import shlex
from pathlib import Path
#runScheduler.sh ['b2b72d36270f4200a9e', '/tmp/docker.b2b72d36270f4200a9e.json', '8', '8096']
def breakpoint(title=None, message=None):
    warning_box = QMessageBox()
    warning_box.setIcon(QMessageBox.Warning)
    warning_box.setWindowTitle(title)
    warning_box.setText(message)
    warning_box.setStandardButtons(QMessageBox.Ok)
    warning_box.exec_()
    
def detect_container():
    # Check for Docker-specific markers
    if (os.path.isfile("/.dockerenv")):
        return True

    # Check for Apptainer/Singularity-specific markers
    if (
        os.environ.get("APPTAINER_NAME") or
        os.environ.get("SINGULARITY_NAME") or
        os.path.isdir("/.singularity.d") or
        os.path.isfile("/.apptainer")
    ):
        return True
    return False

class ConsoleProcess:
    #note that the logBaseDir is hard coded under /data/.bwb for now - this should be changed for whatever the primary volume mapping is
    # subclass that attaches a process and pipes the output to textedit widget console widget
    def __init__(self, console=None, errorHandler=None, finishHandler=None):
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
        #call cleanup just in case there is a previous processDir with same name
        self.cleanup()
        if console:
            self.log.readyReadStandardOutput.connect(
                lambda: self.writeConsole(
                    self.log, console, self.log.readAllStandardOutput, Qt.white
                )
            )
            self.log.readyReadStandardError.connect(
                lambda: self.writeConsole(
                    self.log, console, self.log.readAllStandardError, Qt.red
                )
            )
        if self.process:
            self.process.readyReadStandardOutput.connect(
                lambda: self.writeConsole(
                    self.process, self.console, self.process.readAllStandardOutput, Qt.white
                )
            )
            self.process.readyReadStandardError.connect(
                lambda: self.writeConsole(
                    self.process, self.console, self.process.readAllStandardError, Qt.red
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
        # Create a dict to store env vars for JSON return
        env_json = {}
        
        attrs = []
        groupSizes = []
        ramSizes = []
        maxThreads = 1
        sys.stderr.write("finished method init\n")
        nWorkers = 1
        
        if "nWorkers" in settings:
            nWorkers = settings["nWorkers"]
        env.insert("NWORKERS", "{}".format(nWorkers))
        env_json["NWORKERS"] = str(nWorkers)

        if "iteratedAttrs" not in settings or not settings["iteratedAttrs"]:
            return env_json
        
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
            # Add to JSON object
            env_json["ITERATEDATTRS"] = ":".join(attrs)
            env_json["GROUPSIZES"] = ":".join(groupSizes)
            env_json["RAMSIZES"] = ":".join(ramSizes)
        
        self.process.setProcessEnvironment(env)
        return env_json

    def addServerSettings(self, settings):
        pass

    def writeConsole(self, process, console, read, color):
        console.setTextColor(color)
        console.append(read().data().decode("utf-8", errors="ignore").rstrip())

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
            
    def startLog(self,schedule=False,namespace=None,global_env=None):
        if schedule:
            self.scheduleLog(namespace)
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
            
    def scheduleLog(self,namespace):
        self.log.start('getlog.sh {}'.format(namespace))
        
    def schedule(self,parms):
        job_id=parms[0]
        self.process.start('runScheduler.sh',parms)
        self.process.waitForFinished()
        self.scheduleLog(job_id)
        
    def start(self,cmds,outputFile=None,schedule=False,direct_json=None,global_env=None):
        self.cleanup()
        if schedule:
            sys.stderr.write('runScheduler.sh {}\n'.format(cmds))
            self.process.start('runScheduler.sh',cmds)
        else:
            os.makedirs(self.logDir,exist_ok=True)
            self.startLog()
            if direct_json:
                console_json = json.dumps({
                    "baseLogDir": self.baseLogDir,
                    "processDir": self.processDir,
                    "outputfile": outputFile,
                    "NWORKERS": global_env["NWORKERS"]
                })
#                breakpoint("start", "starting start with direct_json={}".format(direct_json))
#                breakpoint("start", "console_json={}".format(console_json))
                runCommand = "runDirectJob.sh"
                #make a dict out of direct_json and console_json
                total_json = {
                    "console_json": console_json,
                    "direct_json": direct_json
                }
                #combine console_json and direct_json into a single json string
                json_args = json.dumps(total_json)
                self.writeMessage("runDirectJob.sh {}".format(json_args))
                self.process.start(runCommand,[json_args])
            else:
                cmds.insert(0,self.baseLogDir)
                cmds.insert(0,self.processDir)
                cmds.insert(0,outputFile)
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
    def addThreadsRam(self, nThreads, ram):
        self.jsonObj["nThreads"] = nThreads
        self.jsonObj["ram"] = ram

    def addName(self, name):
        self.jsonObj["name"] = name

    def addDescription(self, description):
        self.jsonObj["description"] = description


class DockerJson:
    def __init__(self, tasksJson,namespace,name="",description=""):
        self.jsonObj = {}
        self.jsonObj["namespace"]=namespace
        self.jsonObj["name"]=name
        self.jsonObj["description"]=description
        self.jsonObj['tasks']=[]
        for taskJson in tasksJson:
            self.jsonObj['tasks'].append(taskJson.jsonObj)

class DockerClient:
    def __init__(self, url, name):
        self.url = url
        self.name = name
        #self.cli = APIClient(base_url=url)
        self.cli = None
        if (url):
            self.cli = APIClient(base_url=url)
        self.isContainer = detect_container()
        self.bwb_instance_id = None
        if(self.isContainer):   
            command = "awk -F'/containers/|/resolv.conf' '$2!=\"\" {print $2; exit}' /proc/self/mountinfo"
            outputString=str(subprocess.check_output(
                    command,
                    shell=True,
                    universal_newlines=True,
                ))
            if outputString:
                self.bwb_instance_id = outputString.splitlines()[0]
            else:
                outputString=str(subprocess.check_output(
                    "cat /proc/self/mountinfo | grep -oP '(?<=docker/containers/).*?(?=/resolv)'",
                    shell=True,
                 universal_newlines=True,
                ))
                self.bwb_instance_id = outputString.strip()
                
        print(self.bwb_instance_id)
        self.bwbMounts = {}
        self.shareMountPoint={};
        self.findVolumeMappings()
        self.findShareMountPoint(overwrite=True)        
        #self.findShareMountPoint()
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
        useGpu=False,
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
        for cmd in cmds:
            taskJson = TaskJson(name)
            taskJson.addBaseArgs(cmd)
            for env, var in environment.items():
                taskJson.addEnv(env, var)
            for container_dir, host_dir in hostVolumes.items():
                taskJson.addVolume(host_dir, container_dir, "rw")
            if exportGraphics:
                taskJson.addEnv("DISPLAY", "{}".format(os.getenv("DISPLAY")))
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

        namespace=str(uuid.uuid4().hex)[0:19]
        dockerJson = DockerJson(tasksJson,namespace=namespace)
        #jsonFile = "/data/dockerTest.json"
        jsonFile = "/tmp/docker.{}.json".format(namespace)
        with open(jsonFile, "w") as outfile:
            json.dump(dockerJson.jsonObj, outfile)
        parms=[namespace,jsonFile,cpuCount,memory]
        consoleProc.start(parms,schedule=True)
        
    def prettyEnv(self,var):
        if type(var) is list:
            output="["
            for v in var:
                #strip single quotes
                if v[0] == "'" and v[-1] =="'":
                    v=v[1:-1]
                output+='\\"{}\\",'.format(v)
            #delete extra comma and replace with ]
            output=output[:-1]+"]"
            #check if output is empty
            if output == "]":
                output="[]"
            return output       
        else:
            try:
                if var[0] == "'" and var[-1] =="'":
                    return var[1:-1]
            except TypeError:
                return var
            except IndexError:
                return var
            return var
    
    def build_jobs(self,commands, global_env=None):
        """
        Given a list of command strings, extract any environment variables specified
        using the "-e key=value" syntax and return a JSON object representing a list
        of jobs. Each job is a dict with two keys:
        - "command": the bare command with the -e parts removed.
        - "env": a dict of environment variable assignments.

        Additionally, if a global_env dictionary is provided, its key/value pairs
        will be merged into every job's environment dictionary. In case of conflicts,
        the command-specific environment variables take precedence.

        Example:
        Input:
            commands = ['echo "Hello" -e MY_VAR=some_value -e ANOTHER_VAR=another_value']
            global_env = {'GLOBAL_VAR': 'global_value', 'MY_VAR': 'global_override'}
        Output:
            [
                {
                "command": "echo Hello",
                "env": {
                    "GLOBAL_VAR": "global_value",
                    "MY_VAR": "some_value",         # command-specific value wins
                    "ANOTHER_VAR": "another_value"
                }
                }
            ]
        """
        jobs = []
        for cmd in commands:
            # Use shlex to properly split the command while respecting quotes.
            tokens = shlex.split(cmd)
            bare_tokens = []
            env = {}
            skip_next = False
            for i, token in enumerate(tokens):
                if skip_next:
                    # This token was consumed as the value after "-e"
                    skip_next = False
                    continue
                if token == '-e':
                    # Expect the next token to be in the form key=value.
                    if i + 1 < len(tokens):
                        env_assignment = tokens[i + 1]
                        if '=' in env_assignment:
                            key, value = env_assignment.split('=', 1)
                            env[key] = value
                        else:
                            print(f"Warning: Expected key=value after -e but got '{env_assignment}'")
                        skip_next = True  # Skip the next token as it is part of -e.
                    else:
                        print("Warning: '-e' flag found at the end of the command with no assignment.")
                else:
                    bare_tokens.append(token)
            # Rebuild the bare command as a string.
            bare_cmd = ' '.join(bare_tokens)
            
            # Merge the global environment variables (if provided) with the extracted ones.
            # Command-specific env vars override any conflicting keys from global_env.
            if global_env:
                #will always have NWORKERS in global_env - this ensures that the env in the json file is never empty
                merged_env = dict(global_env)  # Start with the global environment.
                merged_env.update(env)           # Override with the command-specific ones.
            else:
                merged_env = env

            job = {"command": bare_cmd, "env": merged_env}
            jobs.append(job)
        return json.dumps(jobs)   
    def direct_execute_iter(
        self,
        volumes=None,
        cmds=None,
        environment=None,
        hostVolumes=None,
        consoleProc=None,
        exportGraphics=False,
        useGpu=False,
        useContainer=True,
        portMappings=None,
        testMode=False,
        logFile=None,
        outputFile=None,
        scheduleSettings=None,
        iterateSettings=None,
        iterate=False,
        runDockerMap=False,
        nextFlowMap=False
    ):
    
    # reset logFile when it is not None - can be "" though - this allows an active reset
        if logFile is not None:
            self.logFile = logFile
        #build the global env dict
        global_env = {}
        for env, var in environment.items():
            # strip whitespace
            env.strip()
            # strip quotes if present
            if env[0] == env[-1] and env.startswith(("'", '"')):
                env = env[1:-1]
            global_env[env] = self.prettyEnv(var)
        if iterate and iterateSettings:
            sys.stderr.write("adding iterate settings\n")
            console_envs=consoleProc.addIterateSettings(iterateSettings)
            #add console_envs to global envs
            global_env.update(console_envs)
            sys.stderr.write("added iterate settings\n")
        else:
            global_env["NWORKERS"] = "1"
        #call build_jobs
        #breakpoint("build_jobs", "starting build_jobs with cmds={} and global_env={}".format(cmds, global_env))
        jobs_json = self.build_jobs(cmds, global_env)
        cmd_envs = [jobs_json]
        if testMode:
            baseCmd = ""
            echoStr = ""
            for cmd_env in cmd_envs:
                #generate a bash command string where we pass each of the cmd_envs to the command
                env_assignments = " ".join(["{}='{}'".format(key, value) for key, value in global_env.items()])
                job = json.loads(cmd_env)
                fullCmd = "{} {}".format(env_assignments, job["cmd"]).strip()
                echoStr = echoStr + fullCmd + "\n"
            print(echoStr)
            # Do not test for logFile - this may be None if it is not the first widget in testMode
            if self.logFile:
                with open(self.logFile, "a") as f:
                    f.write(echoStr)
            consoleProc.startTest(echoStr)
        else:
            sys.stderr.write("starting runDirectJob.sh\n")
            dockerCmds=[]
            consoleProc.start(dockerCmds,outputFile=outputFile,direct_json=cmd_envs,global_env=global_env)
    def create_container_iter(
        self,
        volumes=None,
        cmds=None,
        environment=None,
        hostVolumes=None,
        consoleProc=None,
        exportGraphics=False,
        useGpu=False,
        useContainer=True,
        portMappings=None,
        testMode=False,
        logFile=None,
        outputFile=None,
        scheduleSettings=None,
        iterateSettings=None,
        iterate=False,
        runDockerMap=False,
        nextFlowMap=False
    ):
        if not useContainer:
            #pass to the direct_execute_iter function
            return self.direct_execute_iter(
                volumes=volumes,
                cmds=cmds,
                environment=environment,
                hostVolumes=hostVolumes,
                consoleProc=consoleProc,
                exportGraphics=exportGraphics,
                useGpu=useGpu,
                portMappings=portMappings,
                testMode=testMode,
                logFile=logFile,
                outputFile=outputFile,
                scheduleSettings=scheduleSettings,
                iterateSettings=iterateSettings,
                iterate=iterate
            )
        # reset logFile when it is not None - can be "" though - this allows an active reset
        if logFile is not None:
            self.logFile = logFile
        volumeMappings = ""
        if nextFlowMap:
            volumeMappings = " -v /var/run/docker.sock:/var/run/docker.sock -v /tmp/.X11-unix:/tmp/.X11-unix {} ".format(self.findNextFlowSelfMounts())
        elif runDockerMap:
            volumeMappings = " -v /var/run/docker.sock:/var/run/docker.sock -v /tmp/.X11-unix:/tmp/.X11-unix --privileged "
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
            envs = envs + "-e {}={} ".format(env, self.prettyEnv(var))
        # create container cmds
        # the runDockerJob.sh script takes care of the first part of the docker command and cidfile
        # docker  run -i --rm --init --cidfile=<lockfile>'
        dockerBaseFlags = ""
        dockerCmds = []
        if useGpu:
            dockerBaseFlags += "--gpus all "
        if exportGraphics:
            dockerBaseFlags += "-e DISPLAY={} -v /tmp/.X11-unix:/tmp/.X11-unix ".format(os.getenv("DISPLAY"))
        if portMappings:
            for mapping in portMappings:
                dockerBaseFlags += "-p {} ".format(mapping)
        for cmd in cmds:
            dockerCmds.append(
                dockerBaseFlags + " {} {} {}".format(volumeMappings, envs, cmd)
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
            consoleProc.start(dockerCmds,outputFile=outputFile)
    
    def findShareMountPoint(self,overwrite=False):
        if not os.getenv('BWBSHARE' or overwrite):
            bwbshare=""
            bwbhostshare=""
            if self.bwbMounts:
                for key in self.bwbMounts:
                    if not bwbshare or "share" in self.bwbMounts[key]:
                        bwbhostshare=key+"/.bwbshare"
                        bwbshare=self.bwbMounts[key]+"/.bwbshare"
            if not bwbshare:
                bwbshare="/tmp/.X11/.bwbshare"
                bwbhostshare="/tmp/.X11/.bwbshare"
            #remove dir if present and make it
            os.system("rm -rf {}".format(bwbshare))
            os.system("mkdir -p {}".format(bwbhostshare))
            os.environ['BWBSHARE']=bwbshare
            os.environ['BWBHOSTSHARE']=bwbhostshare
            
        #check if mountpoint variable exists
    def findVolumeMappings(self):
        if not self.isContainer:
            self.bwbMounts["/data"] = "/data"
            return
        for c in self.cli.containers():
            container_id = c["Id"]
            if container_id == self.bwb_instance_id:
                for m in c["Mounts"]:
                    sys.stderr.write("Container mount points include {}\n".format(m))
                    if not (
                        "/var/run" in m["Source"] or "/tmp/.X11-unix" in m["Source"]
                    ):
                        self.bwbMounts[m["Source"]] = m["Destination"]
                        
    def findNextFlowSelfMounts(self):
        mountString=""
        for key in self.bwbMounts:
            mountString+=" -v {}:{} ".format(key,key)
        return mountString
        
    def to_best_host_directory(self, path, returnNone=False):
        if not detect_container():
            return path
        sys.stderr.write("bwbMounts are {}\n".format(self.bwbMounts))
        if self.bwbMounts == {}:
            self.findVolumeMappings()
            sys.stderr.write(
                "bwbMounts after findVolume are {}\n".format(self.bwbMounts)
            )
            if not hasattr(self,'shareMountPoint') or not self.shareMountPoint:
                self.shareMountPoint["bwb"]="/tmp/.X11/.bwb"
                self.shareMountPoint["host"]="/tmp/.X11/.bwb"
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
        # Must check for equality or the abspath may have multiple leading slashes which are not cleaned up by normpath
        if cleanPath == cleanDestination:
            return cleanSource
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
