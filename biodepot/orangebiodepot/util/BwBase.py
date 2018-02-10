import os
from AnyQt.QtCore import QThread, pyqtSignal, Qt
from Orange.widgets import widget, gui, settings
from orangebiodepot.util.DockerClient import DockerClient, PullImageThread
from PyQt5 import QtWidgets, QtGui

class LocalContainerRunner(QThread):
    progress = pyqtSignal(int)

    def __init__(self, cli, image_name, volumes, commands, environments = None):
        QThread.__init__(self)
        self.docker = cli
        self.image_name = image_name
        self.volumes = volumes
        self.commands = commands
        self.environments = environments
        self.containerId = ""
        self.Flag_isRunning = False

    def __del__(self):
        self.wait()

    def run(self):
        response = self.docker.create_container(self.image_name, hostVolumes=self.volumes, commands=self.commands, environment=self.environments)
        if response['Warnings'] == None:
            self.containerId = response['Id']
            self.docker.start_container(self.containerId)
        else:
            print(response['Warnings'])

        # Keep running until container is exited
        while self.docker.container_running(self.containerId):
            self.sleep(1)
        # Remove the container when it is finished
        self.docker.remove_container(self.containerId)

class ConnectionDict:
    def __init__(self, inputDict):
        self._dict=inputDict #we do want the name not a copy - i.e. the inputDict should change
        
    def add (self, slot, connectionId=None):
        if slot in self._dict:
            self._dict[slot].append(connectionId)
        else:
            self._dict[slot]=[connectionId]
            
    def remove (self, slot, connectionId=None):
        if slot in self._dict:
            try:
                self._dict[slot].remove(connectionId)
            except ValueError:
                pass            
    def isConnected (self, slot):
        if self._dict is None:
            return False
        if (slot in self._dict) and (self._dict[slot]):
            return True
        return False
        
class OWBwBWidget(widget.OWWidget):

    dockerClient = DockerClient('unix:///var/run/docker.sock', 'local')
    defaultFileIcon=QtGui.QIcon('/biodepot/orangebiodepot/icons/file.png')

    def __init__(self, image_name, image_tag):
        super().__init__()
        self._hostDirectories = {}
        self._dockerImageName = image_name
        self._dockerImageTag = image_tag
        self._dockerVolumes = None
        self._dockerCommands = None
        self._dockerEnvironments = None
        self._Flag_isRunning = False
    """
    GUI elements
    """
    def bwbFileEntry(self, widget, button, ledit, icon=defaultFileIcon,label=None,entryType='file'):
        button.setStyleSheet("border: 1px solid #1a8ac6; border-radius: 2px;")
        button.setIcon(icon)
        ledit.setClearButtonEnabled(True)
        ledit.setPlaceholderText("Enter {}".format(entryType))    
        myLayout=QtGui.QHBoxLayout()
        if label:
            myLabel=QtGui.QLabel(label)
            myLabel.setAlignment(Qt.AlignRight | Qt.AlignVCenter) 
            myLayout.addWidget(myLabel)
        myLayout.addWidget(ledit)
        myLayout.addWidget(button)
        widget.layout().addLayout(myLayout)
    
    def Event_OnRunFinished(self):
        raise Exception('Event_OnRunFinished not implemented!')

    def Event_OnRunMessage(self, message):
        raise Exception('Event_OnRunMessage not implemented!')

    def Initialize_InputKeys(self, keys):
        for k in keys:
            self._hostDirectories[k] = None

    def setDirectories(self, key, path, ctrlLabel = None):
        if ctrlLabel:
            labelTitle = ctrlLabel.text().split(':')[0] + ': {}'
        # When a user removes a connected Directory widget,
        # it sends a signal with path=None
        if path is None:
            self._hostDirectories[key] = None
            if ctrlLabel is not None: ctrlLabel.setText(labelTitle.format('Not Set'))
        else:
            self._hostDirectories[key] = path.strip()
            if ctrlLabel is not None: ctrlLabel.setText(labelTitle.format((str(path))))


    def getDirectory(self, key):
        if key in self._hostDirectories:
            return self._hostDirectories[key]
        return None

    """
    Pull image
    """
    def __dockerPullImage__(self):
        self.Event_OnRunMessage('Pulling \'' + self._dockerImageName + ":" + self._dockerImageTag+ '\' from Dockerhub...')
        self.setStatusMessage("Downloading...")
        self.progressBarInit()
        self._Flag_isRunning = True
        # Pull the image in a new thread
        self.pullImageThread = PullImageThread(self.dockerClient, self._dockerImageName, self._dockerImageTag)
        self.pullImageThread.pull_progress.connect(self.__dockerPullProgress__)
        self.pullImageThread.finished.connect(self.__dockerPullDone__)
        self.pullImageThread.start()

    def __dockerPullProgress__(self, val):
        self.progressBarSet(val)

    def __dockerPullDone__(self):
        self.Event_OnRunMessage('Finished pulling \'' + self._dockerImageName + ":" + self._dockerImageTag + '\'')
        self.progressBarFinished()
        self.__dockerRealRun__()

    """
    Run container
    """
    def dockerRun(self, volumes = None, commands = None, environments = None):
        if not self._Flag_isRunning:
            self._dockerVolumes = volumes
            self._dockerCommands = commands
            self._dockerEnvironments = environments
            # Make sure the docker image is downloaded
            if not self.dockerClient.has_image(self._dockerImageName, self._dockerImageTag):
                self.__dockerPullImage__()
            else:
                self.__dockerRealRun__()

    def __dockerRealRun__(self):
        self._Flag_isRunning = True
        self.setStatusMessage('Running...')
        self.Event_OnRunMessage('Running \'' + self._dockerImageName + '\'')
        # Run the container in a new thread
        self.containerThread = LocalContainerRunner(
                                self.dockerClient,
                                '{}:{}'.format(self._dockerImageName, self._dockerImageTag),
                                self._dockerVolumes,
                                self._dockerCommands,
                                self._dockerEnvironments)
        self.containerThread.progress.connect(self.__dockerRunProgress__)
        self.containerThread.finished.connect(self.__dockerRunDone__)
        self.containerThread.start()

    def __dockerRunProgress__(self, val):
        self.progressBarSet(val)

    def __dockerRunDone__(self):
        self._Flag_isRunning = False
        self.setStatusMessage('Finished!')
        self.Event_OnRunMessage('Finished!')
        self.Event_OnRunFinished()
        
    def generateCmd (self, executables = [], flags = {}, args = []):
        #flags are key value pairs - values begin with = if they are to be joined with the key i.e. --name=John and not --name John 
        cmdStr=''
        for executable in executables:
            cmdStr += executable + ' '
        #short flags no args first
        #short flags args
        #long flags no args
        #long flags args
        #flags no dashs
        #flags no args
        for flagName, flagValue in flags.items():
            if (flagName[:1] == '-') and (flagName[:2] != '--') and (not flagValue):
                cmdStr +=  flagName + ' '
        for flagName,flagValue in flags.items():
            if (flagName[:1] == '-') and (flagName[:2] != '--') and flagValue:
                if(flagValue[:1] == '='):
                    cmdStr +=  flagName + flagValue + ' '
                else:
                    cmdStr +=  flagName + ' ' + flagValue + ' '
        for flagName, flagValue in flags.items():
            if (flagName[:2] == '--')  and (not flagValue):
                cmdStr +=  flagName + ' '
        for flagName,flagValue in flags.items():
            if (flagName[:2] == '--')  and flagValue:
                if(flagValue[:1] == '='):
                    cmdStr +=  flagName + flagValue + ' '
                else:
                    cmdStr +=  flagName + ' ' + flagValue + ' '
        for flagName, flagValue in flags.items():
            if (flagName[:1] != '-')  and (not flagValue):
                cmdStr +=  flagName + ' '
        for flagName,flagValue in flags.items():
            if (flagName[:1] != '-')  and flagValue:
                if(flagValue[:1] == '='):
                    cmdStr +=  flagName + flagValue + ' '
                else:
                    cmdStr +=  flagName + ' ' + flagValue + ' '          
        for arg in args:
            cmdStr += arg + ' '
        #remove any extra space
        if cmdStr:
            cmdStr[:-1]
        return cmdStr

