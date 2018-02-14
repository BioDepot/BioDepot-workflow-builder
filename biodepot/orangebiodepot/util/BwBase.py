import os
from functools import partial
from AnyQt.QtCore import QThread, pyqtSignal, Qt
from Orange.widgets import widget, gui, settings
from orangebiodepot.util.DockerClient import DockerClient, PullImageThread
from PyQt5 import QtWidgets, QtGui

#dummy classes
class ContainerPaths():
    pass
    
class BwbGuiElement():
    pass

class BwbGuiValue():
    pass      
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
    defaultFileIcon=QtGui.QIcon('/biodepot/orangebiodepot/icons/bluefile.png')

    def __init__(self, image_name, image_tag):
        super().__init__()
        self._hostDirectories = {}
        self._dockerImageName = image_name
        self._dockerImageTag = image_tag
        self._dockerVolumes = None
        self._dockerCommands = None
        self._dockerEnvironments = None
        self._Flag_isRunning = False
        #drawing layouts for gui
        #file directory 
        self.fileDirRequiredLayout=QtGui.QGridLayout()
        self.fileDirRequiredLayout.setSpacing(5)
        self.fileDirOptionalLayout=QtGui.QGridLayout()
        self.fileDirOptionalLayout.setSpacing(5)
        setattr(self.fileDirOptionalLayout,'nextRow',1)
        setattr(self.fileDirRequiredLayout,'nextRow',1)
        #lineEdits
        self.leditRequiredLayout=QtGui.QGridLayout()
        self.leditRequiredLayout.setSpacing(5)
        self.leditOptionalLayout=QtGui.QGridLayout()
        self.leditOptionalLayout.setSpacing(5)
        setattr(self.leditOptionalLayout,'nextRow',1)
        setattr(self.leditRequiredLayout,'nextRow',1)
    """
    GUI elements
    """
    def drawGUI(self):
        self.nextRow=1
        consoleBox = gui.widgetBox(self.controlArea, "Status")
        self.infoLabel = gui.widgetLabel(consoleBox, 'Waiting...')
        self.infoLabel.setWordWrap(True)
        self.requiredBox = gui.widgetBox(self.controlArea, "Required parameters")
        self.optionalBox = gui.widgetBox(self.controlArea, "Optional parameters")
        self.drawRequiredElements()
        self.drawOptionalElements()
        
        controlBox = gui.vBox(self.controlArea, "Execution controls")
        btnRun = gui.button(controlBox, self, "Run", callback=self.OnRunClicked)
        
    def drawRequiredElements(self):
        for pname in self.data['requiredParameters']:
            pvalue=self.data['parameters'][pname]
            if not hasattr(self,pname):
                setattr(self,pname,None)
            if (getattr(self,pname) is None) and ('default' in pvalue):
                setattr(self,pname,pvalue['default'])
        for pname in self.data['requiredParameters']:
            pvalue=self.data['parameters'][pname]                
            if ('gui' in pvalue and pvalue['gui'] is not 'FileDir') or (pvalue['type'] is  not 'file') and (pvalue['type'] is  not 'directory'):
                continue
            self.drawFileDirElements(pname, pvalue, box=self.requiredBox,layout=self.fileDirRequiredLayout)
                                        
        for pname in self.data['requiredParameters']:
            pvalue=self.data['parameters'][pname]                
            if ('gui' in pvalue and pvalue['gui'] is not 'Ledit') or (pvalue['type'] is  not 'double') and (pvalue['type'] is  not 'text'):
                continue
            self.drawLedit(pname,pvalue,self.requiredBox,layout=self.leditRequiredLayout)           
            
        for pname in self.data['requiredParameters']:
            pvalue=self.data['parameters'][pname]                
            if ('gui' in pvalue and pvalue['gui'] is not 'Spin') or (pvalue['type'] is  not 'int'):
                continue
            self.drawSpin(pname,pvalue,self.requiredBox)
            
        for pname in self.data['requiredParameters']:
            pvalue=self.data['parameters'][pname]                
            if ('gui' in pvalue and pvalue['gui'] is not 'bool') or (pvalue['type'] is  not 'bool'):
                continue
            self.drawCheckbox(pname,pvalue,self.requiredBox)
                                
    def drawOptionalElements(self):
        for pname in self.data['parameters']:
            if pname not in self.data['requiredParameters']:
                pvalue=self.data['parameters'][pname]
                if not hasattr(self,pname):
                    setattr(self,pname,None)
                if (getattr(self,pname) is None) and ('default' in pvalue):
                    setattr(self,pname,pvalue['default'])
                
        for pname in self.data['parameters']:
            if pname not in self.data['requiredParameters']:
                pvalue=self.data['parameters'][pname]                
                if ('gui' in pvalue and pvalue['gui'] is not 'FileDir') or (pvalue['type'] is  not 'file') and (pvalue['type'] is  not 'directory'):
                    continue
                self.drawFileDirElements(pname, pvalue, box=self.optionalBox,layout=self.fileDirOptionalLayout,addCheckbox=True)
                                            
        for pname in self.data['parameters']:
            if pname not in self.data['requiredParameters']:
                pvalue=self.data['parameters'][pname]                
                if ('gui' in pvalue and pvalue['gui'] is not 'Ledit') or (pvalue['type'] is  not 'double') and (pvalue['type'] is  not 'text'):
                    continue
                self.drawLedit(pname,pvalue,self.optionalBox,layout=self.leditOptionalLayout)           
                
        for pname in self.data['parameters']:
            if pname not in self.data['requiredParameters']:
                pvalue=self.data['parameters'][pname]                
                if ('gui' in pvalue and pvalue['gui'] is not 'Spin') or (pvalue['type'] is  not 'int'):
                    continue
                self.drawSpin(pname,pvalue,self.optionalBox)
                
        for pname in self.data['parameters']:
            if pname not in self.data['requiredParameters']:
                pvalue=self.data['parameters'][pname]                
                if ('gui' in pvalue and pvalue['gui'] is not 'bool') or (pvalue['type'] is  not 'bool'):
                    continue
                self.drawCheckbox(pname,pvalue,self.optionalBox)

    def drawCheckbox(self,pname,pvalue,box=None):
        cb=gui.checkBox(box, self, pname, pvalue['label'])
        setattr(self.checkbox,pname,cb)
                      
    def drawSpin(self,pname,pvalue,box=None):
        checkattr=pname+'Checked'
        if not hasattr(self,checkattr):
            setattr(self,checkattr,None)
        if getattr(self,checkattr) is None:
            setattr(self,checkattr, False)
        
        mySpin=gui.spin(box, self, pname, minv=1, maxv=128, label=pvalue['label'], checked=checkattr)
        setattr(self.spin,pname,mySpin)
        
    def drawLedit(self,pname,pvalue,box=None,layout=None):
        checkattr=pname+'Checked'
        if not hasattr(self,checkattr):
            setattr(self,checkattr,None)
        if getattr(self,checkattr) is None:
            setattr(self,checkattr, False)
        setattr(self.checkbox,pname,gui.checkBox(None, self, checkattr,label=None))
        checkbox=getattr(self.checkbox,pname)
        setattr(self.ledit,pname,gui.lineEdit(None, self, pname,disabled=not checkbox.isChecked()))
        ledit=getattr(self.ledit,pname)
        checkbox.stateChanged.connect(lambda : ledit.setEnabled(checkbox.isChecked()))
        
        self.bwbLedit(box,checkbox,ledit,layout=layout, label=pvalue['label'])
               
    def drawFileDirElements(self,pname, pvalue, box=None, layout=None, addCheckbox=None):
        #note that using lambda does not work in this function - truncates string variable so partial used instead
        #draw required elements'
        checkbox=None
        disabled=False
        if addCheckbox:
            disabled=True
            checkattr=pname+'Checked'
            if not hasattr(self,checkattr):
                setattr(self,checkattr,None)
            if getattr(self,checkattr) is None:
                setattr(self,checkattr, False)
            setattr(self.checkbox,pname,gui.checkBox(None, self, checkattr,label=None))
            checkbox=getattr(self.checkbox,pname)
            disabled=not checkbox.isChecked()

        if pname in self.data['inputs']:
            setattr(self.ledit,pname,gui.lineEdit(None, self, pname,
                    disabled=self.inputConnections.isConnected(pname)))
            disabled=self.inputConnections.isConnected(pname)
        else:
            setattr(self.ledit,pname,gui.lineEdit(None, self, pname,disabled=disabled))
            
        setattr(self.button, pname, gui.button(None, self, "", 
                    callback= partial(self.browseFileDir, attr= str(pname) ,filetype=bool(pvalue['type'] is 'file')), autoDefault=True, width=19, height=19,disabled=disabled))
        ledit=getattr(self.ledit,pname)
        button=getattr(self.button,pname)
        if checkbox:
            checkbox.stateChanged.connect(lambda : ledit.setEnabled(checkbox.isChecked() and not self.inputConnections.isConnected(pname) ))
            checkbox.stateChanged.connect(lambda : button.setEnabled(checkbox.isChecked() and not self.inputConnections.isConnected(pname) ))
        self.bwbFileEntry(box,button,ledit,layout=layout,label=pvalue['label']+':', entryType=pvalue['type'], checkbox=checkbox)
 
    def bwbFileEntry(self, widget, button, ledit, icon=defaultFileIcon,layout=None, label=None,entryType='file', checkbox=None):
        button.setStyleSheet("border: 1px solid #1a8ac6; border-radius: 2px;")
        button.setIcon(icon)
        ledit.setClearButtonEnabled(True)
        ledit.setPlaceholderText("Enter {}".format(entryType))
        col=0
        if checkbox:
            layout.addWidget(checkbox,layout.nextRow,col)
            col+=1     
        if label:
            #myLabel=QtGui.QLabel(label)
            #myLabel.setAlignment(Qt.AlignLeft | Qt.AlignVCenter) 
            layout.addWidget(QtGui.QLabel(label),layout.nextRow,col)
        layout.addWidget(ledit,layout.nextRow,col+1)
        layout.addWidget(button,layout.nextRow,col+2)
        if layout.nextRow == 1:
            widget.layout().addLayout(layout)
        layout.nextRow+=1
        
    def bwbLedit(self, widget,checkbox ,ledit, layout=None, label=None):
        ledit.setClearButtonEnabled(True)
        ledit.setPlaceholderText("Enter parameter")
        layout.addWidget(checkbox,layout.nextRow,0)    
        layout.addWidget(QtGui.QLabel(label),layout.nextRow,1)
        layout.addWidget(ledit,layout.nextRow,2)
        
        if layout.nextRow == 1:
            widget.layout().addLayout(layout)
        layout.nextRow+=1
    
    def browseFileDir(self, attr, filetype=False):
        defaultDir = '/root'
        if os.path.exists('/data'):
            defaultDir = '/data'
        if filetype:
            myFile=QtWidgets.QFileDialog.getOpenFileName(self, "Locate file", defaultDir)[0]
            setattr(self,attr,myFile)
            dirAttr=attr[:-4]+"Dir"
            setattr(self,dirAttr,os.path.dirname(myFile))
            self.setDirectoriesAttr(dirAttr)
        else:
            myDir=QtWidgets.QFileDialog.getExistingDirectory(self, caption="Locate directory", directory=defaultDir)
            setattr(self,attr,myDir)
            self.setDirectoriesAttr(attr)

    def Event_OnRunFinished(self):
        raise Exception('Event_OnRunFinished not implemented!')

    def Event_OnRunMessage(self, message):
        raise Exception('Event_OnRunMessage not implemented!')

    def Initialize_InputKeys(self, keys):
        for k in keys:
            self._hostDirectories[k] = None

    def setDirectoriesAttr(self, attr):
        containerPath=getattr(self.con,attr)
        hostPath=getattr(self,attr)
        if hostPath is None:
            self._hostDirectories[containerPath] = None
        else:
            self._hostDirectories[containerPath] = hostPath.strip()
        
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

    def setDir (self, path, attr, sourceId=None):
        #generic function for inputs
        ledit=getattr(self.ledit,attr)
        button=getattr(self.button,attr)
        if path is None:
            self.inputConnections.remove(attr,sourceId)
            setattr(self,attr,None)
            ledit.clear()
        else:
            setattr(self,attr,path)
            self.inputConnections.add(attr,sourceId)
        self.setDirectoriesAttr(attr)
        button.setEnabled(not self.inputConnections.isConnected(attr))
        ledit.setEnabled(not self.inputConnections.isConnected(attr))
        
    def setFile (self, path, attr, sourceId=None):
        #generic function for inputs
        ledit=getattr(self.ledit,attr)
        button=getattr(self.button,attr)
        dirAttr=attr[:-4]+'Dir'
        if path is None:
            self.inputConnections.remove(attr,sourceId)
            setattr(self,attr,None)
            setattr(self,dirAttr,None)
            ledit.clear()
        else:
            setattr(self,attr,path)
            setattr(self,dirAttr,os.path.dirname(path))
            self.inputConnections.add(attr,sourceId)
        self.setDirectoriesAttr(dirAttr)
        button.setEnabled(not self.inputConnections.isConnected(attr))
        ledit.setEnabled(not self.inputConnections.isConnected(attr))
          
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

