import os
import re
import sys
from functools import partial
from AnyQt.QtCore import QThread, pyqtSignal, Qt
from Orange.widgets import widget, gui, settings
from orangebiodepot.util.DockerClient import DockerClient, PullImageThread
from PyQt5 import QtWidgets, QtGui

from AnyQt.QtWidgets import (
    QWidget, QButtonGroup, QGroupBox, QRadioButton, QSlider,
    QDoubleSpinBox, QComboBox, QSpinBox, QListView, QLabel,
    QScrollArea, QVBoxLayout, QHBoxLayout, QFormLayout,
    QSizePolicy, QApplication, QCheckBox
)
#dummy classes
class ContainerPaths():
    pass
    
class BwbGuiElements():
    #keep track of List of Gui elements
    def __init__(self,required=None,active=None):
        self._dict={}
        self.required=required
        self.active={}
    def add(self,attr,guiElement):
        if attr not in self._dict:
            self._dict[attr]=[]
        self._dict[attr].append(guiElement)
        #sys.stderr.write('add item {} {} {}\n'.format(attr,guiElement,self._dict[attr]))
        
    def disable(self,attr): #gray out to disable
        if attr in self._dict:
            for g in self._dict[attr]:
                #sys.stderr.write('ditem {}\n'.format(g))
                g.setDisabled(True)

    def enable(self,attr,value):
        if attr in self._dict:
            for g in self._dict[attr]:
                #sys.stderr.write('en-item {}\n'.format(g))
                g.setEnabled(True)
                if hasattr(g,'isChecked'): #checkbox
                    if type(value) == bool:
                        g.setChecked(value)
                    else:
                        g.setChecked(True)
                elif hasattr(g,'clear'):
                    g.clear()
    
    def checkTheBox(self, attr,state=True): 
        #if there is a checkBox associated with it then we check or uncheck it - gui.spinBox and checkboxes
        if attr in self._dict[attr]:
            for g in self._dict[attr]:
                if hasattr(g,'isChecked'):
                    g.isChecked(state)
                elif hasattr(g,'checked'):
                    if g.checked is not None:
                        g.checked=state
                
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
        if slot in self._dict and connectionId not in self._dict[slot]:
            self._dict[slot].append(connectionId)
        else:
            self._dict[slot]=[connectionId]
            
    def remove (self, slot, connectionId=None):
        if slot in self._dict:
            try:
                self._dict[slot].remove(connectionId)
            except ValueError:
                pass
            #need to update the volumes
        
    def isConnected (self, slot):
        if self._dict is None:
            return False
        if (slot in self._dict) and (self._dict[slot]):
            sys.stderr.write('slot {} dict {}\n'.format(slot,self._dict[slot]))
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
        self.dockerClient.findVolumeMapping()
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
        #for container paths
        self.con=ContainerPaths()
        #gui element sub-attributes
        self.bgui=BwbGuiElements()
        #self.button=BwbGuiElement()
        #self.checkbox=BwbGuiElement()
        #self.checkboxValue=BwbGuiValue()
        #self.spin=BwbGuiElement()
        self.inputConnections=ConnectionDict(self.inputConnectionsStore)
        #volume mapping


    def initVolumes(self):
        #initializes container volumes
        for mapping in self.data['volumeMappings']:
            bwbVolAttr=mapping['attr']
            if not hasattr(self, bwbVolAttr) :
                setattr(self,bwbVolAttr,None)
            bwbVol=getattr(self,bwbVolAttr)
            setattr(self.con, bwbVolAttr, mapping['conVolume'])
                
    def bwbPathToContainerPath(self, path, isFile=False,returnNone=False):
        #converts the path entered relative to the bwb mountpoint
        #will return None if not found or the original path (default) depending on variable
        #first map it to the  path
        pathFile=None
        if isFile:
            dirPath=os.path.dirname(path)
            pathFile=os.path.basename(path)
            hostPath=os.path.normpath(self.dockerClient.to_host_directory(dirPath,returnNone=False))
        else:    
            hostPath=os.path.normpath(self.dockerClient.to_host_directory(path,returnNone=False))
        conPath=None
        #now get all the possible submappings to volumeMappings by comparing the true hostmappings
        #if submap is found convert the common path to the container path
        #return shortest path
        for mapping in self.data['volumeMappings']:
            conVol = mapping['conVolume']  
            bwbVol = self.hostVolumes[conVol]
            hostVol=os.path.normpath(self.dockerClient.to_host_directory(bwbVol,returnNone=False))
            prefix=None
            if hostVol == hostPath:
                prefix=""
            elif len(hostVol) < len(hostPath):
                prefix=os.path.commonpath([hostVol,hostPath])
            if prefix is not None :
                cleanConVol=os.path.normpath(conVol)
                myConPath= os.path.normpath(str.join(os.sep,(cleanConVol, prefix)))
                if conPath is None or len(myConPath) < len(conPath):
                    conPath=myConPath
#            sys.stderr.write("conVol {} hostVol {} path {} host path {} myConPath {} conPath {}\n".format(conVol,hostVol,path,hostPath,myConPath,conPath))
        
        if conPath is not None :
            if isFile:
                return os.path.normpath(str.join(os.sep,(conPath, pathFile)))
            return conPath
        else:
            if returnNone:
                return conPath
            return path
                    
            
    def checkRequiredParms(self):
        for parm in self.data['requiredParameters']:
            if hasattr(self,parm):
                if not getattr (self,parm):
                    return parm
        return None
        
    def getRequiredVols(self):
        for mapping in self.data['volumeMappings']:
            conVol=mapping['conVolume']
            attr=mapping['attr']
            if not hasattr(self,attr):
                return conVol
            if not getattr(self,attr):
                return conVol
            bwbVol=getattr(self,attr)
            if not bwbVol and 'default' in mapping:
                bwbVol= mapping['default']
            if self.data['parameters'][attr]['type'] is 'file':
                bwbVol=os.path.dirname(os.path.normpath(bwbVol))
            elif self.data['parameters'][attr]['type'] is 'files':
                files=re.split(r'[ ,;]',bwbVol)
                bwbVol=os.path.dirname(os.path.normpath(files[0]))
            else:
               bwbVol=os.path.normpath(bwbVol)
            self.hostVolumes[conVol]=bwbVol
        return None
 
    """
    GUI elements
    """
    def drawGUI(self):
        self.scroll_area = QScrollArea(
            verticalScrollBarPolicy=Qt.ScrollBarAlwaysOn
        )
        self.bigBox=gui.widgetBox(self.controlArea)
        self.scroll_area.setWidget(self.bigBox)
        self.scroll_area.setWidgetResizable(True)
        self.controlArea.layout().addWidget(self.scroll_area)
        
        consoleBox = gui.widgetBox(self.bigBox, "Status")
        self.infoLabel = gui.widgetLabel(consoleBox, 'Waiting...')
        self.infoLabel.setWordWrap(True)
        self.requiredBox = gui.widgetBox(self.bigBox, "Required parameters")
        self.optionalBox = gui.widgetBox(self.bigBox, "Optional parameters")
        self.drawRequiredElements()
        self.drawOptionalElements()
        controlBox = gui.vBox(self.bigBox, "Execution controls")
        btnRun = gui.button(controlBox, self, "Start", callback=self.OnRunClicked)
        css = '''
        QPushButton {background-color: #1588c5; color: white; height: 20px; border: 1px solid #1a8ac6; border-radius: 2px;}
        QPushButton:pressed { background-color: #158805; border-style: inset;}
        QPushButton:disabled { background-color: lightGray; border: 1px solid gray; }
        QPushButton:hover {background-color: #1588f5; }
        ''' 
        btnRun.setStyleSheet(css)
        btnRun.setFixedSize(50,20)

    def OnRunClicked(self):
        self.startJob()
                
    def drawRequiredElements(self):
        for pname in self.data['requiredParameters']:
            pvalue=self.data['parameters'][pname]
            if not hasattr(self,pname):
                setattr(self,pname,None)
            if (getattr(self,pname) is None) and ('default' in pvalue):
                setattr(self,pname,pvalue['default'])
        for pname in self.data['requiredParameters']:
            pvalue=self.data['parameters'][pname]
            if ('gui' in pvalue and pvalue['gui'] is not 'FileDir') or ( pvalue['type'] is  not 'file' and pvalue['type'] is  not 'directory' and pvalue['type'] is  not 'files'):
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
                self.drawLedit(pname,pvalue,self.optionalBox,layout=self.leditOptionalLayout,addCheckbox=True)           
                
        for pname in self.data['parameters']:
            if pname not in self.data['requiredParameters']:
                pvalue=self.data['parameters'][pname]                
                if ('gui' in pvalue and pvalue['gui'] is not 'Spin') or (pvalue['type'] is  not 'int'):
                    continue
                self.drawSpin(pname,pvalue,self.optionalBox,addCheckbox=True)
                
        for pname in self.data['parameters']:
            if pname not in self.data['requiredParameters']:
                pvalue=self.data['parameters'][pname]                
                if ('gui' in pvalue and pvalue['gui'] is not 'bool') or (pvalue['type'] is  not 'bool'):
                    continue
                self.drawCheckbox(pname,pvalue,self.optionalBox)

    def drawCheckbox(self,pname,pvalue,box=None):
        #for booleans - their value is the same as the checkbox state
        cb=gui.checkBox(box, self, pname, pvalue['label'])
        checkAttr=pname+'Checked'
        setattr(self,pname,False)
        setattr(self,checkAttr,getattr(self,pname))
        self.bgui.add(pname,cb)    
                      
    def drawSpin(self,pname,pvalue,box=None,addCheckbox=False):
        #for drawSpin - we use the origin version which already has a checkbox connected
        #TODO could change this to the same way we handle ledits with separate cbo
        cbState=False
        checkAttr=None
        if addCheckbox:
            checkAttr=pname+'Checked'
            setattr(self,checkAttr,cbState)
        mySpin=gui.spin(box, self, pname, minv=1, maxv=128, label=pvalue['label'], checked=checkAttr, disabled=cbState)
        if getattr(self,pname) is None:
            mySpin.clear()
        self.bgui.add(pname,mySpin)
        
    def drawLedit(self,pname,pvalue,box=None,layout=None,addCheckbox=False):
        checkAttr=None
        checkbox=None
        ledit=gui.lineEdit(None, self, pname,disabled=addCheckbox)
        if addCheckbox:
            checkAttr=pname+'Checked'
            setattr(self,checkAttr,False)
            checkbox=gui.checkBox(None, self,checkAttr,label=None)
            checkbox.stateChanged.connect(lambda : ledit.setEnabled(checkbox.isChecked()))
            self.bgui.add(pname,checkbox)
        self.bwbLedit(box,checkbox,ledit,layout=layout, label=pvalue['label'])
        #check if the value is none - then we clear it
        if getattr(self,pname) is None:
            ledit.clear()
        self.bgui.add(pname,ledit)

    def drawFileDirElements(self,pname, pvalue, box=None, layout=None, addCheckbox=False):
        #note that using lambda does not work in this function - truncates string variable so partial used instead
        #draw required elements'
        checkbox=None
        ledit=gui.lineEdit(None, self, pname,disabled=addCheckbox)
        button=gui.button(None, self, "", callback= partial(self.browseFileDir, attr= str(pname) ,filetype=pvalue['type']), 
                          autoDefault=True, width=19, height=19,disabled=addCheckbox)
        if getattr(self,pname) is None:
            ledit.clear()
        if addCheckbox:
            checkAttr=pname+'Checked'
            setattr(self,checkAttr,False)
            checkbox=gui.checkBox(None, self,checkAttr,label=None)
            checkbox.stateChanged.connect(lambda : ledit.setEnabled(checkbox.isChecked() and not self.inputConnections.isConnected(pname) ))
            checkbox.stateChanged.connect(lambda : button.setEnabled(checkbox.isChecked() and not self.inputConnections.isConnected(pname) ))
            self.bgui.add(pname,checkbox)
        self.bwbFileEntry(box,button,ledit,layout=layout,label=pvalue['label']+':', entryType=pvalue['type'], checkbox=checkbox)
        self.bgui.add(pname,ledit)
        self.bgui.add(pname,button)
 
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
        colNum=1
        if(checkbox):
            layout.addWidget(checkbox,layout.nextRow,0)    
        layout.addWidget(QtGui.QLabel(label),layout.nextRow,colNum)
        layout.addWidget(ledit,layout.nextRow,colNum+1)
        if layout.nextRow == 1:
            widget.layout().addLayout(layout)
        layout.nextRow+=1
    
    def browseFileDir(self, attr, filetype=None):
        defaultDir = '/root'
        if os.path.exists('/data'):
            defaultDir = '/data'
        if filetype is 'file':
            myFile=QtWidgets.QFileDialog.getOpenFileName(self, "Locate file", defaultDir)[0]
            setattr(self,attr,myFile)
            dirAttr=attr+"Dir"
            setattr(self,dirAttr,os.path.dirname(myFile))
            self.setDirectoriesAttr(dirAttr)
        elif filetype is 'files':
            myFiles=QtWidgets.QFileDialog.getOpenFileNames(self, "Locate file", defaultDir)
            if myFiles:
                setattr(self,attr, ' '.join(myFiles[0]))
            else:
                setattr(self,attr,None)
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
        if not hasattr(self,attr):
            setattr(self,attr,None)
        if not hasattr(self.con,attr):
            setattr(self.con,attr,None) 
        containerPath=getattr(self.con,attr)
        hostPath=getattr(self,attr)
        if hostPath is None:
            self._hostDirectories[containerPath] = None
        else:
            self._hostDirectories[containerPath] = hostPath.strip()
    
    def handleInputs(self, value, attr, sourceId=None):
        if value is None:
            self.inputConnections.remove(attr,sourceId)
            setattr(self,attr,None)
        else:
            self.inputConnections.add(attr,sourceId)
            setattr(self,attr,value)
        self.updateGui(attr,value)
        
    def updateGui(self,attr,value):
        if self.inputConnections.isConnected(attr):
            sys.stderr.write('disabling {}\n'.format(attr))
            self.bgui.disable(attr)
        else:
            sys.stderr.write('enabling {} with {}\n'.format(attr,value))
            self.bgui.enable(attr,value)
            
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
        
    def generateCmdFromData (self):
        flags={}
        args=[]
        for pname, pvalue in self.data['parameters'].items():
            #if required or checked then it is added to the flags
            addFlag=False
            checkattr=pname+'Checked'
            if hasattr(self,checkattr) and getattr(self,checkattr):
                addFlag=True
            if pname in self.data['requiredParameters']:
                sys.stderr.write('{} is required\n'.format(pname))
                addFlag=True
            else:
                sys.stderr.write('{} is not required\n'.format(pname))
            if addFlag:
                if pvalue['flags']:
                    if pvalue['type'] is 'bool':
                        flags[pvalue['flags'][0]] = None
                    elif pvalue['type'] is 'file':
                        filename=str(getattr(self,pname))
                        if filename:
                            hostFilename=self.bwbPathToContainerPath(filename, isFile=True,returnNone=False)
                            flags[pvalue['flags'][0]]=str(hostFilename)
                    elif pvalue['type'] is 'files':
                        files=re.split(r'[ ,;]',getattr(self,pname))
                        if files:
                            hostFiles=[]
                            for f in files:
                                hostFiles.append(self.bwbPathToContainerPath(filename, isFile=True,returnNone=False))
                            flags[pvalue['flags'][0]]=' '.join(hostFiles)                    
                    elif pvalue['type'] is 'directory':
                        path=str(getattr(self,pname))
                        if path:
                            hostPath=self.bwbPathToContainerPath(path, returnNone=False)
                            flags[pvalue['flags'][0]]=str(hostPath)
                    else:                        
                        flags[pvalue['flags'][0]]=str(getattr(self,pname))
                else:
                    if pvalue['type'] is 'file':
                        filename=str(getattr(self,pname))
                        if filename:
                            hostFilename=self.bwbPathToContainerPath(filename, isFile=True,returnNone=False)
                            args.append(hostFilename)
                    elif pvalue['type'] is 'files':
                        files=re.split(r'[ ,;]',getattr(self,pname))
                        for f in files:
                            args.append(self.bwbPathToContainerPath(f, isFile=True,returnNone=False))                   
                    elif pvalue['type'] is 'directory':
                        path=str(getattr(self,pname))
                        if path:
                            hostPath=self.bwbPathToContainerPath(path, returnNone=False)
                            args.append(hostpath)
                    else:
                        args.append(str(getattr(self,pname)))
        return self.generateCmdFromBash([self.data['command']],flags,args=args)
                
    def generateCmdFromBash (self, executables = [], flags = {}, args = []):
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

