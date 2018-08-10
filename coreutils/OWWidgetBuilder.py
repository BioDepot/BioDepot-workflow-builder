import os
import re
import sys
import json
import jsonpickle
import pickle
import csv
import tempfile
import OWImageBuilder
from pathlib import Path
from shutil import copyfile
from orangebiodepot.util.createWidget import mergeWidget, createWidget, findDirectory, findIconFile
from copy import deepcopy
from collections import OrderedDict
from functools import partial
from AnyQt.QtCore import QThread, pyqtSignal, Qt
from Orange.widgets import widget, gui, settings
from orangebiodepot.util.DockerClient import DockerClient, PullImageThread, ConsoleProcess
from PyQt5 import QtWidgets, QtGui
from PyQt5.QtWidgets import QInputDialog, QLineEdit, QMainWindow, QApplication, QPushButton, QWidget, QAction, QTabWidget, QVBoxLayout, QMessageBox

from AnyQt.QtWidgets import (
    QWidget, QButtonGroup, QGroupBox, QRadioButton, QSlider,
    QDoubleSpinBox, QComboBox, QSpinBox, QListView, QLabel,
    QScrollArea, QVBoxLayout, QHBoxLayout, QFormLayout,
    QSizePolicy, QApplication, QCheckBox
)
defaultIconFile='/biodepot/Bwb_core/icons/default.png'
class tabbedWindow(QTabWidget):
    def __init__(self, parent = None):
        super(tabbedWindow, self).__init__(parent)
        
    def add(self,title,minHeight=240):
        tab=QWidget()
        tab.setMinimumHeight(minHeight)
        self.addTab(tab,title)
        box=gui.widgetBox(tab)
        tab.setLayout(QVBoxLayout())
        tab.layout().addWidget(self.getScrollArea(box))    
        return self.getLeditLayout(box)
        
    def addBox(self,title,minHeight=240):
        tab=QWidget()
        tab.setMinimumHeight(minHeight)
        self.addTab(tab,title)
        box=gui.widgetBox(tab)
        tab.setLayout(QVBoxLayout())
        tab.layout().addWidget(self.getScrollArea(box))    
        return box, self.getLeditLayout(box)
               
    def getScrollArea(self, box):
        scroll_area = QScrollArea(
            verticalScrollBarPolicy=Qt.ScrollBarAlwaysOn
        )
        scroll_area.setWidget(box)
        scroll_area.setWidgetResizable(True)
        return scroll_area
    
    def getLeditLayout(self,box):
        layout=QtGui.QGridLayout()
        layout.setSpacing(5)
        setattr(layout,'nextRow',1)
        box.layout().addLayout(layout)
        return layout
        
class DragAndDropList(QtGui.QListWidget):
    #overloads the Drag and dropEvents to emit code
    itemMoved = pyqtSignal(int, int) # Old index, new index, item
    def __init__(self, parent=None, **args):
        super(DragAndDropList, self).__init__(parent, **args)
        self.setAcceptDrops(True)
        self.setDragEnabled(True)
        self.setDragDropMode(QtGui.QAbstractItemView.InternalMove)
        self.drag_item = None
        self.drag_row = None

    def dropEvent(self, event):
        super(DragAndDropList, self).dropEvent(event)
        self.itemMoved.emit(self.drag_row, self.currentRow())
        self.drag_item = None
     
    def startDrag(self, supportedActions):
        self.drag_row = self.currentRow()
        super(DragAndDropList, self).startDrag(supportedActions)
        
class WidgetItem():
    #each widget item has:
    #A set of gui widgets that correspond to the name and each parameter
    #A set of state vectors corresponding to the different elements that make up the item
        
    def __init__(self, guiElements=OrderedDict(), states={}):
        #guielements has parameter as key and qwidget as value
        self.guiElements=guiElements
        self.states=states
        if self.states and self.guiElements:
            for (key,element) in self.guiElements.items():
                element.setState(self.states[key])
        elif self.guiElements:
            self.blankState()
                
    def printValue(self):
        output='{} : {{'.format(self.guiElements['name'].getValue())
        for key, element in self.guiElements.items():
            if key != 'name' and element.getValue() is not None:
                output+=' {} : {}, '.format(key,element.getValue())
        output=output[:-2]+' }'
        return output
    
    def getStateValue(self,states):
        name=None
        myDict={}
        i=0
        if states:
            for key, element in self.guiElements.items():
                if key == 'name' or key == 'Name':
                    name=element.getStateValue(states[i])
                else:
                    myDict[key]=element.getStateValue(states[i])
                i+=1    
            return (name,myDict)
        return None
        
    def getState(self):
        serialStates=[]
        for key, element in self.guiElements.items():
            self.states[key]=element.getState()
            serialStates.append(self.states[key])
        return serialStates
        
    def setState(self,serialStates):
        if serialStates is None:
            self.blankState()
        else:
            i=0
            for key, element in self.guiElements.items():
                self.states[key]=serialStates[i]
                element.setState(self.states[key])
                i+=1
            
    def blankState(self):
        for key, element in self.guiElements.items():
            element.setState(None)
            self.states[key]=None

class OWWidgetBuilder(widget.OWWidget):
    name = "Widget Builder"
    description = "Build a new widget from a set of bash commands and a Docker container"
    category = "Bwb-core"
    icon = "icons/build.png"
    priority = 2

    inputs = []
    outputs = []
    pset=partial(settings.Setting,schema_only=True)
    
    want_main_area = False
    want_control_area = True
    allStates=pset({})
    allAttrs=pset({})
    data={}
    
    def register(self,widgetPath,widgetName):
        sys.stderr.write('register with path {} name {}\n'.format(widgetPath,widgetName))
        jsonFile="{}/{}.json".format(widgetPath,widgetName)
        directory=findDirectory(jsonFile)
        destPath='/widgets/{}'.format(widgetName);
        
        if os.path.realpath(widgetPath) != os.path.realpath(destPath):
            #check if it exists
            if os.path.exists(destPath):
                qm = QtGui.QMessageBox
                ret=qm.question(self,'', "{} exists - OverWrite ?".format(widgetName), qm.Yes | qm.No)
                if ret == qm.No:
                    return
                os.system("rm {} -rf ".format(destPath))    
            os.system("cp -r {} {}".format(widgetPath,destPath))
        #make linkages
        os.system ("ln -sf  /widgets/{}/{}.py /biodepot/{}/OW{}.py".format(widgetName,widgetName,directory,widgetName))
        
    def drawExec(self, layout=None):
        css = '''
        QPushButton {background-color: #1588c5; color: white; height: 20px; border: 1px solid black; border-radius: 2px;}
        QPushButton:hover {background-color: #1555f5; }
        QPushButton:hover:pressed { background-color: #1588c5; color: black; border-style: inset; border: 1px solid white} 
        QPushButton:disabled { background-color: lightGray; border: 1px solid gray; }
        '''
        #pname='github'
        #githubLedit=self.makeLedit(pname,'Enter directory', label='Bwb directory:')
        #self.initAllStates(pname,githubLedit)
        #githubLedit.ledit.textChanged.connect(lambda: self.updateAllStates(pname,githubLedit,githubLedit.getState()))
        #layout.addWidget(githubLedit.label,layout.nextRow,0)
        #layout.addWidget(githubLedit.ledit,layout.nextRow,1,1,1)
        #button=gui.button(None, self, "", callback= partial(self.browseFileDir, attr='gitHub',fileType='Directory'),autoDefault=True, width=19, height=19)
        #button.setIcon(self.browseIcon)
        #layout.addWidget(button,layout.nextRow,2)
        # self.saveJsonBtn = gui.button(None, self, "Save json", callback=self.saveJson)
        # self.saveJsonBtn.setStyleSheet(css)
        # self.saveJsonBtn.setFixedSize(80,20)

        self.saveWidgetBtn = gui.button(None, self, "Save", callback=self.saveWidget)
        self.saveWidgetBtn.setStyleSheet(css)
        self.saveWidgetBtn.setFixedSize(50,20)
        self.saveWidgetAsBtn = gui.button(None, self, "Save as", callback=self.saveWidgetAs)
        self.saveWidgetAsBtn.setStyleSheet(css)
        self.saveWidgetAsBtn.setFixedSize(70,20)
        self.loadWidgetBtn = gui.button(None, self, "Load", callback=self.loadWidget)
        self.loadWidgetBtn.setStyleSheet(css)
        self.loadWidgetBtn.setFixedSize(70,20)
        self.registerBtn = gui.button(None, self, "Register", callback=self.registerWidget)
        self.registerBtn.setStyleSheet(css)
        self.registerBtn.setFixedSize(100,20)
        self.rebuildBtn = gui.button(None, self, "Rebuild", callback=self.rebuildWidget)
        self.rebuildBtn.setStyleSheet(css)
        self.rebuildBtn.setFixedSize(100,20)
        #choose save mode

        saveLabel=QtGui.QLabel('Save mode:')
        saveLabel.setAlignment(Qt.AlignRight | Qt.AlignVCenter)
        box=QtGui.QHBoxLayout()
        box.addWidget(saveLabel)
        box.addWidget(self.saveMode)
        box.addWidget(self.saveWidgetBtn)
        box.addWidget(self.saveWidgetAsBtn)
        box.addWidget(self.loadWidgetBtn)
        box.addWidget(self.registerBtn)
        box.addWidget(self.rebuildBtn)
        box.addStretch(1)
        layout.addLayout(box)

    def onSaveModeChange(self,widget):
        self.saveModeIndex=widget.currentIndex()
        self.allAttrs['saveModeIndex']=self.saveModeIndex
        
    def loadWidget(self,loadWidgetDir=None,loadNameCheck=True):
        if self.widgetDir:
            startDir=self.widgetDir
        else:
            startDir='/templates/Generic'
        if not loadWidgetDir:
            loadWidgetDir=QtWidgets.QFileDialog.getExistingDirectory(self, caption="Choose widget to load", directory=startDir)
        if loadWidgetDir:
            widgetName=os.path.split(loadWidgetDir)[-1]
            allAttrsFile="{}/{}.attrs".format(loadWidgetDir,widgetName)
            allStatesFile="{}/{}.states".format(loadWidgetDir,widgetName)
            if os.path.exists(allAttrsFile):
                self.allAttrs=self.unPickleData(allAttrsFile)
            if os.path.exists(allStatesFile):
                self.allStates=self.unPickleData(allStatesFile)
            #we need to make default files in case the default files change from versions
            if loadNameCheck or not self.widgetName:
                self.widgetName=self.getWidgetName()
            self.makeDefaultFiles()
            if not self.widgetName:
                return
            if self.isDrawn:
                self.updateWidget()
            else:
                self.startWidget()
        return

    def getWidgetName(self):
        niceName=self.widgetName
        if not niceName:
            if 'name' in self.allAttrs and self.allAttrs['name']:
                niceName=self.allAttrs['name']
            else:
                niceName="generic{}".format(os.getpid())
            niceName =re.sub(r"[^\w\s]", '',niceName)
            niceName = re.sub(r"\s+", '_', niceName)
        nameDir, okPressed = QInputDialog.getText(self, "Widget Name ","Enter widget name:", QLineEdit.Normal,niceName)
        if okPressed and nameDir:
            niceName=re.sub(r"[^\w\s]", '',nameDir)
            niceName = re.sub(r"\s+", '_', nameDir)
        else:
                #return what we started with
            return self.widgetName
        return niceName
    
    def rebuildWidget(self):
        myDir=QtWidgets.QFileDialog.getExistingDirectory(self, caption="Locate Bwb directory", directory=self.defaultDir)
        (imageName,accept)=QtGui.QInputDialog.getText(self,"New image name", "Enter image name",QtGui.QLineEdit.Normal,"biodepot/bwb")
        if imageName:
            qm = QtGui.QMessageBox
            ret=qm.question(self,'', "Are you sure you want to rebuild {} from directory {} ?".format(imageName,myDir), qm.Yes | qm.No)
            if ret == qm.Yes:
                self.console=QtGui.QTextEdit()
                self.console.setReadOnly(True)
                pal=QtGui.QPalette()
                pal.setColor(QtGui.QPalette.Base,Qt.black)
                pal.setColor(QtGui.QPalette.Text,Qt.green)
                self.console.setPalette(pal)
                self.console.setAutoFillBackground(True)
                self.console.show()
                self.pConsole=ConsoleProcess(console=self.console,finishHandler=self.finishRebuild)
                cmd='rsync -av /biodepot/ {}/biodepot/ --delete '.format(myDir) 
                cmd+= '&& rsync -av /widgets/ {}/widgets/ --delete '.format(myDir)
                cmd+= '&& rsync -av /workflows/ {}/workflows/ --delete '.format(myDir)
                cmd+= '&& docker build -t {} {}'.format(imageName,myDir)
                self.registerBtn.setEnabled(False)
                self.rebuildBtn.setEnabled(False)
                self.saveWidgetBtn.setEnabled(False)
                self.saveWidgetAsBtn.setEnabled(False)
                self.loadWidgetBtn.setEnabled(False)
                self.pConsole.process.start('/bin/bash',['-c',cmd])
    def finishRebuild(self,stopped=None):
        self.registerBtn.setEnabled(True)
        self.rebuildBtn.setEnabled(True)
        self.saveWidgetBtn.setEnabled(True)
        self.saveWidgetAsBtn.setEnabled(True)
        self.loadWidgetBtn.setEnabled(True)        
        #self.pcconsole.close()
        
    def buildData(self):
        myData={}
        self.data['name']=self.widgetName
        if 'category' not in self.data:
            self.data['category']='User'
        for attr in ('name','description','category','docker_image_name','docker_image_tag',
            'priority','icon','inputs','outputs','volumes','parameters','command','autoMap'):
            if attr in self.data and self.data[attr]:
                myData[attr]=deepcopy(self.data[attr])
            else:
                myData[attr]=None
        #separate multiple commands and give default command
        if 'command' in myData and myData['command']: 
            if '\n' in myData['command']:
                myData['command']=myData['command'].split('\n')
            else:
                myData['command']=[myData['command']]

        
        #add persistance -all for now - add option to for transient data later
        myData['persistentSettings']='all'
        
        #add required elements
        reqList=[]
        if 'parameters' in myData and myData['parameters']:
            for pname,pvalue in myData['parameters'].items():
                if not pvalue['optional']:
                    reqList.append(pname)
        myData['requiredParameters']=reqList
        #add volume mappings
        if 'volumes' in myData and myData['volumes']:
            myMappings=[]
            for attr, volume in myData['volumes'].items():
                myMappings.append({'conVolume':volume['containerVolume'], 'attr' : attr})
            myData['volumeMappings']=myMappings
            myData.pop('volumes',None)
        #replace text str with type(str)
        for pname in ('inputs','outputs'):
            if pname in myData and myData[pname]:
                for key, myDict in myData[pname].items():
                    if myDict['type'] == 'str':
                        myDict['type'] = type('str')

        #for now just keep default, flags, 
        #also make all defaults false for booleans - will fix these kluges after cleaning up BwBase code
        if 'parameters' in myData and myData['parameters']:
            for key, myDict in myData['parameters'].items():
                newDict={}
                
                if  'default' in myDict and myDict['default'] is not None:
                    sys.stderr.write('key is {} default is {}\n'.format(key,myDict['default']))
                    #make sure that these are the correct type
                    if myDict['type'] == 'bool':
                        if myDict['default'] == 'False':
                            newDict['default']=False
                        elif myDict['default'] == 'True':
                            newDict['default']=True
                        else:
                            raise Exception ('{} is boolean - default values must be True or False not {}'.format(pname,myDict['default']))
                    #check for lists - these are treated as column delimited
                    elif 'type' in myDict and 'list' in myDict['type']:
                        reader=csv.reader([myDict['default']],skipinitialspace=True)
                        for row in reader:
                            newDict['default']=row
                            sys.stderr.write('default list is type{} {}\n'.format(type(newDict['default']),newDict['default']))
                    elif 'type' in myDict and myDict['type'] == 'int':
                        newDict['default']=int(myDict['default'])
                        sys.stderr.write('key is {} new default is {}\n'.format(key,myDict['default']))
                    elif 'type' in myDict and (myDict['type'] == 'double' or myDict['type'] == 'float') :
                        newDict['default']=float(myDict['default'])                
                    elif  'type' in myDict and (myDict['type'] == 'str' or myDict['type'] == 'file' or myDict['type'] == 'directory') :
                        newDict['default']=str(myDict['default'])
                    else:
                        newDict['default']=myDict['default']
                
                if 'flag' in myDict:
                    newDict['flag']=myDict['flag']
                #arguments are the same as having a null flag value
                if 'argument' in myDict and myDict['argument']:
                    newDict['argument'] = True
                if 'label' in myDict:
                    newDict['label']=myDict['label']
                else:
                    newDict['label']=None
                if 'type' in myDict and myDict['type']:
                    newDict['type']=myDict['type']
                if 'env' in myDict and myDict['env']:
                    newDict['env']=myDict['env']
                myData['parameters'][key]=newDict
        return myData
    def saveJson(self):
        myData=self.buildData()
        #get json filename
        dataJ=jsonpickle.encode(myData)
        jSaveFile=QtWidgets.QFileDialog.getSaveFileName(self,"QFileDialog.getSaveFileName()","","All Files (*);;json Files (*.json)")[0]
        if jSaveFile:
            sys.stderr.write('save file is {}\n'.format(jSaveFile))
            with open(jSaveFile,"w") as f:
                f.write(dataJ)
            f.close()
    
    def saveWidget(self):
        if not self.widgetDir or not self.widgetName:
            self.saveWidgetAs()
            return
        self.pickleWidget()
        title='Save {}'.format(self.widgetName)
        message='Saved widget to {}'.format(self.widgetDir)
        QtGui.QMessageBox.information(self, title,message,QtGui.QMessageBox.Ok)
                    
    def saveWidgetAs(self):
        self.widgetName=self.getWidgetName()
        if not self.widgetName:
            return
        oldWidgetName=self.widgetName
        self.widgetName=self.getWidgetName()
        self.outputDir=QtWidgets.QFileDialog.getExistingDirectory(self, caption="Choose directory to save the widget in", directory=self.defaultDir)
        if self.outputDir: 
            self.widgetDir=self.outputDir+'/'+ self.widgetName
            self.allAttrs['name']=self.widgetName
            self.makeDefaultFiles()
            self.pickleWidget()
            
            title='Save {}'.format(self.widgetName)
            message='Saved widget to {}'.format(self.widgetDir)
            QtGui.QMessageBox.information(self, title,message,QtGui.QMessageBox.Ok)
        return
            
    def makeDefaultFiles(self):
        iconDir='{}/icon'.format(self.widgetDir)
        os.system('mkdir -p {} '.format(iconDir))
        if not os.listdir('{}/icon'.format(self.widgetDir)):
            copyfile(defaultIconFile,iconDir+ '/'+os.path.basename(defaultIconFile))
        os.system('mkdir -p {}/Dockerfiles'.format(self.widgetDir))
        
    def pickleWidget(self):
        myData=self.buildData()
        Path(self.widgetDir).mkdir(parents=True,exist_ok=True)
        outputWidget="{}/{}.py".format(self.widgetDir,self.widgetName)
        #check if data only
        if self.saveMode.currentIndex() < 2:
            if self.saveMode.currentIndex() == 1:
                #merg
                mergeWidget(None,outputWidget,self.widgetName,inputData=myData)
            else:
                createWidget(None,outputWidget,self.widgetName,inputData=myData)
        allStatesFile="{}/{}.states".format(self.widgetDir,self.widgetName)
        allAttrsFile="{}/{}.attrs".format(self.widgetDir,self.widgetName)
        self.pickleData(self.allStates,allStatesFile)
        self.pickleData(self.allAttrs,allAttrsFile) 
                   
    def getDefaultDir(self):
        defaultDir = '/root'
        if os.path.exists('/data'):
            defaultDir = '/data' 
        return defaultDir
        
    def registerWidget(self):
        if self.widgetDir:
            startDir=self.widgetDir
        else:
            startDir=self.defaultDir
        registerWidgetDir=QtWidgets.QFileDialog.getExistingDirectory(self, caption="Choose widget to register", directory=startDir)
        if registerWidgetDir:
            widgetName=os.path.split(registerWidgetDir)[-1]
            self.register(registerWidgetDir,widgetName)
            
    def pickleData(self,data,filename,jsonFlag=False):
        if jsonFlag:
            myJdata=jsonpickle.encode(data)
            with open(filename,"w") as f:
                f.write(myJdata)
            f.close()
        else:
            with open(filename, 'wb') as f:
                pickle.dump(data, f, protocol=pickle.HIGHEST_PROTOCOL)
            f.close()

        
    def unPickleData(self,filename,jsonFlag=False):
        if jsonFlag:
            with open(filename,'r') as f:
                data=jsonpickle.decode(f)
            f.close()
        else:
            with open(filename, 'rb') as f:
                data=pickle.load(f)
            f.close()
        return data

    def __init__(self,widgetID=None):
        super().__init__()
        css = '''
        QPushButton {background-color: #1588c5; color: white; height: 20px; border: 1px solid black; border-radius: 2px;}
        QPushButton:hover {background-color: #1555f5; }
        QPushButton:hover:pressed { background-color: #1588c5; color: black; border-style: inset; border: 1px solid white} 
        QPushButton:disabled { background-color: lightGray; border: 1px solid gray; } 
        '''  
        
        #self.setStyleSheet(css)
        self.browseIcon=QtGui.QIcon('/biodepot/Bwb_core/icons/bluefile.png')
        self.addIcon=QtGui.QIcon('/biodepot/Bwb_core/icons/add.png')
        self.removeIcon=QtGui.QIcon('/biodepot/Bwb_core/icons/remove.png')
        self.submitIcon=QtGui.QIcon('/biodepot/Bwb_core/icons/submit.png')
        self.reloadIcon=QtGui.QIcon('/biodepot/Bwb_core/icons/reload.png')
        self.outputWidget=""
        self.defaultDir=self.getDefaultDir()
        self.widgetDir=None
        self.widgetName=None
        self.widgetDir=None
        self.isDrawn=False
        self.containerID=None
        #need save mode right at beginning
        self.saveMode=QtGui.QComboBox()
        self.saveMode.addItem('Overwrite')
        self.saveMode.addItem('Merge')
        self.saveMode.addItem('Data')
        self.saveModeIndex=0
        self.saveMode.setCurrentIndex=self.saveModeIndex 
        self.saveMode.currentIndexChanged.connect(lambda: self.onSaveModeChange(self.saveMode))
        sys.stderr.write('widgetID is {}\n'.format(widgetID))
        if widgetID:
            #widgetID is given when this is called from menu and not the initial widget builder
            if widgetID == 'New':
                tmp = tempfile.mkdtemp()
                self.loadWidget()
                self.setWindowTitle(self.widgetName+':Definition')
            else:
                widgetSplit=widgetID.split('.')
                widgetSplit[-1]=widgetSplit[-1][2:]
                self.widgetName=widgetSplit[-1]
                self.setWindowTitle(self.widgetName+':Definition')
                self.widgetDir='/widgets/{}'.format(self.widgetName)
                sys.stderr.write('widgetDir is {} widgetName is {}\n'.format(self.widgetDir,self.widgetName))
                self.loadWidget(loadWidgetDir=self.widgetDir,loadNameCheck=False)
                if 'saveModeIndex' in self.allAttrs:
                    self.saveModeIndex=self.allAttrs['saveModeIndex']
                self.saveMode.setCurrentIndex=self.saveModeIndex 
        else:
            flags = QtGui.QMessageBox.Yes 
            flags |= QtGui.QMessageBox.No
            loadFlag= QtGui.QMessageBox.question(self, "Create new widget","Use existing file as template?",flags)
            if loadFlag == QtGui.QMessageBox.Yes:
                self.loadWidget()
            self.saveWidget()
            if not self.widgetName or not self.widgetDir:
                self.widgetName="generic{}".format(os.getpid())
                self.widgetDir='/tmp/{}'.format(self.widgetName)
            self.setWindowTitle(self.widgetName+':Definition')
            if not self.isDrawn:
                self.startWidget()
    
    def syncNames(self,oldName):
        for myType in ('attrs','states','json','py'):
            oldFile='{}/{}.{}'.format(self.widgetDir,oldName,myType)
            if os.path_exists(oldFile):
                os.system('mv {} []/{}.{}'.format(oldFile,self.widgetName,myType))
        
    def clearLayout(self,layout):
        while layout.count():
            child = layout.takeAt(0)
            if child.widget():
                child.widget().deleteLater()
    
    def startWidget(self):
        self.isDrawn=True
        self.setWindowTitle(self.widgetName+':Definition')
        for attr in ('name','description','category','docker_image_name','docker_image_tag',
                     'priority','icon','inputs','outputs','volumes','parameters','command','autoMap','buildCommand','saveModeIndex'):
            if attr in self.allAttrs:
                self.data[attr]=self.allAttrs[attr]
            else:
                self.data[attr]=None
                self.allAttrs[attr]=None
        self.tabs = tabbedWindow()
        self.controlArea.layout().addWidget(self.tabs)
        self.setStyleSheet(":disabled { color: #282828}")
        
        
        
        #self.clearLayout(self.controlArea.layout())
        #self.controlArea.layout().addWidget(self.scroll_area)
        
        
        #controlBox = gui.vBox(self.generalBox)
        #self.drawExec(box=controlBox)ne)
        #requiredBox = gui.widgetBox(self.generalBox, "Widget entries")
        #draw Ledits for the frequired elements
        leditGeneralLayout=self.tabs.add('General')
        nameLabel=QtGui.QLabel('Name: '+self.widgetName)
        leditGeneralLayout.addWidget(nameLabel,leditGeneralLayout.nextRow,0)
        leditGeneralLayout.nextRow=leditGeneralLayout.nextRow+1
        for pname in ['description','category','docker_image_name','docker_image_tag']:
            self.drawLedit(pname,layout=leditGeneralLayout)
        self.drawLedit('priority',layout=leditGeneralLayout)
        #file entry for icon
        self.drawLedit('icon',layout=leditGeneralLayout,addBrowseButton=True, fileType=None)

        #listwidgets for inputs and outputs 
        #define widgetItems for the different widgetLists
        #top level widgets are drawXXX - these can have multiple substituents
        #lower level widgets are makeXXX - these can also have multiple substituents
        
        self.drawIListWidget('inputs',layout=self.tabs.add('Inputs'))            
        self.drawOListWidget('outputs',layout=self.tabs.add('Outputs'))
        self.drawVolumeListWidget('volumes',layout=self.tabs.add('Volumes'))
        self.drawParamsListWidget('parameters',layout=self.tabs.add('Parameters'))
        self.drawCommand('command',layout=self.tabs.add('Command'))
        self.drawDocker('buildCommand',layout=self.tabs.add('Docker'))
        self.drawExec(self.controlArea.layout())
    
    def updateWidget(self):
        self.clearLayout(self.controlArea.layout())
        self.startWidget()
        
                          
    def updateCheckBox(self,checkBox,widget=None):
        if(checkBox.isEnabled()):
            widget.setEnabled(checkBox.isChecked())
        
    def addListWidget(self,addBtn,qWidgetList,qWidgetItem):
        if qWidgetList.selectedItems():
            index=qWidgetList.row(qWidgetList.selectedItems()[0])
            qWidgetList.states[index]=qWidgetItem.getState()
            qWidgetList.selectedItems()[0].setText(qWidgetItem.printValue())
        else:    
            item=qWidgetItem.printValue()
            qWidgetList.addItem(item)
            qWidgetList.states.append(qWidgetItem.getState())
            sys.stderr.write('adding state {}\n'.format(qWidgetItem.getState()))
            qWidgetItem.blankState()
        qWidgetList.updateAllStates()

    def removeListWidget(self,removeBtn,qWidgetList,qWidgetItem):
        if qWidgetList.selectedItems():
            for item in qWidgetList.selectedItems():
                del qWidgetList.states[qWidgetList.row(item)]
                qWidgetList.takeItem(qWidgetList.row(item))
        if not qWidgetList.count():
            removeBtn.setEnabled(False)
        qWidgetItem.blankState()
        qWidgetList.updateAllStates()
        
    def qwUpdateAllStates(self,pname,qWidgetList,qWidgetItem):
        serialState=[]
        widgetList=[]
        for state in qWidgetList.states:
            serialState.append(state)
            widgetList.append(qWidgetItem.getStateValue(state))
            sys.stderr.write('append state {} value {}\n'.format(state,qWidgetItem.getStateValue(state)))
        if serialState:
            self.allStates[pname]=serialState
            self.data[pname]=OrderedDict(widgetList)
        else:
            self.allStates[pname]=None
            self.data[pname]=qWidgetItem.getStateValue(None)
            
    def qwInitAllStates(self,pname,qWidgetList,qWidgetItem):
        widgetList=[]
        if pname in self.allStates and self.allStates[pname]:
            for serialState in self.allStates[pname]:
                sys.stderr.write('init pname {} state {}\n'.format(pname, serialState))
                qWidgetItem.blankState()
                qWidgetItem.setState(serialState)
                item=qWidgetItem.printValue()
                qWidgetList.addItem(item)
                qWidgetList.states.append(qWidgetItem.getState())
                widgetList.append(qWidgetItem.getStateValue(serialState))
        if widgetList:
            self.data[pname]=OrderedDict(widgetList)
        else:
            self.data[pname]=None
        qWidgetItem.blankState()
    
    def onListWidgetSelect(self,qWidgetList,addBtn,removeBtn,qWidgetItem):
        if qWidgetList.selectedItems():
            if len(qWidgetList.selectedItems()) > 1:
                addBtn.setEnabled(False)
            else:
                item=qWidgetList.selectedItems()[0]
                myState=qWidgetList.states[qWidgetList.row(item)]
                sys.stderr.write('loading {}\n'.format(myState))
                qWidgetItem.setState(myState)
                addBtn.setEnabled(True)
        else:
            removeBtn.setEnabled(False)
            addBtn.setEnabled(True)
            
    def onItemMoved(self,oldRow,newRow,qWidgetList):
        sys.stderr.write("oldRow is {} newRow is {}\n".format(oldRow,newRow))
        temp=qWidgetList.states[oldRow]
        qWidgetList.states[newRow+1:oldRow+1]=qWidgetList.states[newRow:oldRow]
        qWidgetList.states[newRow]=temp
        qWidgetList.updateAllStates()

    def drawCommand(self,pname,layout=None):
        labelTextBox=self.makeTextBox(pname,label='Enter command:')
        self.initAllStates(pname,labelTextBox)
        labelTextBox.textBox.textChanged.connect(lambda: self.updateAllStates(pname,labelTextBox,labelTextBox.getState()))
        #layout.addWidget(labelTextBox.label,layout.nextRow,0)
        layout.addWidget(labelTextBox.textBox,layout.nextRow,1,1,4)
        layout.nextRow = layout.nextRow + 1
        
    def drawDocker(self,pname,layout=None):
        css = '''
        QPushButton {background-color: #1588c5; color: white; height: 20px; border: 1px solid black; border-radius: 2px;}
        QPushButton:hover {background-color: #1555f5; }
        QPushButton:hover:pressed { background-color: #1588c5; color: black; border-style: inset; border: 1px solid white} 
        QPushButton:disabled { background-color: lightGray; border: 1px solid gray; }
        '''
        #self.drawLedit('Add file to Dockerfiles',layout=layout,addBrowseButton=True, fileType=None)
        #self.drawLedit('Add directory to Dockerfiles',layout=layout,addBrowseButton=True, fileType='Directory')
        addDateCb=self.makeCheckBox ('addBuildDate','Add date to docker tag',default=True,persist=True,track=True)
        containerIDLabel=QtGui.QLabel('Container ID: {}'.format(self.containerID))
        imageBuilderLabel=QtGui.QLabel('Launch Image Builder')
        imageBuilderBtn = gui.button(None, self, "Launch", callback=self.startImageBuilder)
        imageBuilderBtn.setStyleSheet(css)
        imageBuilderBtn.setFixedSize(80,20)
        updateLabel=QtGui.QLabel('Import Dockerfiles directory')
        updateBtn = gui.button(None, self, "Import", callback=self.updateDockerfiles)
        updateBtn.setStyleSheet(css)
        updateBtn.setFixedSize(80,20)        
        buildCommandBox=self.makeTextBox(pname,label='Docker build command:')
        self.initAllStates(pname,buildCommandBox)
        buildCommandBox.textBox.textChanged.connect(lambda: self.updateAllStates(pname,buildCommandBox,buildCommandBox.getState()))
        layout.addWidget(updateLabel,layout.nextRow,0)
        layout.addWidget(updateBtn,layout.nextRow,1)
        layout.addWidget(imageBuilderLabel,layout.nextRow+1,0)
        layout.addWidget(imageBuilderBtn,layout.nextRow+1,1)
        layout.addWidget(addDateCb,layout.nextRow+2,0)
        layout.addWidget(containerIDLabel,layout.nextRow+3,0)
        layout.addWidget(addDateCb,layout.nextRow+4,0)
        layout.addWidget(buildCommandBox.label,layout.nextRow+5,0)
        layout.addWidget(buildCommandBox.textBox,layout.nextRow+5,1,1,4)
    def startImageBuilder(self):
        widget=OWImageBuilder.OWImageBuilder()
        widget.showNormal()
        widget.raise_()
        widget.activateWindow()
        
    def updateDockerfiles(self):
        pass
    def makeTextBox(self,attr,label):
        box=QHBoxLayout()
        textLabel=None
        if(label):
            textLabel=QtGui.QLabel(label)
            textLabel.setAlignment(Qt.AlignTop)
        textBox=QtGui.QPlainTextEdit()
        textBox.setStyleSheet(":disabled { color: #282828}")
        setattr(box,'textBox',textBox)
        setattr(box,'label',textLabel)
        setattr(box,'getValue',textBox.toPlainText())
        setattr(box,'getStateValue',lambda state: self.getTextBoxStateValue(box,state))
        setattr(box,'getState',lambda: self.getTextBoxState(box))
        setattr(box,'setState',lambda state : self.setTextBoxState(box,state) )
        return box

    def getTextBoxState(self,widget):
        return [widget.textBox.isEnabled(), widget.textBox.toPlainText()]
        
    def getTextBoxStateValue(self,widget,state):
        if state and state[0] and state[1]:
            return state[1]
        return None
                    
    def setTextBoxState(self,widget,state):
        if state is None:
            widget.label.setEnabled(True)
            widget.textBox.setEnabled(True)
            widget.textBox.clear()
        else:
            widget.textBox.setEnabled(state[0])
            widget.textBox.setPlainText(state[1])

    #workhorse widget
    def makeLedit(self,leditAttr,text=None,label=None,addCheckBox=False,initialValue=None):
        leditLabel=None
        checkBox=None
        if(label):
            leditLabel=QtGui.QLabel(label)
        ledit = QtGui.QLineEdit(self)
        ledit.setClearButtonEnabled(True)
        ledit.setPlaceholderText(text)
        ledit.setStyleSheet(":disabled { color: #282828}")
        if not initialValue:
            ledit.clear()
        if addCheckBox:
            checkBox=QtGui.QCheckBox(None,self)
            checkBox.stateChanged.connect(lambda : self.updateCheckBox(checkBox,ledit))
            ledit.setEnabled(checkBox.isChecked())            
        box=QHBoxLayout()
        if checkBox:
            box.addWidget(checkBox)
        if leditLabel:
            box.addWidget(leditLabel)
            box.addWidget(ledit)
        setattr(box,'ledit',ledit)
        setattr(box,'label',leditLabel)
        setattr(box,'checkBox',checkBox)
        setattr(box,'getValue',lambda : self.getLeditValue(ledit))
        setattr(box,'getStateValue',lambda state : self.getLeditStateValue(checkBox,leditLabel,ledit,state))
        setattr(box,'getState',lambda : self.getLeditState(checkBox,leditLabel,ledit))
        setattr(box,'setState',lambda state : self.setLeditState(checkBox,leditLabel,ledit,state))
        return box
        
    def getLeditValue(self,ledit):
        if ledit.isEnabled():
            return ledit.text()
        return None
    
    def getLeditStateValue(self,checkBox,label,ledit,state):
        sys.stderr.write('getledit state vector is {}\n'.format(state)) 
        if(state is None):
            return None
        if checkBox:
            if not state[0][1] :
                return None
        return state[2][1]
            
    def getLeditState(self,checkBox,label,ledit):
        #can't serialize qt objects
        #so we just concatenate the data as a list instead of a dict 
        state=[]
        if checkBox:
            state.append([checkBox.isEnabled(),checkBox.isChecked()])
        else:
            state.append([None,None])
        if label:
            state.append([label.isEnabled(),label.text()])
        else:
            state.append([None,None])
        if ledit:
            state.append([ledit.isEnabled(),ledit.text()])
        else:
            state.append([None,None])
        return state
        
    def setLeditState(self,checkBox,label,ledit,state):
        if state is None:
            #initialize it
            if checkBox:
                checkBox.setEnabled(True)
                checkBox.setChecked(False)
                ledit.clear()
                ledit.setEnabled(False)
                label.setEnabled(True)
            else:
                ledit.setEnabled(True)
                ledit.clear()
                label.setEnabled(True)
            return
        if checkBox is not None:
            checkBox.setEnabled(state[0][0])
            checkBox.setChecked(state[0][1])
        
        if label is not None:
            sys.stderr.write('labelEnable is {} labelState is {}\n'.format(state[1][0],state[1][1]))
            label.setEnabled(state[1][0])
            label.setText(state[1][1])
                  
        if ledit is not None:
            ledit.setEnabled(state[2][0])
            if state[1][1]:
                ledit.setText(state[2][1])
            else:
                ledit.clear()
                
    def makeComboBox (self,pname,label, elements):
        comboBoxLabel=QtGui.QLabel(label)
        comboBox=QtGui.QComboBox()
        comboBox.addItems(elements)
        comboBox.currentIndex=0
        box=QHBoxLayout()
        box.addWidget(comboBoxLabel)
        box.addWidget(comboBox)
        setattr(box,'getValue',lambda : self.getComboValue(comboBox))
        setattr(box,'getStateValue',lambda state: self.getComboStateValue(comboBox,state))
        setattr(box,'getState',lambda : self.getComboState(comboBox,comboBoxLabel))
        setattr(box,'setState',lambda state : self.setComboState(comboBox,comboBoxLabel,state))
        return box
        
    def getComboValue(self,comboBox):
        if comboBox.isEnabled():
            return comboBox.currentText()  
        return None 
        
    def getComboStateValue(self,comboBox,state):
        if not state or not state[0]:
            return None
        elif state[0][0]: #check if enabled
            index=state[0][1]
            if index >= 0:
                allItems=state[0][2]
                return allItems[index]
            else:
                return None
        return None    
                 
    def getComboState(self,comboBox,label):
        state=[]
        if comboBox:
            allItems=[comboBox.itemText(i) for i in range (comboBox.count())]
            enabled=comboBox.isEnabled()
            index=comboBox.findText(comboBox.currentText())
            state.append([enabled,index,allItems])
        else:
            state.append(None)
        if label:
            state.append([label.isEnabled(),label.text()])
        else:
            state.append(None)
        return state 
    
    def setComboState(self,comboBox,label,state):
        if state is None:
            #intialize
            comboBox.setEnabled(True)
            comboBox.setCurrentIndex(0)
            return None
        if state[0] is None:
            comboBox=None
        else:
           comboBox.setEnabled(state[0][0])
           comboBox.clear()
           comboBox.addItems(state[0][2])
           comboBox.setCurrentIndex(state[0][1])
        if state[1] is None:
            label=None
        else:
            label.setEnabled(state[1][0])
            label.setText(state[1][1])
        
    def makeCheckBox (self, attr,label,default=False,persist=False,track=False):
        #checkbox for binary options
        #not used as part of other elements
        if attr not in self.data:
            self.data[attr]=default
        checkBox=QtGui.QCheckBox(label,self)
        setattr(checkBox,'getValue',lambda : checkBox.isEnabled() and checkBox.isChecked())
        setattr(checkBox,'getStateValue',lambda state: self.getCheckBoxStateValue(checkBox,state,default))
        setattr(checkBox,'getState',lambda : self.getCheckBoxState(checkBox))
        setattr(checkBox,'setState',lambda state : self.setCheckBoxState(checkBox,state,default))
        if persist:
            self.initAllStates(attr,checkBox)
            checkBox.stateChanged.connect(lambda : self.updateAllStates(attr,checkBox,checkBox.getState()))
        if track:
            if not persist:
                checkBox.stateChanged.connect(lambda : self.updateAllStates(attr,checkBox,checkBox.getState()))
        return checkBox
            
    def getCheckBoxStateValue(self,checkBox,state,default):
        if not state:
            return default
        if state[0] and state[1]:
            return True
        return False
        
    def getCheckBoxState(self,checkBox):
        state=[]
        if checkBox:
            state=[checkBox.isEnabled() ,checkBox.isChecked()]
        else:
            state=None
        return state
        
    def setCheckBoxState(self,checkBox,state,default=False):
        if state is None:
            #initialize
            checkBox.setEnabled(True)
            checkBox.setChecked(default)
            return state
        if state[0] is None:
            checkBox=None
        else:
            checkBox.setEnabled(state[0])
            checkBox.setChecked(state[1])

    def makeListWidget(self,pname,boxEdit):
        #setup boxEdit
        #logic is handled by add remove buttons
        if pname not in self.data:
            self.data[pname]=None; 
        boxEdit=DragAndDropList(self)
        boxEdit.setDragDropMode(QtWidgets.QAbstractItemView.InternalMove)
        boxEdit.setSelectionMode(QtWidgets.QAbstractItemView.ExtendedSelection)
        boxEdit.setStyleSheet(":disabled { color: #282828}")
        boxEdit.setMinimumHeight(60)
        return boxEdit
        
    #multiple widget elements
    def makeListWidgetUnit (self, pname,layout=None, lineWidgets=None, otherWidgets=None, elementsPerLine=4):
        #widgets is ordered dict of widgets managed by listWidget
        boxEdit=self.makeListWidget(pname,None)
        setattr(boxEdit,'states',[])
        lineItem=WidgetItem(OrderedDict(lineWidgets),{})
        setattr(boxEdit,'updateAllStates',lambda: self.qwUpdateAllStates(pname,boxEdit,lineItem))
        self.qwInitAllStates(pname,boxEdit,lineItem)
        setattr(boxEdit,'oldRow',boxEdit.currentRow())
        #buttons
        addBtn=gui.button(None, self, "", callback=lambda: self.addListWidget(addBtn,boxEdit,lineItem), autoDefault=False)
        removeBtn=gui.button(None, self, "", callback=lambda: self.removeListWidget(removeBtn,boxEdit,lineItem), autoDefault=False)
        removeBtn.setEnabled(bool(boxEdit.selectedItems()))
        boxEdit.itemSelectionChanged.connect(lambda: removeBtn.setEnabled(bool(boxEdit.selectedItems())))
        boxEdit.itemSelectionChanged.connect(lambda: self.onListWidgetSelect(boxEdit,addBtn,removeBtn,lineItem))
        boxEdit.itemMoved.connect(lambda oldRow,newRow : self.onItemMoved(oldRow,newRow,boxEdit))
        addBtn.setIcon(self.addIcon)
        removeBtn.setIcon(self.removeIcon)
        #set button styles
        buttonStyle='background: None; border: None ; border-radius: 0;'
        addBtn.setStyleSheet(buttonStyle)
        removeBtn.setStyleSheet(buttonStyle)
        
        #layout
        #basic layout is other widgets on top with filesBox below and the lineItem underneath
        
        filesBoxLeditLayout=QtGui.QVBoxLayout()
        #add to the main parameters box
        myBox=gui.vBox(None)
        myBox.layout().addLayout(filesBoxLeditLayout)
        #label=QtGui.QLabel(pname+':')
        #label.setAlignment(Qt.AlignTop)
        #layout.addWidget(label,layout.nextRow,0)
        #otherwidgets here
        if otherWidgets:
            for widget in otherWidgets:
                layout.addWidget(widget,layout.nextRow,1,1,2)
                layout.nextRow = layout.nextRow + 1
        #box with list here
        layout.addWidget(myBox,layout.nextRow,1,1,2)
        layout.nextRow = layout.nextRow + 1
        #input with line layout here     
        lineLayout=QtGui.QGridLayout()
        col=0
        row=1
        for item in lineWidgets:
            layout=item[1]
            #The Widgets are actually layouts except for the checkboxes
            if col >= elementsPerLine:
                row+=1
                col=0
            if isinstance(layout,QtGui.QCheckBox):
                lineLayout.addWidget(layout,row,col)
            else:
                lineLayout.addLayout(layout,row,col)
            col+=1
        lineLayout.addWidget(addBtn,row,col)
        lineLayout.addWidget(removeBtn,row,col+1)
        
        #now add the two layouts to the bigBox layout
        filesBoxLeditLayout.addWidget(boxEdit)
        filesBoxLeditLayout.addLayout(lineLayout)
        
    def drawIListWidget (self, pname, layout=None):
        nameBox=self.makeLedit(pname+'nameLedit','Enter name','Name')
        callbackBox=self.makeLedit(pname+'callbackLedit','Enter callback', 'callback',addCheckBox=True)
        comboBox=self.makeComboBox(pname,'Type:',['str','Orange.data.Table'])   
        widgetList=[('name',nameBox),('callback',callbackBox),('type',comboBox)]
        self.makeListWidgetUnit (pname, layout=layout, lineWidgets=widgetList)
        
    def drawOListWidget (self, pname, layout=None):
        nameBox=self.makeLedit(pname+'nameLedit','Enter name','Name')
        defaultBox=self.makeLedit(pname+'defaultLedit','Enter default', 'Default value',addCheckBox=True)
        comboBox=self.makeComboBox(pname,'Type:',['str','Orange.data.Table'])       
        widgetList=[('name',nameBox),('default',defaultBox),('type',comboBox)]
        self.makeListWidgetUnit (pname, layout=layout, lineWidgets=widgetList)


    def drawVolumeListWidget (self,pname,layout=None):
        nameBox=self.makeLedit(pname+'Name','Enter name','Name')
        volumeBox=self.makeLedit(pname+'volumeLedit','Enter volume','Additional volume')
        widgetList=[('name',nameBox),('containerVolume',volumeBox)]
        autoMapCb=self.makeCheckBox ('autoMap','Pass current Bwb volumes to container',default=True,persist=True,track=True)
        self.makeListWidgetUnit (pname, layout=layout, lineWidgets=widgetList,otherWidgets=[autoMapCb])
        
    def drawParamsListWidget (self, pname, layout=None):
        nameBox=self.makeLedit(pname+'nameLedit','Enter name','Name')
        flagBox=self.makeLedit(pname+'flagLedit','Enter flag', 'flag',addCheckBox=True)
        labelBox=self.makeLedit(pname+'labelLedit','Enter label', 'label',addCheckBox=True)
        envBox=self.makeLedit(pname+'envLedit','Enter ENV variable', 'env',addCheckBox=True)
        defaultBox=self.makeLedit(pname+'defaultLedit','Enter default', 'default',addCheckBox=True)
        setattr(self,'optional',False)
        setattr(self,'argument',False)
        optionalCb=self.makeCheckBox ('optional','Optional')
        argumentCb=self.makeCheckBox ('argument','Argument')
        #connect argument checkbox to disabling the flag checkbox
        flagBox.checkBox.stateChanged.connect(lambda : (not flagBox.checkBox.isChecked()) or (argumentCb.setChecked(False)))
        argumentCb.stateChanged.connect(lambda : (not argumentCb.isChecked()) or (flagBox.checkBox.setChecked(False) or flagBox.ledit.clear()))       
        comboBox=self.makeComboBox(pname,'Type:',['str','file','file list','directory','directory list','bool','bool list','text','text list','int','int list','double','double list'])
        widgetList=[('name',nameBox),
                    ('type',comboBox),
                    ('flag',flagBox), 
                    ('argument',argumentCb),
                    ('env',envBox),
                    ('label',labelBox), 
                    ('default',defaultBox),
                    ('optional',optionalCb)
                   ]
        self.makeListWidgetUnit (pname, layout=layout, lineWidgets=widgetList)

    def drawLedit(self,pname,layout=None,addBrowseButton=False,fileType=None, callback=None):
        #make labeledLedit combo layout
        labelLedit=self.makeLedit(pname,'Enter '+ pname, label=pname+':')
        self.initAllStates(pname,labelLedit)
        labelLedit.ledit.textChanged.connect(lambda: self.updateAllStates(pname,labelLedit,labelLedit.getState()))
        layout.addWidget(labelLedit.label,layout.nextRow,0)
        if addBrowseButton:
            layout.addWidget(labelLedit.ledit,layout.nextRow,1,1,1)
            if callback is None:
                callback=self.browseFileDir
            button=gui.button(None, self, "", callback= lambda : callback(pname,ledit=labelLedit.ledit,fileType=fileType),autoDefault=True, width=19, height=19)
            button.setIcon(self.browseIcon)
            layout.addWidget(button,layout.nextRow,2)
        else:
            layout.addWidget(labelLedit.ledit,layout.nextRow,1,1,2)
        layout.nextRow+=1
        return labelLedit
    
    def updateAllStates(self,pname,widget,state):
        self.allStates[pname]=state
        self.data[pname]=widget.getStateValue(state)
        sys.stderr.write('updating all - pname {} state {}\n'.format(pname,state))
        
    def initAllStates(self, attr ,widget):
        if attr in self.allStates:
            widget.setState(self.allStates[attr])
        else:
            widget.setState(None)
            self.allStates[attr]=widget.getState()
        self.data[attr]=widget.getStateValue(self.allStates[attr])
            
    def browseFileDir(self, attr,ledit=None,fileType=None):
        self.defaultDir = self.getDefaultDir()
        if fileType == 'Directory':
            myFileDir=QtWidgets.QFileDialog.getExistingDirectory(self, caption="Locate directory", directory=self.defaultDir)
        else:
            myFileDir=QtWidgets.QFileDialog.getOpenFileName(self, "Locate file", self.defaultDir)[0]
        if myFileDir:
            self.data[pname]=myFileDir
            self.defaultDir=myFileDir
        if ledit:
            ledit.setText(myFileDir)
            
    def  updateCheckBoxLayout(self, cb, layout):
        if cb.isChecked():
            for i in reversed(range(layout.count())):
                widget=layout.itemAt(i).widget()
                if isinstance(widget,QCheckBox):
                    widget.setChecked(False)
                    widget.setEnabled(False)
                elif isinstance(widget,QLabel):
                    continue
                elif isinstance(widget,QtWidgets.QLineEdit):
                    widget.clear()
                    widget.setEnabled(False)
                else:
                    widget.setEnabled(False)
        else:
            for i in reversed(range(layout.count())):
                widget=layout.itemAt(i).widget()
                if isinstance(widget,QCheckBox):
                    widget.setChecked(False)
                    widget.setEnabled(True)
                elif isinstance(widget,QLabel):
                    continue
                elif isinstance(widget,QtWidgets.QLineEdit):
                    widget.setEnabled(False)
                else:
                    widget.setEnabled(True)

