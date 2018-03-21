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


#for selecting multiple directories
class getExistingDirectories(QtWidgets.QFileDialog):
    def __init__(self, *args):
        super(getExistingDirectories, self).__init__(*args)
        self.setOption(self.DontUseNativeDialog, True)
        self.setFileMode(self.Directory)
        self.setOption(self.ShowDirsOnly, True)
        self.findChildren(QtWidgets.QListView)[0].setSelectionMode(QtWidgets.QAbstractItemView.ExtendedSelection)
        self.findChildren(QtWidgets.QTreeView)[0].setSelectionMode(QtWidgets.QAbstractItemView.ExtendedSelection)
#dummy class
class ContainerPaths():
    pass

class BwbGuiElements():
    #keep track of List of Gui elements
    #add support for initialization callbacks
    #keep track of values
    
    #need to check for tuple elements from Orange gui
    def disableElement(self,element):
        if isinstance(element,tuple):
            for g in element:
                g.setDisabled(True)
        else:
            element.setDisabled(True)
            
    def enableElement(self,element):
        if isinstance(element,tuple):
            for g in element:
                g.setEnabled(True)
        else:
            element.setEnabled(True)
                   
    def __init__(self,required=None,active=None):
        self._enableCallbacks={}
        self._updateCallbacks={}
        self._dict={}
        self.required=required
        self.active={}

    def add(self,attr,guiElement,enableCallback=None,updateCallback=None):
        if attr not in self._dict:
            self._dict[attr]=[]
        self._dict[attr].append(guiElement)
        if enableCallback is not None:
            self._enableCallbacks[attr]=enableCallback
        if updateCallback is not None:
            self._updateCallbacks[attr]=updateCallback
            
    def addList(self,attr,guiElements,enableCallback=None,updateCallback=None):
        if attr not in self._dict:
            self._dict[attr]=guiElements
        else:
            self._dict[attr].extend(guiElements)
        if enableCallback is not None:
            self._enableCallbacks[attr]=enableCallback
        if updateCallback is not None:
            self._updateCallbacks[attr]=updateCallback
                
    def disable(self,attr,ignoreCheckbox=False): #gray out to disable
        if attr in self._updateCallbacks:
            self._updateCallbacks[attr]()
        if attr in self._dict:
            if ignoreCheckbox:
                for g in self._dict[attr]:
                    if not isinstance(g, QtWidgets.QCheckBox):
                        sys.stderr.write('ignore cb, {} disabling {}\n'.format(attr,g))
                        self.disableElement(g)
            else:
                for g in self._dict[attr]:
                    sys.stderr.write('{} disabling {}\n'.format(attr,g)),
                    self.disableElement(g)
            return True
        return False

    def enable(self,attr,value):
        clearLedit=False
        if value is None:
            clearLedit=True
        if attr in self._updateCallbacks:
            self._updateCallbacks[attr]() 
        if attr in self._dict:
            if attr in self._enableCallbacks:
                self._enableCallbacks[attr](value,clearLedit)
            else:
                for g in self._dict[attr]:
                    self.enableElement(g)
        return True
        
    def disableAll(self):
        for attr in self._dict.keys():
            self.disable(attr)

            
    def reenableAll(self,OWself):
        for attr in self._dict.keys():
            if not OWself.inputConnections.isConnected(attr) :
                self.enable(attr,getattr(OWself,attr))

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
            if connectionId is None:
                del self._dict[slot]
            else:
                try:
                    self._dict[slot].remove(connectionId)
                except ValueError:
                    pass
            #need to update the volumes

    def isConnected (self, slot):
        if self._dict is None:
            return False
        if (slot in self._dict):
            return True
        return False

    def isSet (self, slot): #same as is connected but checks for non null value
        if self._dict is None:
            return False
        if (slot in self._dict) and (self._dict[slot] is not None):
            return True
        return False

class OWBwBWidget(widget.OWWidget):
    dockerClient = DockerClient('unix:///var/run/docker.sock', 'local')
    defaultFileIcon=QtGui.QIcon('/biodepot/orangebiodepot/icons/bluefile.png')
    browseIcon=QtGui.QIcon('/biodepot/orangebiodepot/icons/bluefile.png')
    addIcon=QtGui.QIcon('/biodepot/orangebiodepot/icons/add.png')
    removeIcon=QtGui.QIcon('/biodepot/orangebiodepot/icons/remove.png')
    submitIcon=QtGui.QIcon('/biodepot/orangebiodepot/icons/submit.png')
    reloadIcon=QtGui.QIcon('/biodepot/orangebiodepot/icons/reload.png')
    
#Initialization
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
        self.filesBoxLayout=QtGui.QVBoxLayout()
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
        #keep track of gui elements associated with each attribute
        self.bgui=BwbGuiElements()
        #keep track of all the connections - n
        self.inputConnections=ConnectionDict(self.inputConnectionsStore)

    def initVolumes(self):
        #initializes container volumes
        for mapping in self.data['volumeMappings']:
            bwbVolAttr=mapping['attr']
            if not hasattr(self, bwbVolAttr) :
                setattr(self,bwbVolAttr,None)
            bwbVol=getattr(self,bwbVolAttr)
            setattr(self.con, bwbVolAttr, mapping['conVolume'])

#Drawing the GUI
    def drawGUI(self):
        self.setStyleSheet(":disabled { color: #282828}")
        self.scroll_area = QScrollArea(
            verticalScrollBarPolicy=Qt.ScrollBarAlwaysOn
        )
        self.bigBox=gui.widgetBox(self.controlArea)
        self.scroll_area.setWidget(self.bigBox)
        self.scroll_area.setWidgetResizable(True)
        self.controlArea.layout().addWidget(self.scroll_area)
        consoleBox = gui.widgetBox(self.bigBox,None)
        self.infoLabel = gui.widgetLabel(consoleBox, None)
        self.infoLabel.setWordWrap(True)
        if 'requiredParameters' in self.data and self.data['requiredParameters']:
            self.requiredBox = gui.widgetBox(self.bigBox, "Required parameters")
            self.drawRequiredElements()
            
        if self.findOptionalElements():
            self.optionalBox = gui.widgetBox(self.bigBox, "Optional parameters")
            self.drawOptionalElements()
            
        #disable connected elements
        for i in self.inputs:
            attr=i.name
            if self.inputConnections.isConnected(attr):
               self.bgui.disable(attr)
            
        #check if the requirements to be run are met
        self.checkTrigger()
        controlBox = QtGui.QVBoxLayout()
        self.drawExec(box=self.controlArea.layout())

        
    def drawRequiredElements(self):
        for pname in self.data['requiredParameters']:
            if not ('parameters' in self.data) or not (pname in self.data['parameters']):
                continue
            pvalue=self.data['parameters'][pname]
            if not hasattr(self,pname):
                setattr(self,pname,None)
            if (getattr(self,pname) is None) and ('default' in pvalue):
                setattr(self,pname,pvalue['default'])

        for pname in self.data['requiredParameters']:
            if not ('parameters' in self.data) or not ( pname in self.data['parameters']):
                continue
            pvalue=self.data['parameters'][pname]
            if ('gui' in pvalue and pvalue['gui'] != 'file' and pvalue['gui'] != 'directory') or ( pvalue['type'] != 'file' and pvalue['type'] != 'directory'):
                continue
            self.drawFileDirElements(pname, pvalue, box=self.requiredBox,layout=self.fileDirRequiredLayout)
            
        for pname in self.data['requiredParameters']:
            if not ('parameters' in self.data) or not ( pname in self.data['parameters']):
                continue
            pvalue=self.data['parameters'][pname]
            if ('gui' in pvalue and pvalue['gui'][-4:] != 'list') or ( pvalue['type'][-4:] != 'list' ):
                continue
            self.drawFilesBox(pname, pvalue, box=self.requiredBox,layout=self.fileDirRequiredLayout)

        for pname in self.data['requiredParameters']:
            if not ('parameters' in self.data) or not( pname in self.data['parameters']):
                continue
            pvalue=self.data['parameters'][pname]
            if ('gui' in pvalue and pvalue['gui'] != 'Ledit') or (pvalue['type'] != 'double') and (pvalue['type'] != 'text'):
                continue
            self.drawLedit(pname,pvalue,self.requiredBox,layout=self.leditRequiredLayout)

        for pname in self.data['requiredParameters']:
            if not ('parameters' in self.data) or not ( pname in self.data['parameters']):
                continue
            pvalue=self.data['parameters'][pname]
            if ('gui' in pvalue and pvalue['gui'] != 'Spin') or (pvalue['type'] != 'int'):
                continue
            self.drawSpin(pname,pvalue,self.requiredBox)

        for pname in self.data['requiredParameters']:
            if not ('parameters' in self.data) or not ( pname in self.data['parameters']):
                continue
            pvalue=self.data['parameters'][pname]
            if ('gui' in pvalue and pvalue['gui'] != 'bool') or (pvalue['type'] != 'bool'):
                continue
            self.drawCheckbox(pname,pvalue,self.requiredBox)
    
    def findOptionalElements(self):
       #checks that there are optional elements
        if not 'parameters' in self.data:
            return False
        for pname in self.data['parameters']:
            if pname not in self.data['requiredParameters']:
                return True
        return False
 
    def drawOptionalElements(self):
        for pname in self.data['parameters']:
            if pname not in self.data['requiredParameters']:
                pvalue=self.data['parameters'][pname]
                if not hasattr(self,pname):
                    setattr(self,pname,None)
                if (getattr(self,pname) is None) and ('default' in pvalue):
                    setattr(self,pname,pvalue['default'])
                if not (pname in self.optionsChecked):
                    self.optionsChecked[pname]=False

        for pname in self.data['parameters']:
            if not ('parameters' in self.data) or not ( pname in self.data['parameters']):
                continue
            if pname not in self.data['requiredParameters']:
                pvalue=self.data['parameters'][pname]
                if ('gui' in pvalue and pvalue['gui'] != 'file' and pvalue['gui'] != 'directory') or (pvalue['type'] != 'file') and (pvalue['type'] != 'directory'):
                    continue
                self.drawFileDirElements(pname, pvalue, box=self.optionalBox,layout=self.fileDirOptionalLayout,addCheckbox=True)
                
        for pname in self.data['parameters']:
            if not ('parameters' in self.data) or not ( pname in self.data['parameters']):
                continue
            if pname not in self.data['requiredParameters']:
                pvalue=self.data['parameters'][pname]
                if ('gui' in pvalue and pvalue['gui'][-4:] != 'list') or ( pvalue['type'][-4:] != 'list' ):                    
                    continue
                self.drawFilesBox (pname, pvalue, box=self.optionalBox,layout=self.fileDirOptionalLayout,addCheckbox=True)
                
        for pname in self.data['parameters']:
            if not ('parameters' in self.data) or not ( pname in self.data['parameters']):
                continue
            if pname not in self.data['requiredParameters']:
                pvalue=self.data['parameters'][pname]
                if ('gui' in pvalue and pvalue['gui'] != 'Ledit') or (pvalue['type'] != 'double') and (pvalue['type'] != 'text'):
                    continue
                self.drawLedit(pname,pvalue,self.optionalBox,layout=self.leditOptionalLayout,addCheckbox=True)

        for pname in self.data['parameters']:
            if not ('parameters' in self.data) or not ( pname in self.data['parameters']):
                continue
            if pname not in self.data['requiredParameters']:
                pvalue=self.data['parameters'][pname]
                if ('gui' in pvalue and pvalue['gui'] != 'Spin') or (pvalue['type'] != 'int'):
                    continue
                self.drawSpin(pname,pvalue,self.optionalBox,addCheckbox=True)

        for pname in self.data['parameters']:
            if not ('parameters' in self.data) or not ( pname in self.data['parameters']):
                continue
            if pname not in self.data['requiredParameters']:
                pvalue=self.data['parameters'][pname]
                if ('gui' in pvalue and pvalue['gui'] != 'bool') or (pvalue['type'] != 'bool'):
                    continue
                sys.stderr.write('drawOpt pname {} pvalue {}\n'.format(pname, pvalue))
                self.drawCheckbox(pname,pvalue,self.optionalBox)

    def drawCheckbox(self,pname,pvalue, box=None):
        #for booleans - their value is the same as the checkbox state
        sys.stderr.write('drawCB pname {} pvalue {} label {}\n'.format(pname, pvalue,pvalue['label']))
        if not hasattr(self,pname):
            if 'default' in pvalue:
                setattr(self,pname,pvalue['default'])
            else:
                setattr(self,pname,False)
        cb=gui.checkBox(box, self, pname, pvalue['label'])
        checkAttr=pname+'Checked'
        setattr(self,checkAttr,getattr(self,pname))
        #check if inactive
        self.bgui.add(pname,cb)
            
    def drawSpin(self,pname,pvalue, box=None,addCheckbox=False):
        #for drawSpin - we use the origin version which already has a checkbox connected
        #TODO could change this to the same way we handle ledits with separate cbo
        checkAttr=None
        if addCheckbox:
            checkAttr=pname+'Checked'
            setattr(self,checkAttr,self.optionsChecked[pname])
        mySpin=gui.spin(box, self, pname, minv=1, maxv=128, label=pvalue['label'], checked=checkAttr, checkCallback=lambda : self.updateSpinCheckbox(pname))
        if getattr(self,pname) is None:
            mySpin.clear()
        self.bgui.add(pname,mySpin,enableCallback=lambda value,clearLedit  : self.enableSpin(value,clearLedit,mySpin))
        
    def drawLedit(self,pname,pvalue,box=None,layout=None,addCheckbox=False):
        checkAttr=None
        checkbox=None
        ledit=gui.lineEdit(None, self, pname,disabled=addCheckbox)
        if addCheckbox:
            checkAttr=pname+'Checked'
            setattr(self,checkAttr,self.optionsChecked[pname])
            checkbox=gui.checkBox(None, self,checkAttr,label=None)
            checkbox.stateChanged.connect(lambda : ledit.setEnabled(checkbox.isChecked()))
            checkbox.stateChanged.connect(lambda : self.updateCheckbox(pname,checkbox.isChecked(),ledit.text()))
            self.bgui.add(pname,checkbox)
        self.bwbLedit(box,checkbox,ledit,layout=layout, label=pvalue['label'])
        #check if the value is none - then we clear it
        if getattr(self,pname) is None:
            ledit.clear()
        self.bgui.add(pname,ledit,enableCallback=lambda value,clearLedit: self.enableLedit(value,clearLedit,checkbox,ledit))
        if addCheckbox:
            self.updateCheckbox(pname,checkbox.isChecked(),ledit.text())
    
    def enableLedit(self,value,clearLedit,cb,ledit):
        if cb:
            cb.setEnabled(True)
            if cb.isChecked():
                ledit.setEnabled(True)
            else:
                ledit.setEnabled(False)
        else:
            ledit.setEnabled(True)
        if value is None or clearLedit:
            ledit.clear()
        
    def drawFileDirElements(self,pname, pvalue, box=None, layout=None, addCheckbox=False):
        checkbox=None
        ledit=gui.lineEdit(None, self, pname,disabled=addCheckbox)
        #note that using lambda does not work in this function - truncates string variable so partial used instead
        button=gui.button(None, self, "", callback= partial(self.browseFileDir, attr= pname ,filetype=pvalue['type']),
                          autoDefault=True, width=19, height=19,disabled=addCheckbox)
        if getattr(self,pname) is None:
            ledit.clear()
        if addCheckbox:
            checkAttr=pname+'Checked'
            setattr(self,checkAttr,self.optionsChecked[pname])
            checkbox=gui.checkBox(None, self,checkAttr,label=None)
            sys.stderr.write('updating filedir {}\n'.format(pname))
            checkbox.stateChanged.connect(lambda : self.updateCheckbox(pname,checkbox.isChecked(),getattr(self,pname)))
            self.bgui.add(pname,checkbox)
        self.bwbFileEntry(box,button,ledit,layout=layout,label=pvalue['label']+':', entryType=pvalue['type'], checkbox=checkbox)
        self.bgui.add(pname,ledit)
        self.bgui.add(pname,button,enableCallback=lambda value, clearLedit: self.enableFileDir(value, clearLedit, checkbox,ledit,button))
        if addCheckbox:
            self.updateCheckbox(pname,checkbox.isChecked(),ledit.text())
            
    def enableFileDir(self,value, clearLedit, cb,ledit,btn):
        if cb:
            cb.setEnabled(True)
            if cb.isChecked():
                ledit.setEnabled(True)
                btn.setEnabled(True)
            else:
                ledit.setEnabled(False)
                btn.setEnabled(False)                
        else:
            ledit.setEnabled(True)
            btn.setEnabled(True)
        if value is None or clearLedit:
            ledit.clear()
            
    def drawFilesBoxBtnRules(self,boxEdit,ledit,addBtn,removeBtn):
        if not ledit.text():
            ledit.clear()
            addBtn.setEnabled(False)
        if not boxEdit.selectedItems():
            removeBtn.setEnabled(False)
            
    def enableFilesBox(self, value, clearLedit, checkbox,browseBtn,boxEdit,ledit,addBtn,removeBtn):
        #first element is checkbox
        #last element is browseBtn if it exists
        if checkbox:
            checkbox.setEnabled(True)
        ledit.clear()
        boxEdit.clear()
        if value:
            boxEdit.addItems(str.splitlines(value))
        if not checkbox or checkbox.isChecked():
            for g in  [boxEdit,ledit,browseBtn,addBtn,removeBtn]:
                if g:
                    g.setEnabled(True)          
        self.drawFilesBoxBtnRules(boxEdit,ledit,addBtn,removeBtn)
        
    def updateFilesBox(self,attr,ledit,boxEdit): #updates for input - called before and after addition and deletion of input
        if hasattr(self,attr):
            value=getattr(self,attr)
            ledit.clear()
            if value is None: #explicitly check for this to avoid text None from appearing
                boxEdit.clear()
            else:
                boxEdit.clear()
                boxEdit.addItems(str.splitlines(value))
                self.updateBoxEditValue(attr,boxEdit)
            
    def updateBoxEditValue(self,attr,boxEdit):
        myItems=[]
        sys.stderr.write('{} objects in {}\n'.format(boxEdit.count(),boxEdit))
        for i in range(boxEdit.count()):
            myItems.append(boxEdit.item(i).text())
            sys.stderr.write('add item number {} value {}\n'.format(i,boxEdit.item(i).text()))
        if myItems:
            setattr(self,attr,"\n".join(myItems))
        
    def drawFilesBox (self,pname, pvalue, box=None, layout=None, addCheckbox=False):
        #for multiple files or directories draw a widget list with a line edit below and buttons to browse, add, delete, submit and reload
        #is used for multiple entry of any type
        
        disabledFlag=False
        checkbox=None
        browseBtn=None
        
        if addCheckbox or self.inputConnections.isConnected(pname):
            disabledFlag=True
        
        elements=[]
        #add checkbox if necessary
        if addCheckbox:
            value=None
            if hasattr(self,pname):
                value=getattr(self,pname)
            checkAttr=pname+'Checked' #this is not actually used but needed for orange gui checkbox element
            setattr(self,checkAttr,self.optionsChecked[pname])
            checkbox=gui.checkBox(None, self,checkAttr,label=None)
            checkbox.stateChanged.connect(lambda : self.updateCheckbox(pname,checkbox.isChecked(),value=value))
            elements.append(checkbox)
            
        #line edit for manual file entry
        leditAttr=pname+'Ledit'
        if not hasattr(self,leditAttr):
            setattr(self,leditAttr,None);
        ledit = gui.lineEdit(None, self,leditAttr,disabled=addCheckbox)
        ledit.setClearButtonEnabled(True)
        if pvalue['type'] == 'file list':
            ledit.setPlaceholderText("Enter file")
        elif pvalue['type'] == 'directory list':
            ledit.setPlaceholderText("Enter directory")
        else:
            ledit.setPlaceholderText("Enter parameter")
        
        ledit.setStyleSheet(":disabled { color: #282828}")
        ledit.clear()
        elements.append(ledit)

        #setup boxEdit
        boxEdit=QtGui.QListWidget(self)
        boxEdit.setDragDropMode(QtWidgets.QAbstractItemView.InternalMove)
        boxEdit.setSelectionMode(QtWidgets.QAbstractItemView.ExtendedSelection)
        boxEdit.setStyleSheet(":disabled { color: #282828}")
        boxEdit.setMinimumHeight(60)
        #fill boxEdit
        if hasattr(self,pname):
            fileStr=getattr(self,pname)
            if fileStr:
                boxEdit.addItems(str.splitlines(fileStr))
        boxEdit.setDisabled(disabledFlag)
        elements.append(boxEdit)
        
        #buttons
        addBtn=gui.button(None, self, "", callback=lambda: self.addLedit(pname,ledit,boxEdit,addBtn), autoDefault=False,disabled=disabledFlag)
        removeBtn=gui.button(None, self, "", callback=lambda: self.removeItem(pname,boxEdit,removeBtn), autoDefault=False,disabled=disabledFlag)
        if  pvalue['type'] == 'file list':
            browseBtn = gui.button(None, self, "", callback=partial(self.browseFiles, boxEdit=boxEdit,attr=pname), autoDefault=False,disabled=disabledFlag)
        elif pvalue['type'] == 'directory list':
            browseBtn = gui.button(None, self, "", callback=partial(self.browseDirs, boxEdit=boxEdit,attr=pname), autoDefault=False,disabled=disabledFlag)
        #set styles
        buttonStyle='background: None; border: None ; border-radius: 0;'
        
        
        for btn in (addBtn,removeBtn,browseBtn):
            if btn:
                btn.setStyleSheet(buttonStyle)
                elements.append(btn)      
     
        #set icons
        self.bgui.addList(pname,elements,enableCallback=lambda value, clearLedit: self.enableFilesBox(value,clearLedit, checkbox,browseBtn,boxEdit,ledit,addBtn,removeBtn),updateCallback=lambda: self.updateFilesBox(pname,ledit,boxEdit))
        if browseBtn:
            browseBtn.setIcon(self.browseIcon)
        addBtn.setIcon(self.addIcon)
        removeBtn.setIcon(self.removeIcon)
        
        #check rules for buttons    
        
        self.drawFilesBoxBtnRules(boxEdit,ledit,addBtn,removeBtn)
        
        #connects from ledit and boxEdit to buttons
        ledit.textChanged.connect(lambda: addBtn.setEnabled(bool(ledit.text())))
        boxEdit.itemSelectionChanged.connect(lambda: removeBtn.setEnabled(bool(boxEdit.selectedItems())))
        
        #layout section
        filesBoxLeditLayout=QtGui.QVBoxLayout()
        #add to the main parameters box
        myBox=gui.vBox(None)
        myBox.layout().addLayout(filesBoxLeditLayout)
        startCol=0
        label=QtGui.QLabel(pvalue['label']+':')
        label.setAlignment(Qt.AlignTop)
        layout.addWidget(label,layout.nextRow,startCol)
        layout.addWidget(myBox,layout.nextRow,1,1,2)
        
        layout.nextRow = layout.nextRow + 1
        #line layout     
        lineLayout=QtGui.QGridLayout()
        
        if addCheckbox:
            lineLayout.addWidget(checkbox,1,0)
            startCol=1
        
        lineLayout.addWidget(ledit,1,startCol+1)
        if browseBtn:
            lineLayout.addWidget(browseBtn,1,startCol+2)
        else:
            startCol=startCol-1
        lineLayout.addWidget(addBtn,1,startCol+3)
        lineLayout.addWidget(removeBtn,1,startCol+4)
       
        #now add the two layouts to the bigBox layout
        filesBoxLeditLayout.addWidget(boxEdit)
        filesBoxLeditLayout.addLayout(lineLayout)
        self.updateBoxEditValue(pname,boxEdit)
    
    def addLedit(self,attr,ledit,boxEdit,addBtn):
        #adds text in ledit to items in boxEdit
        if ledit.text():
            boxEdit.addItem(ledit.text())
            sys.stderr.write('adding {} to {}\n'.format(ledit.text(),attr))
            ledit.clear()
            addBtn.setEnabled(False)
            self.updateBoxEditValue(attr,boxEdit)
            sys.stderr.write('boxEdit values is {}\n'.format(getattr(self,attr)))
            

    def removeItem(self,attr,boxEdit,removeBtn):
        if boxEdit.selectedItems():
            for item in boxEdit.selectedItems():
                boxEdit.takeItem(boxEdit.row(item))
            self.updateBoxEditValue(attr,boxEdit)
        if not boxEdit.count():
            removeBtn.setEnabled(False)

        
    def drawExec(self, box=None):
        #find out if there are triggers
        self.candidateTriggers=[]
        for pname in self.data['inputs']:
            if not(pname in self.data['requiredParameters']):
                self.candidateTriggers.append(pname)

        #initialize the exec state
        self.execLayout=QtGui.QGridLayout()
        self.execLayout.setSpacing(5)
        self.cboRunMode=QtGui.QComboBox()
        self.cboRunMode.addItem('Manual')
        self.cboRunMode.addItem('Automatic')
        if self.candidateTriggers:
            self.cboRunMode.addItem('Triggered')
        elif self.runMode == 2: #reset in case there is a change in an older workflow
            self.runMode=0
        self.cboRunMode.setCurrentIndex(self.runMode)
        if self.candidateTriggers:
            self.execBtn=QtGui.QToolButton(self)
            self.execBtn.setText('Select Triggers')
            self.execMenu=QtGui.QMenu(self)
            self.triggerMenu={}
            for attr in self.candidateTriggers:
                action=self.execMenu.addAction(attr)
                action.setCheckable(True)
                action.setChecked( bool(attr in self.runTriggers))
                action.changed.connect(self.chooseTrigger)
                self.triggerMenu[action]=attr
            self.execBtn.setMenu(self.execMenu)
            self.execBtn.setPopupMode(QtGui.QToolButton.InstantPopup)
            if self.runMode ==2:
                self.execBtn.setEnabled(True)
            else:
                self.execBtn.setEnabled(False)
        self.cboRunMode.currentIndexChanged.connect(self.runModeChange)
        myLabel=QtGui.QLabel('RunMode:')
        myLabel.setAlignment(Qt.AlignRight | Qt.AlignVCenter)
        css = '''
        QPushButton {background-color: #1588c5; color: white; height: 20px; border: 1px solid #1a8ac6; border-radius: 2px;}
        QPushButton:pressed { background-color: #158805; border-style: inset;}
        QPushButton:disabled { background-color: lightGray; border: 1px solid gray; }
        QPushButton:hover {background-color: #1588f5; }
        '''
        self.btnRun = gui.button(None, self, "Start", callback=self.OnRunClicked)
        self.btnRun.setStyleSheet(css)
        self.btnRun.setFixedSize(60,20)
        self.execLayout.addWidget(self.btnRun,1,0)
        self.execLayout.addWidget(myLabel,1,1)
        self.execLayout.addWidget(self.cboRunMode,1,2)
        if self.candidateTriggers:
            self.execLayout.addWidget(self.execBtn,1,3)
        box.layout().addLayout(self.execLayout)

    def runModeChange(self):
        self.runMode=self.cboRunMode.currentIndex()
        if self.candidateTriggers:
            if self.runMode == 2:
                self.execBtn.setEnabled(True)
            else:
                self.execBtn.setEnabled(False)
        self.checkTrigger()

    def chooseTrigger(self):
        action=self.execMenu.sender()
        attr=self.triggerMenu[action]
        checked=action.isChecked()
        if attr is None or checked is None:
            return
        if checked and attr not in self.runTriggers:
            self.runTriggers.append(attr)
        if not checked and attr in self.runTriggers:
            self.runTriggers.remove(attr)

    def checkTrigger(self):
        #this should be checked any time there is a change
        if self.runMode ==0: #manual - only go on start button
            return
        elif self.runMode ==1: #automatic same as pushing start button
            self.OnRunClicked()
        elif self.candidateTriggers:
            #check if the input triggers are set
            for trigger in self.runTriggers:
                if not inputConnections.isSet(trigger):
                    return
            self.OnRunClicked()

    def bwbFileEntry(self, widget, button, ledit, icon=browseIcon,layout=None, label=None,entryType='file', checkbox=None):
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
   
    def browseFiles(self,boxEdit,attr=None):
        defaultDir = '/root'
        if os.path.exists('/data'):
            defaultDir = '/data'
        boxEdit.addItems(QtWidgets.QFileDialog.getOpenFileNames(self, caption="Choose files", directory=defaultDir, filter="Any file (*.*)")[0])
        if attr:
            self.updateBoxEditValue(attr,boxEdit)
        
    def browseDirs(self,boxEdit,attr=None):
        dlg=getExistingDirectories()
        if dlg.exec_() == QtWidgets.QDialog.Accepted:
            self.boxEdit.addItems(dlg.selectedFiles())
            if attr:
                self.updateBoxEditValue(attr,boxEdit)

    def browseFileDir(self, attr, filetype=None):
        defaultDir = '/root'
        if os.path.exists('/data'):
            defaultDir = '/data'
        if filetype =='file':
            myFile=QtWidgets.QFileDialog.getOpenFileName(self, "Locate file", defaultDir)[0]
            if myFile:
                setattr(self,attr,myFile)
                dirAttr=attr+"Dir"
                setattr(self,dirAttr,os.path.dirname(myFile))
        else:
            myDir=QtWidgets.QFileDialog.getExistingDirectory(self, caption="Locate directory", directory=defaultDir)
            if myDir:
                setattr(self,attr,myDir)
    
    def updateCheckbox(self, pname, state, value=None):
        self.optionsChecked[pname]=state
        sys.stderr.write('updating checkbox pname {} connect {} isChecked {}\n'.format(pname,self.inputConnections.isConnected(pname),state))
        if not (self.inputConnections.isConnected(pname)) and state:
            self.bgui.enable(pname,value)
            sys.stderr.write('enabled\n')
        elif not (self.inputConnections.isConnected(pname)) and not state:
            self.bgui.disable(pname,ignoreCheckbox=True)
            sys.stderr.write('disabled\n')
    
    def enableSpin(self,value,clearLedit, guiSpin):
        (cb,spin)=guiSpin
        if cb.isChecked():
            spin.setEnabled(True)
        else:
            spin.setEnabled(False)
        cb.setEnabled(True)
            
    def updateSpinCheckbox(self,pname):
        checkAttr=pname+'Checked'
        self.optionsChecked[pname]=getattr(self,checkAttr)
        
#Handle inputs
    def handleInputs(self, value, attr, sourceId=None):
        sys.stderr.write('handler: attr {} value {}\n'.format(attr,value))
        if value is None:
            self.inputConnections.remove(attr,sourceId)
            sys.stderr.write('removing {} disabled {}\n'.format(attr,self.inputConnections.isConnected(attr)))
            setattr(self,attr,None) #this gives text None in ledit for some reason
        else:
            self.inputConnections.add(attr,sourceId)
            sys.stderr.write('handler adding input: attr {} value {}\n'.format(attr,value))
            setattr(self,attr,value)
            self.checkTrigger()
        self.updateGui(attr,value)

    def updateGui(self,attr,value,removeFlag=False):
        if self.inputConnections.isConnected(attr):
            sys.stderr.write('disabling {}\n'.format(attr))
            self.bgui.disable(attr,value)
        else:
            sys.stderr.write('enabling {} with {}\n'.format(attr,value))
            self.bgui.enable(attr,value)

#Generate commands and run job

    def startJob(self):
        self.hostVolumes = {}
        #check for missing parameters and volumes
        missingParms=self.checkRequiredParms()
        if missingParms:
            self.infoLabel.setText("missing required parameters: {}".format(missingParms))
            return
        missingVols=self.getRequiredVols()
        if missingVols:
            self.infoLabel.setText("missing or incorrect volume mappings to: {}".format(missingVols))
            return
        #get ready to start
        attrList=self.__dict__.keys()
        self.bgui.disableAll()
        self.disableExec()
        cmd=self.generateCmdFromData()
        self.envVars={}
        self.getEnvironmentVariables()
        sys.stderr.write('envs {}\n'.format(self.envVars))
        sys.stderr.write('volumes {}\n'.format(self.hostVolumes))
        sys.stderr.write('cmd {}\n'.format(cmd))
        try:
            self.dockerRun(self.hostVolumes,cmd,environments=self.envVars)
        except:
            self.bgui.reenableAll(self)
            self.reenableExec()

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
            if self.data['parameters'][attr]['type'] =='file':
                bwbVol=os.path.dirname(os.path.normpath(bwbVol))
            elif self.data['parameters'][attr]['type'] =='files' or self.data['parameters'][attr]['type'] =='file list':
                files=str.splitlines(bwbVol)
                bwbVol=self.findTopDirectory(files)
            else:
               bwbVol=os.path.normpath(bwbVol)
            self.hostVolumes[conVol]=bwbVol
        return None

    def generateCmdFromData (self):
        flags={}
        args=[]
        for pname, pvalue in self.data['parameters'].items():
           
            #possible to have an requirement or parameter that is not in the executable line
            #this is indicated by the absence of a flags field
            if 'flags' not in pvalue:
                continue
            #if required or checked then it is added to the flags
            addFlag=False
            #checkattr is needed for the orange gui checkboxes but is not otherwise updated 
            if pname in self.optionsChecked and self.optionsChecked[pname]:
                addFlag=True
            #also need to check for booleans which are not tracked by optionsChecked
            
            if pvalue['type'] == 'bool':
                if hasattr(self,pname) and getattr(self,pname):
                    addFlag=True
            if pname in self.data['requiredParameters']:
                sys.stderr.write('{} with is required\n'.format(pname))
                if hasattr(self,pname):
                    sys.stderr.write('{} value is {}\n'.format(pname,getattr(self,pname)))
                addFlag=True
            else:
                sys.stderr.write('{} is not required\n'.format(pname))
                if hasattr(self,pname):
                    sys.stderr.write('{} value is {}\n'.format(pname,getattr(self,pname)))
            sys.stderr.write('gencmd name {} value {} addFlag {} type {} flags {}\n'.format(pname,pvalue,addFlag,pvalue['type'],pvalue['flags']))
            if addFlag:
                if pvalue['flags'] and pvalue['flags'][0]:
                    if pvalue['type'] == 'bool':
                        flags[pvalue['flags'][0]] = None
                    elif pvalue['type'] == 'file':
                        filename=str(getattr(self,pname))
                        if filename:
                            hostFilename=self.bwbPathToContainerPath(filename, isFile=True,returnNone=False)
                            flags[pvalue['flags'][0]]=str(hostFilename)
                    elif pvalue['type'] == 'file list':
                        files=str.splitlines(getattr(self,pname))
                        if files:
                            hostFiles=[]
                            for f in files:
                                hostFiles.append(self.bwbPathToContainerPath(f, isFile=True,returnNone=False))
                            flags[pvalue['flags'][0]]="\n".join(hostFiles)
                    elif pvalue['type'] == 'directory':
                        path=str(getattr(self,pname))
                        if path:
                            hostPath=self.bwbPathToContainerPath(path, returnNone=False)
                            flags[pvalue['flags'][0]]=str(hostPath)
                    elif pvalue['type'][-4:] =='list':
                        flags[pvalue['flags'][0]]=getattr(self,pname)
                    else:
                        flags[pvalue['flags'][0]]=str(getattr(self,pname))
                else:
                    if pvalue['type'] == 'file':
                        filename=str(getattr(self,pname))
                        if filename:
                            hostFilename=self.bwbPathToContainerPath(filename, isFile=True,returnNone=False)
                            args.append(hostFilename)
                    elif pvalue['type'] =='file list':
                        files=str.splitlines(getattr(self,pname))
                        sys.stderr.write('files are {}\n'.format(files))
                        for f in files:
                            args.append(self.bwbPathToContainerPath(f, isFile=True,returnNone=False))
                    elif pvalue['type'] =='directory':
                        path=str(getattr(self,pname))
                        if path:
                            hostPath=self.bwbPathToContainerPath(path, returnNone=False)
                            args.append(hostPath)
                    elif pvalue['type'][-4:] =='list' or pvalue['type'][-4:] =='List':
                        myItems=str.splitlines(getattr(self,pname))
                        sys.stderr.write('list value name {} values {}\n'.format(pname,myItems))
                        args.extend(myItems)
                    else:
                        args.append(str(getattr(self,pname)))
        return self.generateCmdFromBash([self.data['command']],flags,args=args)

    def generateCmdFromBash (self, executables = [], flags = {}, args = []):
        sys.stderr.write('executables {} flags {} args {}\n'.format(executables,flags,args))
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
            if (flagName and flagName[:1] == '-') and (flagName[:2] != '--') and (not flagValue):
                cmdStr +=  flagName + ' '
        for flagName,flagValue in flags.items():
            if (flagName and flagName[:1] == '-') and (flagName[:2] != '--') and flagValue:
                if(flagValue[:1] == '='):
                    cmdStr +=  flagName + flagValue + ' '
                else:
                    cmdStr +=  flagName + ' ' + flagValue + ' '
        for flagName, flagValue in flags.items():
            if (flagName and flagName[:2] == '--')  and (not flagValue):
                cmdStr +=  flagName + ' '
        for flagName,flagValue in flags.items():
            if (flagName and flagName[:2] == '--')  and flagValue:
                if(flagValue[:1] == '='):
                    cmdStr +=  flagName + flagValue + ' '
                else:
                    cmdStr +=  flagName + ' ' + flagValue + ' '
        for flagName, flagValue in flags.items():
            if (flagName and flagName[:1] != '-')  and (not flagValue):
                cmdStr +=  flagName + ' '
        for flagName,flagValue in flags.items():
            if (flagName and flagName[:1] != '-')  and flagValue:
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

    def getEnvironmentVariables(self):
        #dynamic environment variables
        for pname in self.data['parameters']:
            pvalue=self.data['parameters'][pname]
            sys.stderr.write('checking var {} with value {} for env {}\n'.format(pname,getattr(self,pname),list(pvalue.keys())))
            if 'env' in pvalue and getattr(self,pname) is not None:
                self.envVars[pvalue['env']] = getattr(self,pname)
                sys.stderr.write('var {} env {} assigned to {}\n'.format(pname,pvalue['env'],getattr(self,pname)))
        #now assign static environment variables
        if 'env' in self.data:
            for e in self.data['env']:
                if e not in self.envVars:
                    self.envVars[e]=self.data['env'][e]

#Docker run/pull (Jiamings code)
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

#Event handlers
    def OnRunClicked(self):
        if hasattr (self,'userStartJob'):
            self.userStartJob()
        else:
            self.startJob()

    def Event_OnRunFinished(self):
        self.infoLabel.setText("Finished")
        self.bgui.reenableAll(self)
        self.reenableExec()
        self.handleOutputs()

    def Event_OnRunMessage(self, message):
        self.infoLabel.setText(message)

#Utilities
    def bwbPathToContainerPath(self, path, isFile=False,returnNone=False):
        #converts the path entered relative to the bwb mountpoint
        #will return None if not found or the original path (default) depending on variable
        #first map it to the  path
        pathFile=None
        if isFile:
            dirPath=os.path.dirname(path)
            pathFile=os.path.basename(path)
            hostPath=os.path.normpath(self.dockerClient.to_best_host_directory(dirPath,returnNone=False))
        else:
            hostPath=os.path.normpath(self.dockerClient.to_best_host_directory(path,returnNone=False))
        conPath=None
        #now get all the possible submappings to volumeMappings by comparing the true hostmappings
        #if submap is found convert the common path to the container path
        #return shortest path
        for mapping in self.data['volumeMappings']:
            conVol = mapping['conVolume']
            bwbVol = self.hostVolumes[conVol]
            hostVol=os.path.normpath(self.dockerClient.to_best_host_directory(bwbVol,returnNone=False))
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
        if conPath is not None :
            if isFile:
                return os.path.normpath(str.join(os.sep,(conPath, pathFile)))
            return conPath
        else:
            if returnNone:
                return conPath
            return path

    def findTopDirectory(self,files):
        #return the top directory in a set of files
        #actually returns the shortest directory
        bestPath=None
        for f in files:
            bwbVol=os.path.dirname(os.path.normpath(f))
            if bestPath is None:
                bestPath = bwbVol
            elif len(bwbVol) < len(bestPath):
                bestPath=bwbVol
        return bestPath

    def disableExec(self):
        self.btnRun.setText('Running')
        self.btnRun.setEnabled(False)
        self.cboRunMode.setEnabled(False)
            
    def reenableExec(self):
        self.btnRun.setText('Start')
        self.btnRun.setEnabled(True)
        self.cboRunMode.setEnabled(True)
