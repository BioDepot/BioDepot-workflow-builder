import os
import re
import sys
import logging
from functools import partial
from pathlib import Path
from AnyQt.QtCore import QThread, pyqtSignal, Qt
from Orange.widgets import widget, gui, settings
from orangebiodepot.util.DockerClient import DockerClient, PullImageThread, ConsoleProcess
from PyQt5 import QtWidgets, QtGui, QtCore


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
        self.findChildren(QtWidgets.QTreeView)[0].setSelectionMode(QtWidgets.QAbstractItemView.ExtendedSelection)
class BwbGridLayout():
    #adds methods to keep track of rows and columns
    def __init__(self, spacing=5, startCol=0,startRow=0):
        self._layout=QtGui.QGridLayout()
        self._layout.setSpacing(spacing)
        self.nextCol=startCol
        self.nextRow=startRow

    def addWidget(self,widget,space=None,width=None,height=1,linefeed=None):
        if space is not None:
            self.nextCol+=space
        if width is not None:
            self._layout.addWidget(widget,self.nextRow,self.nextCol,width,height)
            self.nextCol+=width
        else:
            self._layout.addWidget(widget,self.nextRow,self.nextCol)
            self.nextCol+=1
        if linefeed is not None:
            self.nextRow+=linefeed
            self.nextCol=0
            
    def addSpace(self,space=1):
        self.nextCol+=space
        
    def addLinefeed(self,lf=1,startCol=0):
        self.nextCol=startCol
        self.nextRow+=lf    
    
    def layout(self):
        return self._layout
            
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
        sys.stderr.write('checking attr {}\n'.format(attr))
        clearLedit=False
        if value is None:
            clearLedit=True
        if attr in self._updateCallbacks:
            sys.stderr.write('applying update for attr {}\n'.format(attr))
            self._updateCallbacks[attr]() 
        if attr in self._dict:
            sys.stderr.write('found attr in dict {}\n'.format(attr))
            if attr in self._enableCallbacks:
                sys.stderr.write('enable callback for {}\n'.format(attr))
                self._enableCallbacks[attr](value,clearLedit)
            else:
                sys.stderr.write('enable for {}\n'.format(attr))
                for g in self._dict[attr]:
                    sys.stderr.write('enable element {}\n'.format(g))
                    self.enableElement(g)
        return True
        
    def disableAll(self):
        for attr in self._dict.keys():
            self.disable(attr)

            
    def reenableAll(self,OWself):
        for attr in self._dict.keys():
            if not OWself.inputConnections.isConnected(attr) :
                self.enable(attr,getattr(OWself,attr))

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

    def isConnected (self, slot, connectionId=None):
        if self._dict is None:
            return False
        if slot in self._dict:
            if connectionId:
                if connectionId in self._dict[slot]:
                    return True
            else:
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
        self.inputConnections=ConnectionDict(self.inputConnectionsStore)
        self._dockerImageName = image_name
        self._dockerImageTag = image_tag
        #drawing layouts for gui
        #file directory
        self.filesBoxLayout=QtGui.QVBoxLayout()
        self.fileDirRequiredLayout=BwbGridLayout()
        self.fileDirOptionalLayout=BwbGridLayout()

        #lineEdits
        self.leditRequiredLayout=BwbGridLayout()
        self.leditOptionalLayout=BwbGridLayout()
        
        #keep track of gui elements associated with each attribute
        self.bgui=BwbGuiElements()
        #For compatibility if triggers are not being kept
        if not hasattr(self,'triggerReady'):
            self.triggerReady={}
        
    def initVolumes(self):
        #initializes container volumes
        if 'volumeMappings' in self.data and self.data['volumeMappings']:
            for mapping in self.data['volumeMappings']:
                bwbVolAttr=mapping['attr']
                if not hasattr(self, bwbVolAttr) :
                    setattr(self,bwbVolAttr,None)


#Drawing the GUI
    def drawGUI(self):
        
        self.setStyleSheet(":disabled { color: #282828}")
        self.topScrollArea = QScrollArea(
            verticalScrollBarPolicy=Qt.ScrollBarAlwaysOn
        )
        self.topBox=gui.widgetBox(self.controlArea)
        self.topScrollArea.setWidget(self.topBox)
        self.topScrollArea.setWidgetResizable(True)
        self.controlArea.layout().addWidget(self.topScrollArea)
        if 'requiredParameters' in self.data and self.data['requiredParameters']:
            self.requiredBox = gui.widgetBox(self.topBox, "Required parameters")
            self.requiredBox.layout().addLayout(self.fileDirRequiredLayout.layout())
            setattr(self.fileDirRequiredLayout.layout(),'added',True)
            self.requiredBox.layout().addLayout(self.leditRequiredLayout.layout())
            setattr(self.leditRequiredLayout.layout(),'added',True)              
            self.drawRequiredElements()
            
        if self.findOptionalElements():
            self.optionalBox = gui.widgetBox(self.topBox, "Optional parameters")
            self.optionalBox.layout().addLayout(self.fileDirOptionalLayout.layout())
            setattr(self.fileDirOptionalLayout.layout(),'added',True)
            self.optionalBox.layout().addLayout(self.leditOptionalLayout.layout())
            setattr(self.leditOptionalLayout.layout(),'added',True)              
            self.drawOptionalElements()
            
        #disable connected elements
        for i in self.inputs:
            attr=i.name
            if self.inputConnections.isConnected(attr):
               self.bgui.disable(attr)
        
        #make a box for the console and console control - not abs necessary but it might help with future org and with the size hinting system
        self.consoleBox = gui.widgetBox(self.controlArea)
        consoleControlLayout=BwbGridLayout()
        self.consoleBox.layout().addLayout(consoleControlLayout.layout())
        setattr(consoleControlLayout,'added',True)        
        self.drawConsoleControl(box=self.consoleBox,layout=consoleControlLayout)
        self.console=QtGui.QTextEdit()
        self.console.setReadOnly(True)
        pal=QtGui.QPalette()
        pal.setColor(QtGui.QPalette.Base,Qt.black)
        pal.setColor(QtGui.QPalette.Text,Qt.green)
        self.console.setPalette(pal)
        self.console.setAutoFillBackground(True)
        self.consoleBox.layout().addWidget(self.console)
        controlBox = QtGui.QVBoxLayout()
        self.pConsole=ConsoleProcess(console=self.console,finishHandler=self.onRunFinished)
        self.drawExec(box=self.controlArea.layout())
        self.checkTrigger()

    def drawConsoleControl(self,box=None,layout=None):
        css = '''
        QPushButton {background-color: #1588c5; color: white; height: 20px; border: 1px solid #1a8ac6; border-radius: 2px;}
        QPushButton:pressed { background-color: #158805; border-style: inset;}
        QPushButton:disabled { background-color: lightGray; border: 1px solid gray; }
        QPushButton:hover {background-color: #1588f5; }
        '''
        layout.addWidget(QtGui.QLabel('Console: '))
        pname='saveLog'
        pvalue={ 'type' : 'directory', 'label' : 'AutoLog'}
        if not hasattr(self,pname):
            setattr(self,pname,None)

        self.btnConsoleClear = gui.button(None, self, "Clear", callback=self.clearConsole)
        self.btnConsoleClear.setStyleSheet(css)
        self.btnConsoleClear.setFixedSize(60,20)
        self.btnConsoleSave = gui.button(None, self, "Save", callback=self.saveConsole)
        self.btnConsoleSave.setStyleSheet(css)
        self.btnConsoleSave.setFixedSize(60,20)
        layout.addWidget(self.btnConsoleClear)
        layout.addWidget(self.btnConsoleSave)
        self.drawFileDirElements(pname, pvalue, box=box,layout=layout, addCheckbox=True)        
        #keep track of buttons separately because the drawFileDir routine has its own callback for enabling the elements
        #should fix this at some point to check which elements have been dealt with using a dict or work out a recursive scheme with elements 
        btnPname=pname+'Btns'
        if not hasattr(self,btnPname):
            setattr(self,btnPname,None)
        self.bgui.add(btnPname,self.btnConsoleClear)
        self.bgui.add(btnPname,self.btnConsoleSave)
        
    def clearConsole(self):
        self.console.clear()
        
    def saveConsole(self):
        defaultDir = '/root'
        if os.path.exists('/data'):
            defaultDir = '/data'
        fileName = QtWidgets.QFileDialog.getSaveFileName(self,"QFileDialog.getSaveFileName()","","Text files (*.txt);;All Files (*)")[0]
        with open(fileName, 'w') as myFile:
            myFile.write(str(self.console.toPlainText()))
            myFile.close()


        #consoleControlLayout.addWidget(outputLabel,0,0) 
    def drawElements(self, elementList,isOptional=False):
        for pname in elementList:
            if not ('parameters' in self.data) or not ( pname in self.data['parameters']):
                continue
            pvalue=self.data['parameters'][pname]
            if ('gui' in pvalue and pvalue['gui'] != 'file' and pvalue['gui'] != 'directory') or ( pvalue['type'] != 'file' and pvalue['type'] != 'directory'):
                continue
            if isOptional:
                self.drawFileDirElements(pname, pvalue, box=self.optionalBox,layout=self.fileDirOptionalLayout,addCheckbox=True)
            else:
                self.drawFileDirElements(pname, pvalue, box=self.requiredBox,layout=self.fileDirRequiredLayout)
            
        for pname in elementList:
            if not ('parameters' in self.data) or not ( pname in self.data['parameters']):
                continue
            pvalue=self.data['parameters'][pname]
            if ('gui' in pvalue and pvalue['gui'][-4:] != 'list') or ( pvalue['type'][-4:] != 'list' ):
                continue
            sys.stderr.write('drawing textBox for  pname {} pvalue {}\n'.format(pname,pvalue))
            if isOptional:
                self.drawTextBox (pname, pvalue, box=self.optionalBox,layout=self.fileDirOptionalLayout,addCheckbox=True)
            else:
                self.drawTextBox(pname, pvalue, box=self.requiredBox,layout=self.fileDirRequiredLayout)

        for pname in elementList:
            if not ('parameters' in self.data) or not( pname in self.data['parameters']):
                continue
            pvalue=self.data['parameters'][pname]
            if ('gui' in pvalue and pvalue['gui'] != 'Ledit') or (pvalue['type'] != 'double' and pvalue['type'] != 'str' and pvalue['type'] != type('text')):
                sys.stderr.write('type is {} {}\n'.format(pvalue['type'],type('text')))
                continue
            if isOptional:
                self.drawLedit(pname,pvalue,self.optionalBox,layout=self.leditOptionalLayout,addCheckbox=True)
            else:
                self.drawLedit(pname,pvalue,self.requiredBox,layout=self.leditRequiredLayout)
                
        for pname in elementList:
            if not ('parameters' in self.data) or not ( pname in self.data['parameters']):
                continue
            pvalue=self.data['parameters'][pname]
            if ('gui' in pvalue and pvalue['gui'] != 'Spin') or (pvalue['type'] != 'int'):
                continue
            if isOptional:
                self.drawSpin(pname,pvalue,self.optionalBox,addCheckbox=True)
            else:
                self.drawSpin(pname,pvalue,self.requiredBox)

        for pname in elementList:
            if not ('parameters' in self.data) or not ( pname in self.data['parameters']):
                continue
            pvalue=self.data['parameters'][pname]
            if ('gui' in pvalue and pvalue['gui'] != 'bool') or (pvalue['type'] != 'bool'):
                continue
            if isOptional:
                self.drawCheckbox(pname,pvalue,self.optionalBox)
            else:
                self.drawCheckbox(pname,pvalue,self.requiredBox)
        
    def drawRequiredElements(self):
        for pname in self.data['requiredParameters']:
            if not ('parameters' in self.data) or not (pname in self.data['parameters']):
                continue
            pvalue=self.data['parameters'][pname]
            if not hasattr(self,pname):
                setattr(self,pname,None)
            if (getattr(self,pname) is None) and ('default' in pvalue):
                setattr(self,pname,pvalue['default'])
        self.drawElements(self.data['requiredParameters'])


    
    def findOptionalElements(self):
       #checks that there are optional elements
        if not 'parameters' in self.data or not self.data['parameters']:
            return False
        for pname in self.data['parameters']:
            if pname not in self.data['requiredParameters']:
                return True
        return False
 
    def drawOptionalElements(self):
        optionalList=[]
        for pname in self.data['parameters']:
            if pname not in self.data['requiredParameters']:
                pvalue=self.data['parameters'][pname]
                if not hasattr(self,pname):
                    setattr(self,pname,None)
                if (getattr(self,pname) is None) and ('default' in pvalue):
                    setattr(self,pname,pvalue['default'])
                if not (pname in self.optionsChecked):
                    self.optionsChecked[pname]=False
                optionalList.append(pname)    
        self.drawElements(optionalList,isOptional=True)
            

    def drawCheckbox(self,pname,pvalue, box=None):
        #for booleans - their value is the same as the checkbox state
        sys.stderr.write('drawCB pname {} pvalue {} label {}\n'.format(pname, pvalue,pvalue['label']))
        if not hasattr(self,pname):
            if 'default' in pvalue:
                if type (pvalue['default']) is str:
                    if pvalue['default'] ==  'True':
                       setattr(self,pname,True)
                    if pvalue['default'] ==  'False':
                       setattr(self,pname,False)
                    else:
                        raise Exception ('{} is boolean - default values must be True or False not {}'.format(pname,pvalue['default']))
                elif type (pvalue['default']) is bool:
                    setattr(self,pname,pvalue['default'])
            else:
                setattr(self,pname,False)
        sys.stderr.write('draw CB pname {} value {}\n'.format(pname,getattr(self,pname)))
        cb=gui.checkBox(box, self, pname, pvalue['label'])
        checkAttr=pname+'Checked'
        setattr(self,checkAttr,getattr(self,pname))
        #check if inactive
        self.bgui.add(pname,cb)
            
    def drawSpin(self,pname,pvalue, box=None,addCheckbox=False):
        #for drawSpin - we use the origin version which already has a checkbox connected
        #TODO could change this to the same way we handle ledits with separate cbo
        #the gui spin box returns either (cb,sbox) or just sbox depending on whether there is a checkabox
        checkBox=None
        checkAttr=None
        if addCheckbox:
            if pname not in self.optionsChecked:
                self.optionsChecked[pname]=False
            checkAttr=pname+'Checked'  
            setattr(self,checkAttr,self.optionsChecked[pname]) 
        default=0
        if 'default' in pvalue:
            default=pvalue['default']
        if pvalue['type'] == 'int':
            if not hasattr(self,pname):
                setattr(self,pname,int(default))
            elif getattr(self,pname):
                setattr(self,pname,int(getattr(self,pname)))
            else:
                setattr(self,pname,int(default))
        else:
            if not hasattr(self,pname):
                setattr(self,pname,float(default))
            elif getattr(self,pname):
                setattr(self,pname,int(getattr(self,pname)))
            else:
                setattr(self,pname,float(default))
        if addCheckbox:
            (checkBox,mySpin)=gui.spin(box, self, pname, minv=1, maxv=128, label=pvalue['label'], checked=checkAttr, checkCallback=lambda : self.updateSpinCheckbox(pname))
        else:
             mySpin=gui.spin(box, self, pname, minv=1, maxv=128, label=pvalue['label'], checked=checkAttr, checkCallback=lambda : self.updateSpinCheckbox(pname))
           
        if getattr(self,pname) is None:
            mySpin.clear()
        self.bgui.add(pname,mySpin,enableCallback=lambda value,clearLedit  : self.enableSpin(value,clearLedit,checkBox,mySpin))
        
    def drawLedit(self,pname,pvalue,box=None,layout=None,addCheckbox=False):
        checkAttr=None
        checkbox=None
        ledit=gui.lineEdit(None, self, pname,disabled=addCheckbox)
        if addCheckbox:
            if pname not in self.optionsChecked:
                self.optionsChecked[pname]=False
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
        return {'ledit' : ledit, 'checkbox' : checkbox}
    
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
            if pname not in self.optionsChecked:
                self.optionsChecked[pname]=False
            setattr(self,checkAttr,self.optionsChecked[pname])
            checkbox=gui.checkBox(None, self,checkAttr,label=None)
            sys.stderr.write('updating filedir {}\n'.format(pname))
            checkbox.stateChanged.connect(lambda : self.updateCheckbox(pname,checkbox.isChecked(),getattr(self,pname)))
            self.bgui.add(pname,checkbox)
        labelValue=pvalue['label']
        if labelValue is None:
            labelValue=""
        sys.stderr.write('adding filedir for pname {} using layout {}\n'.format(pname,layout))    
        self.bwbFileEntry(box,button,ledit,layout=layout,label=labelValue+':', entryType=pvalue['type'], checkbox=checkbox)
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
            
    def drawTextBoxBtnRules(self,boxEdit,ledit,addBtn,removeBtn):
        if not ledit.text():
            ledit.clear()
            addBtn.setEnabled(False)
        if not boxEdit.selectedItems():
            removeBtn.setEnabled(False)
            
    def enableTextBox(self, value, clearLedit, checkbox,browseBtn,boxEdit,ledit,addBtn,removeBtn):
        #first element is checkbox
        #last element is browseBtn if it exists
        if checkbox:
            checkbox.setEnabled(True)
        ledit.clear()
        boxEdit.clear()
        if value:
            boxEdit.addItems(value)
        if not checkbox or checkbox.isChecked():
            for g in  [boxEdit,ledit,browseBtn,addBtn,removeBtn]:
                if g:
                    g.setEnabled(True)          
        self.drawTextBoxBtnRules(boxEdit,ledit,addBtn,removeBtn)
        
    def updateTextBox(self,attr,ledit,boxEdit): #updates for input - called before and after addition and deletion of input
        if hasattr(self,attr):
            value=getattr(self,attr)
            ledit.clear()
            if value is None: #explicitly check for this to avoid text None from appearing
                boxEdit.clear()
            else:
                boxEdit.clear()
                boxEdit.addItems(value)
                self.updateBoxEditValue(attr,boxEdit)
            
    def updateBoxEditValue(self,attr,boxEdit):
        myItems=[]
        sys.stderr.write('{} objects in {}\n'.format(boxEdit.count(),boxEdit))
        for i in range(boxEdit.count()):
            myItems.append(boxEdit.item(i).text())
            sys.stderr.write('add item number {} value {}\n'.format(i,boxEdit.item(i).text()))
        if myItems:
            setattr(self,attr,myItems)
        
    def drawTextBox (self,pname, pvalue, box=None, layout=None, addCheckbox=False):
        #for multiple files or directories draw a widget list with a line edit below and buttons to browse, add, delete, submit and reload
        #is used for multiple entry of any type
        
        sys.stderr.write('adding text box for pname{} pvalue {} box {} layout {} addCheckbox {}\n'.format(pname,pvalue,box,layout,addCheckbox))
        
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
        
        #add boxEdit layout
        layoutAttr=pname+'Layout'
        setattr(self,layoutAttr,self.addBoxEdit(pname,pvalue,layout,ledit,checkbox,elements=elements,disabledFlag=disabledFlag))
    
        
    def addBoxEdit(self,pname,pvalue,layout,ledit,checkbox,elements=None,disabledFlag=False):
        #setup boxEdit - boxEdit element values other than the list itself are not tracked and not saved in settings 
        if not hasattr(self,pname):
            setattr(self,pname,None); 
        boxEdit=DragAndDropList(self)
        boxEdit.setDragDropMode(QtWidgets.QAbstractItemView.InternalMove)
        boxEdit.setSelectionMode(QtWidgets.QAbstractItemView.ExtendedSelection)
        boxEdit.setStyleSheet(":disabled { color: #282828}")
        boxEdit.setMinimumHeight(60)
        
        #fill boxEdit - ONLY part that is tracked
        if hasattr(self,pname):
            entryList=getattr(self,pname)
            if type(entryList) is not list:
                boxEdit.addItems([entryList])
            else:
                boxEdit.addItems(entryList)
        boxEdit.setDisabled(disabledFlag)
        if elements:
            elements.append(boxEdit)
        #buttons
        browseBtn=None
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
        self.bgui.addList(pname,elements,enableCallback=lambda value, clearLedit: self.enableTextBox(value,clearLedit, checkbox,browseBtn,boxEdit,ledit,addBtn,removeBtn),updateCallback=lambda: self.updateTextBox(pname,ledit,boxEdit))
        if browseBtn:
            browseBtn.setIcon(self.browseIcon)
        addBtn.setIcon(self.addIcon)
        removeBtn.setIcon(self.removeIcon)
        
        #check rules for buttons    
        
        self.drawTextBoxBtnRules(boxEdit,ledit,addBtn,removeBtn)
        
        #connects from ledit and boxEdit to buttons
        ledit.textChanged.connect(lambda: addBtn.setEnabled(bool(ledit.text())))
        boxEdit.itemSelectionChanged.connect(lambda: removeBtn.setEnabled(bool(boxEdit.selectedItems())))
        
        #connect to changes in drag and drop
        boxEdit.itemMoved.connect(lambda : self.updateBoxEditValue(pname,boxEdit))
        
        #layout section
        filesBoxLeditLayout=QtGui.QVBoxLayout()
        #add to the main parameters box
        myBox=gui.vBox(None)
        myBox.layout().addLayout(filesBoxLeditLayout)
        startCol=0
        label=QtGui.QLabel(pvalue['label']+':')
        label.setAlignment(Qt.AlignTop)
        layout.addWidget(label)
        layout.addWidget(myBox,width=2,linefeed=1)
        #line layout     
        lineLayout=BwbGridLayout()
        
        if checkbox:
            lineLayout.addWidget(checkbox)
        lineLayout.addWidget(ledit)
        if browseBtn:
            lineLayout.addWidget(browseBtn)
        lineLayout.addWidget(addBtn)
        lineLayout.addWidget(removeBtn)
        
        filesBoxLeditLayout.addWidget(boxEdit)
        filesBoxLeditLayout.addLayout(lineLayout.layout())
        
        #just in case the state has changed while drawing this - do an update
        self.updateBoxEditValue(pname,boxEdit)
        return filesBoxLeditLayout
        
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
        self.btnRun = gui.button(None, self, "Start", callback=self.onRunClicked)
        self.btnRun.setStyleSheet(css)
        self.btnRun.setFixedSize(60,20)
        self.btnStop = gui.button(None, self, "Stop", callback=self.onStopClicked)
        self.btnStop.setStyleSheet(css)
        self.btnStop.setFixedSize(60,20)
        self.btnStop.setEnabled(False)
        self.execLayout.addWidget(self.btnRun,1,0)
        self.execLayout.addWidget(self.btnStop,1,1)
        self.execLayout.addWidget(myLabel,1,2)
        self.execLayout.addWidget(self.cboRunMode,1,3)
        if self.candidateTriggers:
            self.execLayout.addWidget(self.execBtn,1,4)
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
            self.triggerReady[attr]=True
        if not checked and attr in self.runTriggers:
            self.runTriggers.pop(attr)
            self.triggerReady.pop(attr)
    
    def checkTrigger(self):
        #this should be checked any time there is a change
        if self.runMode ==0: #manual - only go on start button
            return
        elif self.runMode ==1: #automatic same as pushing start button
            self.onRunClicked()
        elif self.runTriggers:            
            #check if the input triggers are set
            for trigger in self.runTriggers:
                if trigger not in self.triggerReady:
                    return
                if not self.triggerReady[trigger] or not self.inputConnections.isSet(trigger):
                    return
            self.onRunClicked()

    def bwbFileEntry(self, widget, button, ledit, icon=browseIcon,layout=None, label=None,entryType='file', checkbox=None):
        button.setStyleSheet("border: 1px solid #1a8ac6; border-radius: 2px;")
        button.setIcon(icon)
        ledit.setClearButtonEnabled(True)
        ledit.setPlaceholderText("Enter {}".format(entryType))
        if checkbox:
            layout.addWidget(checkbox)
        if label:
            #myLabel=QtGui.QLabel(label)
            #myLabel.setAlignment(Qt.AlignLeft | Qt.AlignVCenter)
            layout.addWidget(QtGui.QLabel(label))
        layout.addWidget(ledit)
        layout.addWidget(button)
        if not hasattr (layout,'added') or not getattr(layout,'added'):
            widget.layout().addLayout(layout.layout())
            setattr(layout,'added',True)
        layout.addLinefeed()

    def bwbLedit(self, widget,checkbox ,ledit, layout=None, label=None):
        ledit.setClearButtonEnabled(True)
        ledit.setPlaceholderText("Enter parameter")
        if(checkbox):
            layout.addWidget(checkbox)
            layout.addWidget(QtGui.QLabel(label))
        else:
            layout.addWidget(QtGui.QLabel(label),space=1)
        layout.addWidget(ledit)
        if not hasattr (layout,'added') or not getattr(layout,'added'):
            widget.layout().addLayout(layout.layout())
            setattr(layout,'added',True)
        layout.addLinefeed()
   
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
    
    def enableSpin(self,value,clearLedit,cb,spin):        
        if cb is None:
            spin.setEnabled(True)
        else:
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
        sys.stderr.write('checking inputs handler: attr {} value {}\n'.format(attr,value))
        if value is None:
            self.inputConnections.remove(attr,sourceId)
            sys.stderr.write('sig handler removing {} disabled {}\n'.format(attr,self.inputConnections.isConnected(attr)))
            setattr(self,attr,None) #this gives text None in ledit for some reason
            if attr in self.runTriggers:
                self.triggerReady[attr]=False
        else:
            self.inputConnections.add(attr,sourceId)
            sys.stderr.write('sig handler adding input: attr {} value {}\n'.format(attr,value,))
            setattr(self,attr,value)
            self.checkTrigger()
            if attr in self.runTriggers:
                sys.stderr.write("trigger value for attr {} is {}\n".format(attr,self.triggerReady[attr])) 
                self.triggerReady[attr]=True
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
            self.console.append("missing required parameters: {}".format(missingParms))
            return
        missingVols=self.getRequiredVols()
        if missingVols:
            self.console.append("missing or incorrect volume mappings to: {}".format(missingVols))
            return
        #get ready to start
        attrList=self.__dict__.keys()
        self.bgui.disableAll()
        self.disableExec()
        cmd=self.generateCmdFromData()
        self.envVars={}
        self.getEnvironmentVariables()
        try:
            imageName='{}:{}'.format(self._dockerImageName, self._dockerImageTag)
            self.pConsole.writeMessage('Generating Docker command from image {}\nVolumes {}\nCommands {}\nEnvironment {}\n'.format(imageName, self.hostVolumes, cmd , self.envVars))
            self.status='running'
            self.setStatusMessage('Running...')
            self.dockerClient.create_container_cli(imageName, hostVolumes=self.hostVolumes, commands=cmd, environment=self.envVars,consoleProc=self.pConsole)
        except BaseException as e:
            self.bgui.reenableAll(self)
            self.reenableExec()
            self.pConsole.writeMessage("unable to start Docker command "+ str(e))
            self.setStatusMessage('Error')

    def checkRequiredParms(self):
        for parm in self.data['requiredParameters']:
            if hasattr(self,parm):
                if not getattr (self,parm):
                    return parm
        return None

    def getRequiredVols(self):
        #get all the autoMaps
        #the mountpoint is passed because it will be converted later into the global hostpath
        bwbDict={}
        if 'autoMap' in self.data and self.data['autoMap']:
            for bwbVolume, containerVolume in self.dockerClient.bwbMounts.items():
                self.hostVolumes[containerVolume]=bwbVolume
                bwbDict[containerVolume]=bwbVolume
        if 'volumeMappings' in self.data and self.data['volumeMappings']:
            for mapping in self.data['volumeMappings']:
                conVol=mapping['conVolume']
                attr=mapping['attr']
                if not hasattr(self,attr) or not getattr(self,attr):
                    if 'autoMap' in self.data and self.data['autoMap'] and conVol in bwbDict:
                        setattr(self,attr,bwbDict[conVol])
                    else:
                        return conVol
                self.hostVolumes[conVol]=getattr(self,attr)
        return None
    
    def joinFlagValue(self,flag,value):
        if flag:
            if flag.strip()[-1] == '=':
                return flag.strip()+value.strip()
            return flag.strip()+ ' ' + str(value).strip()
        return value
    
    def flagString(self, pname):
        if 'parameters' in self.data and pname in self.data['parameters']:
            pvalue= self.data['parameters'][pname]
            sys.stderr.write('pname {} pvalue {}\n'.format(pname,pvalue))
            if('flag' in pvalue  and hasattr(self,pname)):
                flagName=pvalue['flag']
                flagValue=getattr(self,pname)
                if pvalue['type'] == 'bool':
                    if flagValue:
                        return flagName
                    return None
                elif pvalue['type'] == 'file':
                    filename=str(flagValue)
                    if filename:
                        hostFilename=self.bwbPathToContainerPath(filename, isFile=True,returnNone=False)
                        sys.stderr.write('orig file {} convert to container {}\n'.format(filename,hostFilename))
                        return self.joinFlagValue(flagName,hostFilename)
                    return None
                elif pvalue['type'] == 'file list':
                    #files=str.splitlines(flagValue)
                    files=flagValue
                    sys.stderr.write('flagnName {} files are {}\n'.format(flagName,files))
                    if files:
                        hostFiles=[]
                        for f in files:
                            hostFiles.append(self.bwbPathToContainerPath(f, isFile=True,returnNone=False))
                            sys.stderr.write('flagnName {} adding file {} hostFiles {}\n'.format(flagName,f,hostFiles))
                        if flagName:
                            return " ".join([flagName] + hostFiles)
                        else:
                            return " ".join(hostFiles)
                    return None
                elif pvalue['type'] == 'directory':
                    path=str(flagValue)
                    if path:
                        hostPath=self.bwbPathToContainerPath(path, returnNone=False)
                    return self.joinFlagValue(flagName,str(hostPath))
                elif pvalue['type'][-4:] =='list':
                    if flagName:
                        return flagName +" " + " ".join(flagValue)
                    else:
                        return " ".join(flagValue)
                elif flagValue:
                    return self.joinFlagValue(flagName,flagValue)
            
        return None
            
    def generateCmdFromData (self):
        flags=[]
        args=[]
        for pname, pvalue in self.data['parameters'].items():
            #possible to have an requirement or parameter that is not in the executable line
            #environment variables can have a Null value in the flags field
            #arguments are the only type that have no flag
            if 'argument' in pvalue:
                fStr=self.flagString(pname)
                args.append(fStr)
                continue
            if 'flag' not in pvalue or pvalue['flag'] is None:
                continue
            #if required or checked then it is added to the flags
            addParms=False
            #checkattr is needed for the orange gui checkboxes but is not otherwise updated 
            if pname in self.optionsChecked and self.optionsChecked[pname]:
                addParms=True
                
            #also need to check for booleans which are not tracked by optionsChecked
            if pvalue['type'] == 'bool'and hasattr(self,pname) and getattr(self,pname):
                addParms=True
                
            if pname in self.data['requiredParameters'] and hasattr(self,pname):
                addParms=True

            if addParms:
                fStr=self.flagString(pname)
                if pvalue['flag']:
                    if fStr:
                        flags.append(fStr)
                        
        return self.generateCmdFromBash(self.data['command'],flags=flags,args=args)

    def replaceVars(self,cmd,pnames,varSeen):
        pattern = r'\_bwb\{([^\}]+)\}'
        regex = re.compile(pattern)
        subs=[]
        sys.stderr.write('command is {}\n'.format(cmd))
        for match in regex.finditer(cmd):
            sys.stderr.write('matched {}\n'.format(match.group(1)))
            sub=match.group(1)
            if sub not in subs:
                subs.append(sub)
        for sub in subs:
            #remove front and trailing spaces
            #create match
            matchStr='_bwb{' + sub +'}'
            psub=sub.strip()
            sys.stderr.write('looking for {} in {}\n'.format(psub,pnames))
            if psub in pnames:
                fStr=self.flagString(psub)
                sys.stderr.write('psub {} is  in {}\n'.format(psub,pnames))
                sys.stderr.write('fstr {}\n'.format(fStr))
                if fStr:
                    sys.stderr.write('replace matchStr {} in cmd {} with fStr {}\n'.format(matchStr,cmd,fStr))
                    cmd=cmd.replace(matchStr,fStr)
                    sys.stderr.write('to give new cmd {}\n'.format(cmd))
                    varSeen[fStr]=True
                else:
                    cmd=cmd.replace(matchStr,'')
            else:
                cmd=cmd.replace(matchStr,'')
        return cmd
        
    def generateCmdFromBash (self, executables, flags = [], args = []):
        #by default all flags and arguments are applied to final command
        #use positional _bwb{} variables so specify flags and arguments if there are multiple commands
        #unused args and flags are applied to the final command
        
        #multi executable commands need to use bash -c 'cmd1 && cmd2' type syntax - note this can cause problems when stopping container
        cmdStr="bash -c '"
        if len(executables) == 1:
            cmdStr=""
        varSeen={}
        sys.stderr.write('execs {}\n'.format(executables))
        sys.stderr.write('executables {}\n'.format(executables))
        for executable in executables:
            lastExecutable= (executable == executables[-1])
            sys.stderr.write('Orig executable {} flags {} args {}\n'.format(executable,flags,args))
            pnames=None
            if 'parameters' in self.data:
                pnames=self.data['parameters'].keys()
            executable=self.replaceVars(executable,pnames,varSeen)
            sys.stderr.write('New executable {} flags {} args {}\n'.format(executable,flags,args))
            cmdStr += executable +' '
            #extra args and flags go after last executable
            if lastExecutable:
                for flag in flags:
                    if flag not in varSeen:
                        cmdStr+=flag + ' '
                for arg in args:
                    if arg not in varSeen:
                        cmdStr+= str(arg) + ' '
                if len(executables) >1:
                    cmdStr+="'"
            else:
                cmdStr+=" && "
        return cmdStr

    def getEnvironmentVariables(self):
        #dynamic environment variables
        for pname in self.data['parameters']:
            pvalue=self.data['parameters'][pname]
            if 'env' in pvalue and getattr(self,pname) is not None:
                checkAttr=pname+'Checked'
                setenv=False
                if pname in self.optionsChecked and self.optionsChecked[pname]:
                    self.envVars[pvalue['env']] = getattr(self,pname)
                    sys.stderr.write('optional var {} env {} assigned to {}\n'.format(pname,pvalue['env'],getattr(self,pname)))
                elif hasattr(self,checkAttr):
                    if getattr(self,checkAttr):
                        self.envVars[pvalue['env']] = getattr(self,pname)
                        sys.stderr.write('optional var {} env {} assigned to {}\n'.format(pname,pvalue['env'],getattr(self,pname)))
                else:
                    self.envVars[pvalue['env']] = getattr(self,pname)
                    sys.stderr.write('var {} env {} assigned to {}\n'.format(pname,pvalue['env'],getattr(self,pname)))
                    
        #now assign static environment variables
        if 'env' in self.data:
            for e in self.data['env']:
                if e not in self.envVars:
                    self.envVars[e]=self.data['env'][e]

#Event handlers
    def onRunClicked(self):
        if hasattr (self,'userStartJob'):
            self.userStartJob()
        else:
            self.startJob()
            
    def onStopClicked(self):
        self.pConsole.stop('Stopped by user')
        self.setStatusMessage('Stopped')
        self.status='stopped'
        self.bgui.reenableAll(self)
        self.reenableExec()
        
    def onRunFinished(self,code=None,status=None):
        self.pConsole.writeMessage("Finished")
        if code is not None:
           self.pConsole.writeMessage("Exit code is {}".format(code))
        if status is not None:
           self.pConsole.writeMessage("Exit status is {}".format(status))
        self.bgui.reenableAll(self)
        self.reenableExec()
        if self.status != 'stopped' and self.status != 'finished':
            if code or status:
                self.setStatusMessage("Error")
                self.status='error'
            else:
                self.setStatusMessage('Finished')
                self.satus='finished'
                self.handleOutputs()
    
    def onRunError(self,error):
        self.bgui.reenableAll(self)
        self.reenableExec()
        self.console.writeMessage("Error occurred {}\n".format(error),color=Qt.red)
        
    def onRunMessage(self, message):
        self.pConsole.writeMessage(message)

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
            sys.stderr.write('dirPath {} pathFile {} hostPath {}\n'.format(dirPath,pathFile,hostPath))
        else:
            hostPath=os.path.normpath(self.dockerClient.to_best_host_directory(path,returnNone=False))
            sys.stderr.write('hostPath {}\n'.format(hostPath))
        conPath=None
        #now get all the possible submappings to volumeMappings by comparing the true hostmappings
        #if submap is found convert the common path to the container path
        #return shortest path
        for conVol, bwbVol in self.hostVolumes.items():
            hostVol=os.path.normpath(self.dockerClient.to_best_host_directory(bwbVol,returnNone=False))
            sys.stderr.write('checking conVol {} bwbVol {} hostVol {}\ncommon {}\n'.format(conVol,bwbVol,hostVol,os.path.commonpath([hostVol,hostPath])))
            prefix=None
            if hostVol == hostPath:
                prefix=""
            elif Path(hostVol) in Path(hostPath).parents:
                prefix=str(Path(hostPath).relative_to(hostVol))
                sys.stderr.write('prefix is {}\n'.format(prefix))
                if prefix == '.':
                    prefix=''
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
        self.btnStop.setEnabled(True)
        self.cboRunMode.setEnabled(False)
            
    def reenableExec(self):
        self.btnRun.setText('Start')
        self.btnRun.setEnabled(True)
        self.btnStop.setEnabled(False)
        self.cboRunMode.setEnabled(True)
        

