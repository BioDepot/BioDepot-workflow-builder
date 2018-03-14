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

class NamedDict():
    def __init__(self):
        name=None
        _dict={}
        
    def setName(self,inputName):
        self.name=inPutname
        
    def setElement(self,key,value):
        self._dict[key]=value
        
    def delElement(self,key):
        del self._dict[key]
    
    def printItem(self):
        outputStr="{}, {{ ".format(name)
        for key, value in sorted(self._dict.items()):
            outputStr+="{} : {} ".format(key,value)
        outputStr+= '}}'
        return outputStr

class OWWidgetBuilder(widget.OWWidget):
    name = "WidgetBuilder"
    description = "Set one or more data files"
    category = "Data"
    icon = "icons/build.png"
    priority = 2

    inputs = []
    outputs = []

    want_main_area = False
    want_control_area = True
    filesList = settings.Setting([], schema_only=True)
    inputConnections = settings.Setting(None,schema_only=True)
    checked=settings.Setting(False,schema_only=True)

    def __init__(self):
        super().__init__()
        css = '''
        QPushButton {background-color: #1588c5; color: white; height: 20px; border: 1px solid #1a8ac6; border-radius: 2px;}
        QPushButton:pressed { background-color: #158805; border-style: inset;}
        QPushButton:disabled { background-color: lightGray; border: 1px solid gray; }
        QPushButton:hover {background-color: #1588f5; }
        '''  
        self.setStyleSheet(css)
        self.browseIcon=QtGui.QIcon('/biodepot/orangebiodepot/icons/bluefile.png')
        self.addIcon=QtGui.QIcon('/biodepot/orangebiodepot/icons/add.png')
        self.removeIcon=QtGui.QIcon('/biodepot/orangebiodepot/icons/remove.png')
        self.submitIcon=QtGui.QIcon('/biodepot/orangebiodepot/icons/submit.png')
        self.reloadIcon=QtGui.QIcon('/biodepot/orangebiodepot/icons/reload.png')
        

        #set up basic gui 
        self.setStyleSheet(":disabled { color: #282828}")
        self.scroll_area = QScrollArea(
            verticalScrollBarPolicy=Qt.ScrollBarAlwaysOn
        )
        self.bigBox=gui.widgetBox(self.controlArea)
        self.scroll_area.setWidget(self.bigBox)
        self.scroll_area.setWidgetResizable(True)
        self.controlArea.layout().addWidget(self.scroll_area)
        controlBox = gui.vBox(self.bigBox)
        #self.drawExec(box=controlBox)ne)
        requiredBox = gui.widgetBox(self.bigBox, "Required entries")
        #draw Ledits for the frequired elements
        leditRequiredLayout=QtGui.QGridLayout()
        leditRequiredLayout.setSpacing(5)
        setattr(leditRequiredLayout,'nextRow',1)        
        requiredBox.layout().addLayout(leditRequiredLayout)
        for pname in ['Name','Description','Category','Docker_image_name','Docker_image_tag']:
            self.drawLedit(pname,requiredBox,layout=leditRequiredLayout)
        self.Priority=10
        self.drawLedit('Priority',requiredBox,layout=leditRequiredLayout)
        #file entry for icon
        self.drawFileLedit('Icon',requiredBox,layout=leditRequiredLayout)
        #listwidgets for inputs and outputs 
        self.drawListWidget('Inputs',requiredBox,layout=leditRequiredLayout)
        self.drawListWidget('Outputs',requiredBox,layout=leditRequiredLayout)
        self.drawVolumeListWidget('Volumes',requiredBox,layout=leditRequiredLayout)
        self.drawParamsListWidget('Parameters',requiredBox,layout=leditRequiredLayout)
        self.drawCommand('Commands',requiredBox,layout=leditRequiredLayout)
    
    #logic resides here
    
    def updateCheckBox(self,checkBox,widget=None):
        if(checkBox.isEnabled()):
            widget.setEnabled(checkBox.isChecked())
        
    def addListWidget(self,pname,boxEdit,comboBox,addBtn):
        pass
    def removeListWidget(self,pname,boxEdit,comboBox,addBtn):
        pass 
    
    def drawCommand(self,pname,box=None,layout=None):
        label=QtGui.QLabel(pname+":")
        label.setAlignment(Qt.AlignTop)
        textBox=QtGui.QPlainTextEdit()
        layout.addWidget(label,layout.nextRow,0)
        layout.addWidget(textBox,layout.nextRow,1,1,4)
        layout.nextRow = layout.nextRow + 1
    
    #workhorse widget
    def makeLedit(self,leditAttr,text=None,label=None,addCheckBox=False,addBrowse=False,initialValue=None):
        leditLabel=None
        checkBox=None
        if(label):
            leditLabel=QtGui.QLabel(label)
        if not hasattr(self,leditAttr):
            setattr(self,leditAttr,initialValue);
        ledit = gui.lineEdit(None, self,leditAttr)
        ledit.setClearButtonEnabled(True)
        ledit.setPlaceholderText(text)
        ledit.setStyleSheet(":disabled { color: #282828}")
        if not initialValue:
            ledit.clear()
        if addCheckBox:
            checkAttr=leditAttr+'Checked'
            if not hasattr(self,checkAttr):
                setattr(self,checkAttr,False)
            checkBox=gui.checkBox(None, self,checkAttr,label=None)
            checkBox.setChecked(getattr(self,checkAttr))
            checkBox.stateChanged.connect(lambda : self.updateCheckBox(checkBox,ledit))
            ledit.setEnabled(checkBox.isChecked())
            
        box=QHBoxLayout()
        if checkBox:
            box.addWidget(checkBox)
        if leditLabel:
            box.addWidget(leditLabel)
            box.addWidget(ledit)
        return box
        
    def makeComboBox (self,pname,label, elements):
        comboBoxLabel=QtGui.QLabel(label)
        comboBox=QtGui.QComboBox()
        comboBox.addItems(elements)
        comboBox.currentIndex=0
        box=QHBoxLayout()
        box.addWidget(comboBoxLabel)
        box.addWidget(comboBox)
        return box
        
        
    def makeListWidget(self,pname,boxEdit):
        #setup boxEdit
        #logic is handled by add remove buttons
        if not hasattr(self,pname):
            setattr(self,pname,None); 
        boxEdit=QtGui.QListWidget(self)
        boxEdit.setDragDropMode(QtWidgets.QAbstractItemView.InternalMove)
        boxEdit.setSelectionMode(QtWidgets.QAbstractItemView.ExtendedSelection)
        boxEdit.setStyleSheet(":disabled { color: #282828}")
        boxEdit.setMinimumHeight(60)
        return boxEdit
        
    #multiple widget elements       
    def drawListWidget (self, pname, box=None, layout=None):
        
        #save items to qlist for ordered Dict
        #each item in the list is an attribute and a dict
            
        #setup boxEdit
        boxEdit=self.makeListWidget(pname,None)

        #ledit
        checkBox=None
        nameBox=self.makeLedit(pname+'nameLedit','Enter name','Name')
        callbackBox=self.makeLedit(pname+'callbackLedit','Enter callback', 'Callback',addCheckBox=True)

        #buttons
        addBtn=gui.button(None, self, "", callback=lambda: self.addListWidget(pname,boxEdit,comboBox,addBtn), autoDefault=False)
        removeBtn=gui.button(None, self, "", callback=lambda: self.removeListWidget(pname,boxEdit,removeBtn), autoDefault=False)
        addBtn.setIcon(self.addIcon)
        removeBtn.setIcon(self.removeIcon)
        
        #comboBox
        comboBox=self.makeComboBox(pname,'Type:',['str'])
        
        #set styles
        buttonStyle='background: None; border: None ; border-radius: 0;'
        addBtn.setStyleSheet(buttonStyle)
        removeBtn.setStyleSheet(buttonStyle)
        
        #layout section
        filesBoxLeditLayout=QtGui.QVBoxLayout()
        
        #add to the main parameters box
        myBox=gui.vBox(None)
        myBox.layout().addLayout(filesBoxLeditLayout)
        startCol=0
        label=QtGui.QLabel(pname+':')
        label.setAlignment(Qt.AlignTop)
        layout.addWidget(label,layout.nextRow,startCol)
        layout.addWidget(myBox,layout.nextRow,1,1,2)
        layout.nextRow = layout.nextRow + 1
        
        #line layout     
        lineLayout=QtGui.QGridLayout()
        lineLayout.addLayout(nameBox,1,0)
        lineLayout.addLayout(callbackBox,1,1)
        lineLayout.addLayout(comboBox,1,2)
        lineLayout.addWidget(addBtn,1,3)
        lineLayout.addWidget(removeBtn,1,4)
        
        #now add the two layouts to the bigBox layout
        filesBoxLeditLayout.addWidget(boxEdit)
        filesBoxLeditLayout.addLayout(lineLayout)

    def drawVolumeListWidget (self, pname, box=None, layout=None):
        
        #save items to qlist for ordered Dict
        #each item in the list is an attribute and a dict
            
        #setup boxEdit
        boxEdit=self.makeListWidget(pname,None)
        
        #ledit and buttons
        volumeBox=self.makeLedit(pname+'volumenLedit','Enter additional volume','Additional volume')
        
     
        #buttons
        addBtn=gui.button(None, self, "", callback=lambda: self.addListWidget(pname,boxEdit,comboBox,addBtn), autoDefault=False)
        removeBtn=gui.button(None, self, "", callback=lambda: self.removeListWidget(pname,boxEdit,removeBtn), autoDefault=False)
        addBtn.setIcon(self.addIcon)
        removeBtn.setIcon(self.removeIcon)
        
        #checkBox
        if not hasattr(self,'AutoMap'):
            setattr(self,'AutoMap',True)
        checkBox=gui.checkBox(box, self,'AutoMap', 'AutoMap')

        #set styles
        buttonStyle='background: None; border: None ; border-radius: 0;'
        addBtn.setStyleSheet(buttonStyle)
        removeBtn.setStyleSheet(buttonStyle)
        
        #layout section
        filesBoxLeditLayout=QtGui.QVBoxLayout()
        
        #add to the main parameters box
        myBox=gui.vBox(None)
        myBox.layout().addLayout(filesBoxLeditLayout)
        startCol=0
        label=QtGui.QLabel(pname+':')
        label.setAlignment(Qt.AlignTop)
        layout.addWidget(label,layout.nextRow,startCol)
        layout.addWidget(myBox,layout.nextRow,1,1,2)
        layout.nextRow = layout.nextRow + 1
        
        #line layout     
        lineLayout=QtGui.QGridLayout()
        lineLayout.addWidget(checkBox,1,0)
        lineLayout.addLayout(volumeBox,1,1)
        lineLayout.addWidget(addBtn,1,2)
        lineLayout.addWidget(removeBtn,1,3)
        
        #now add the two layouts to the bigBox layout
        filesBoxLeditLayout.addWidget(boxEdit)
        filesBoxLeditLayout.addLayout(lineLayout)
        
    def drawParamsListWidget (self, pname, box=None, layout=None):
        #save items to qlist for ordered Dict
        #each item in the list is an attribute and a dict
            
        #setup boxEdit
        boxEdit=self.makeListWidget(pname,None)

        #ledit
        checkBox=None
        nameBox=self.makeLedit(pname+'nameLedit','Enter name','Name')
        flagBox=self.makeLedit(pname+'flagLedit','Enter flag', 'flag',addCheckBox=True)
        labelBox=self.makeLedit(pname+'labelLedit','Enter label', 'label',addCheckBox=True)
        envBox=self.makeLedit(pname+'envLedit','Enter ENV variable', 'env',addCheckBox=True)
        defaultBox=self.makeLedit(pname+'defaultLedit','Enter default', 'default',addCheckBox=True)
        
        #buttons
        addBtn=gui.button(None, self, "", callback=lambda: self.addListWidget(pname,boxEdit,comboBox,addBtn), autoDefault=False)
        removeBtn=gui.button(None, self, "", callback=lambda: self.removeListWidget(pname,boxEdit,removeBtn), autoDefault=False)
        addBtn.setIcon(self.addIcon)
        removeBtn.setIcon(self.removeIcon)
        
        #checkBox

        setattr(self,'optional',False)
        optionalCb=gui.checkBox(box, self,'optional', 'Optional')
        
        setattr(self,'argument',False)
        argumentCb=gui.checkBox(box, self,'argument', 'Argument')
        #connect argument checkbox to disabling the flag checkbox       
        argumentCb.stateChanged.connect(lambda : self.updateCheckBoxLayout(argumentCb,flagBox))
        
        
        #comboBox
        comboBox=self.makeComboBox(pname,'Type:',['str','file','directory','directoryList','fileList','int','double'])
        
        #set styles
        buttonStyle='background: None; border: None ; border-radius: 0;'
        addBtn.setStyleSheet(buttonStyle)
        removeBtn.setStyleSheet(buttonStyle)
        
        #layout section
        filesBoxLeditLayout=QtGui.QVBoxLayout()
        
        #add to the main parameters box
        myBox=gui.vBox(None)
        myBox.layout().addLayout(filesBoxLeditLayout)
        startCol=0
        label=QtGui.QLabel(pname+':')
        label.setAlignment(Qt.AlignTop)
        layout.addWidget(label,layout.nextRow,startCol)
        layout.addWidget(myBox,layout.nextRow,1,1,2)
        layout.nextRow = layout.nextRow + 1
        
        #line layout     
        lineLayout=QtGui.QGridLayout()
        lineLayout.addLayout(nameBox,1,0)
        lineLayout.addLayout(comboBox,1,1)
        lineLayout.addLayout(flagBox,1,2)
        lineLayout.addWidget(argumentCb,1,3)
        lineLayout.addLayout(envBox,2,0)
        lineLayout.addLayout(labelBox,2,1)
        lineLayout.addLayout(defaultBox,2,2)
        lineLayout.addWidget(optionalCb,2,3)        
        lineLayout.addWidget(addBtn,2,5)
        lineLayout.addWidget(removeBtn,2,6)
        
        #now add the two layouts to the bigBox layout
        filesBoxLeditLayout.addWidget(boxEdit)
        filesBoxLeditLayout.addLayout(lineLayout)
                    
    def drawFileLedit(self,pname, box=None, layout=None):
        if not hasattr(self,pname):
            setattr(self,pname,None)
        ledit=gui.lineEdit(None, self, pname)
        if getattr(self,pname) is None:
            ledit.clear()
        #note that using lambda does not work in this function - truncates string variable so partial used instead
        button=gui.button(None, self, "", callback= partial(self.browseFileDir, attr= pname ),
                          autoDefault=True, width=19, height=19)
        button.setIcon(self.browseIcon)
        layout.addWidget(QtGui.QLabel(pname+":"),layout.nextRow,0)
        ledit.setPlaceholderText("Enter {}".format(pname))
        ledit.setClearButtonEnabled(True)
        layout.addWidget(ledit,layout.nextRow,1)
        layout.addWidget(button,layout.nextRow,2)
        layout.nextRow+=1        

    def drawLedit(self,pname,box=None,layout=None):
        if not hasattr(self,pname):
            setattr(self,pname,None)
        ledit=gui.lineEdit(None, self, pname)
        if getattr(self,pname) is None:
            ledit.clear()
        layout.addWidget(QtGui.QLabel(pname+":"),layout.nextRow,0)
        ledit.setPlaceholderText("Enter {}".format(pname))
        ledit.setClearButtonEnabled(True)
        layout.addWidget(ledit,layout.nextRow,1,1,2)
        layout.nextRow+=1
     
    def browseFileDir(self, attr):
        defaultDir = '/root'
        if os.path.exists('/data'):
            defaultDir = '/data'
        myFile=QtWidgets.QFileDialog.getOpenFileName(self, "Locate file", defaultDir)[0]
        if myFile:
            setattr(self,attr,myFile)

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
