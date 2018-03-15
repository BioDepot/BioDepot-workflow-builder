import os
import re
import sys
from collections import OrderedDict
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
        
    def getState(self):
        states={}
        for key, element in self.guiElements.items():
           states[key]=element.getState()
        return states
        
    def loadState(self,states={}):
        self.states=states
        for key, element in self.guiElements.items():
            element.setState(self.states[key])
            
    def blankState(self):
        for key, element in self.guiElements.items():
            element.setState(None)
            self.states[key]=None
        
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
        #checkBox for automap Volumens

        #auto
        
        #listwidgets for inputs and outputs 
        #define widgetItems for the different widgetLists
        
        self.drawIOListWidget('Inputs',requiredBox,layout=leditRequiredLayout)
        self.drawIOListWidget('Outputs',requiredBox,layout=leditRequiredLayout)
        self.drawCheckBox('AutoMap',requiredBox,layout=leditRequiredLayout,default=True, label='Pass current Bwb volumes to container')
        self.drawVolumeListWidget('Volumes',requiredBox,layout=leditRequiredLayout)
        self.drawParamsListWidget('Parameters',requiredBox,layout=leditRequiredLayout)
        self.drawCommand('Commands',requiredBox,layout=leditRequiredLayout)
    
    #logic resides here
    
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
            qWidgetItem.blankState()
        
    def removeListWidget(self,removeBtn,qWidgetList,qWidgetItem):
        if qWidgetList.selectedItems():
            for item in qWidgetList.selectedItems():
                del qWidgetList.states[qWidgetList.row(item)]
                qWidgetList.takeItem(qWidgetList.row(item))
        if not qWidgetList.count():
            removeBtn.setEnabled(False)
        qWidgetItem.blankState()
    
    def onListWidgetSelect(self,qWidgetList,addBtn,removeBtn,qWidgetItem):
        if qWidgetList.selectedItems():
            if len(qWidgetList.selectedItems()) > 1:
                addBtn.setEnabled(False)
            else:
                item=qWidgetList.selectedItems()[0]
                myState=qWidgetList.states[qWidgetList.row(item)]
                qWidgetItem.loadState(myState)
                addBtn.setEnabled(True)
        else:
            removeBtn.setEnabled(False)
            addBtn.setEnabled(True)
    
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
        setattr(box,'getValue',lambda : self.getLeditValue(ledit))
        setattr(box,'getState',lambda : self.getLeditState(checkBox,leditLabel,ledit))
        setattr(box,'setState',lambda state : self.setLeditState(checkBox,leditLabel,ledit,state))
        return box
        
    def getLeditValue(self,ledit):
        if ledit.isEnabled():
            return ledit.text()
        return None
        
    def getLeditState(self,checkBox,label,ledit):
        state={}
        if checkBox:
            state[checkBox]=[checkBox.isEnabled() ,checkBox.isChecked()]
        else:
            state[checkBox]=None
        if label:
            state[label]=[label.isEnabled(), label.text() ]
        else:
            state[label]=None
        if ledit:
            state[ledit]=[ledit.isEnabled(), ledit.text() ]
        else:
            state[label]=None
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
            
        if state[checkBox] is None:
            checkBox=None
        else:
            checkBox.setEnabled(state[checkBox][0])
            checkBox.setChecked(state[checkBox][1])
        if state[ledit] is None:
            ledit=None
        else:
            ledit.setEnabled(state[ledit][0])
            if state[ledit][1]:
                ledit.setText(state[ledit][1])
            else:
                ledit.clear()
        if state[label] is None:
            label=None
        else:
            label.setEnabled(state[label][0])
            label.setText(state[label][1])
                
       
    def makeComboBox (self,pname,label, elements):
        comboBoxLabel=QtGui.QLabel(label)
        comboBox=QtGui.QComboBox()
        comboBox.addItems(elements)
        comboBox.currentIndex=0
        box=QHBoxLayout()
        box.addWidget(comboBoxLabel)
        box.addWidget(comboBox)
        setattr(box,'getValue',lambda : self.getComboValue(comboBox))
        setattr(box,'getState',lambda : self.getComboState(comboBox,comboBoxLabel))
        setattr(box,'setState',lambda state : self.setComboState(comboBox,comboBoxLabel,state))
        return box
        
    def getComboValue(self,comboBox):
        if comboBox.isEnabled():
            return str(comboBox.currentText())
        return None 
        
    def getComboState(self,comboBox,label):
        state={}
        if comboBox:
            allItems=[comboBox.itemText(i) for i in range (comboBox.count())]
            enabled=comboBox.isEnabled()
            index=comboBox.currentIndex
            state[comboBox]=[enabled,index,allItems]
        else:
            state[comboBox]=None
        if label:
            state[label]=[label.isEnabled(),label.text()]
        else:
            state[label]=None
        return state 
    
    def setComboState(self,comboBox,label,state):
        if state is None:
            #intialize
            comboBox.setEnabled(True)
            comboBox.setCurrentIndex(0)
        return
        if state[comboBox] is None:
            comboBox=None
        else:
           comboBox.setEnabled(state[comboBox][0])
           comboBox.setItems(state[comboBox][2])
           comboBox.setCurrentIndex(state[comboBox][1])
           
        if state[label] is None:
            label=None
        else:
            label.setEnabled(state[label][0])
            label.setText(state[label][1])
        
    def makeCheckBox (self, pname,label):
        #checkbox for binary options
        #not used as part of other elements
        setattr(self,'pname',False)
        checkBox=gui.checkBox(None, self, pname, label)
        setattr(checkBox,'getValue',lambda : checkBox.isEnabled() and checkBox.isChecked())
        setattr(checkBox,'getState',lambda : self.getCheckBoxState(checkBox,label))
        setattr(checkBox,'setState',lambda state : self.setCheckBoxState(checkBox,label,state))
        return checkBox
     
    def getCheckBoxState(self,checkBox,label):
        state={}
        if checkBox:
            state[checkBox]=[checkBox.isEnabled() ,checkBox.isChecked()]
        else:
            state[checkBox]=None
        return state
        
    def setCheckBoxState(self,checkBox,label,state):
        if state is None:
            #initialize
            checkBox.setEnabled(True)
            checkBox.setChecked(False)
            return
            
        if state[checkBox] is None:
            checkBox=None
        else:
            checkBox.setEnabled(state[checkBox][0])
            checkBox.setChecked(state[checkBox][1])   
    
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
    def drawIOListWidget (self, pname, box=None, layout=None, attrValues=None):
        
        #the structure for values element of the widget list is:
        # (name,{}) i.e. a pair of str and dict
        # we have a definition list to tell us which values are associated with which elements
        # onchange of the element, the value gets changed
        # This is an orderdict of elements and keyvalues
        # all of this is stored in a class for widgetItem
    
        
        #setup boxEdit
        boxEdit=self.makeListWidget(pname,None)
        setattr(boxEdit,'states',[])


        #ledit
        checkBox=None
        nameBox=self.makeLedit(pname+'nameLedit','Enter name','Name')
        callbackBox=self.makeLedit(pname+'callbackLedit','Enter callback', 'Callback',addCheckBox=True)
        
        #comboBox
        comboBox=self.makeComboBox(pname,'Type:',['str'])
        
        lineItem=WidgetItem(OrderedDict([('name',nameBox),('type',comboBox),('Callback',callbackBox)]))
        
        #buttons
        addBtn=gui.button(None, self, "", callback=lambda: self.addListWidget(addBtn,boxEdit,lineItem), autoDefault=False)
        removeBtn=gui.button(None, self, "", callback=lambda: self.removeListWidget(removeBtn,boxEdit,lineItem), autoDefault=False)
        removeBtn.setEnabled(bool(boxEdit.selectedItems()))
        boxEdit.itemSelectionChanged.connect(lambda: removeBtn.setEnabled(bool(boxEdit.selectedItems())))
        boxEdit.itemSelectionChanged.connect(lambda: self.onListWidgetSelect(boxEdit,addBtn,removeBtn,lineItem))
        
        addBtn.setIcon(self.addIcon)
        removeBtn.setIcon(self.removeIcon)
        
        #set up base item

        
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
            
        #setup boxEdit
        boxEdit=self.makeListWidget(pname,None)
        setattr(boxEdit,'states',[])
        
        #ledit and buttons
        nameBox=self.makeLedit(pname+'Name','Enter name','Name')
        volumeBox=self.makeLedit(pname+'volumeLedit','Enter volume','Additional volume')
     
        #define the line item being edited
        lineItem=WidgetItem(OrderedDict([('name',nameBox),('volume',volumeBox)]),{})

        
        #buttons
        addBtn=gui.button(None, self, "", callback=lambda: self.addListWidget(addBtn,boxEdit,lineItem), autoDefault=False)
        removeBtn=gui.button(None, self, "", callback=lambda: self.removeListWidget(removeBtn,boxEdit,lineItem), autoDefault=False)
        removeBtn.setEnabled(bool(boxEdit.selectedItems()))
        boxEdit.itemSelectionChanged.connect(lambda: removeBtn.setEnabled(bool(boxEdit.selectedItems())))
        boxEdit.itemSelectionChanged.connect(lambda: self.onListWidgetSelect(boxEdit,addBtn,removeBtn,lineItem))

        addBtn.setIcon(self.addIcon)
        removeBtn.setIcon(self.removeIcon)
        
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
        setattr(boxEdit,'states',[])
        
        #ledit
        checkBox=None
        nameBox=self.makeLedit(pname+'nameLedit','Enter name','Name')
        flagBox=self.makeLedit(pname+'flagLedit','Enter flag', 'flag',addCheckBox=True)
        labelBox=self.makeLedit(pname+'labelLedit','Enter label', 'label',addCheckBox=True)
        envBox=self.makeLedit(pname+'envLedit','Enter ENV variable', 'env',addCheckBox=True)
        defaultBox=self.makeLedit(pname+'defaultLedit','Enter default', 'default',addCheckBox=True)

        #checkBox

        setattr(self,'optional',False)
        setattr(self,'argument',False)
        optionalCb=self.makeCheckBox ('optional','Optional')
        argumentCb=self.makeCheckBox ('argument','Argument')
        #connect argument checkbox to disabling the flag checkbox       
        argumentCb.stateChanged.connect(lambda : self.updateCheckBoxLayout(argumentCb,flagBox))
        
        #comboBox
        comboBox=self.makeComboBox(pname,'Type:',['str','file','directory','directoryList','fileList','int','double'])        
        #define the line item being edited
        lineItem=WidgetItem(OrderedDict([('name',nameBox),('flag',flagBox), ('label',labelBox), ('env',envBox),
                                         ('default',defaultBox),('optional',optionalCb),('argument',argumentCb),
                                         ('type',comboBox)]),{})
        
        #buttons
        addBtn=gui.button(None, self, "", callback=lambda: self.addListWidget(addBtn,boxEdit,lineItem), autoDefault=False)
        removeBtn=gui.button(None, self, "", callback=lambda: self.removeListWidget(removeBtn,boxEdit,lineItem), autoDefault=False)
        removeBtn.setEnabled(bool(boxEdit.selectedItems()))
        boxEdit.itemSelectionChanged.connect(lambda: removeBtn.setEnabled(bool(boxEdit.selectedItems())))
        boxEdit.itemSelectionChanged.connect(lambda: self.onListWidgetSelect(boxEdit,addBtn,removeBtn,lineItem))

        addBtn.setIcon(self.addIcon)
        removeBtn.setIcon(self.removeIcon)
        
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
        
    def drawCheckBox(self,pname,box=None,layout=None, label=None, default=False):
        if not hasattr(self,pname):
            setattr(self,pname,default)
        cb=gui.checkBox(None, self,pname, label)
        layout.addWidget(cb,layout.nextRow,0,1,2)
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
