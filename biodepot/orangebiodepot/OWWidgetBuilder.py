import os
import re
import sys
import json
import jsonpickle
from orangebiodepot.util.createWidget import createWidget, register
from copy import deepcopy
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
    name = "WidgetBuilder"
    description = "Build a new widget from a set of bash commands and a Docker container"
    category = "Data"
    icon = "icons/build.png"
    priority = 2

    inputs = []
    outputs = []
    pset=partial(settings.Setting,schema_only=True)
    
    want_main_area = False
    want_control_area = True
    allStates=pset({})
    allAttrs=pset({})
    
    def drawExec(self, layout=None):
        css = '''
        QPushButton {background-color: #1588c5; color: white; height: 20px; border: 1px solid #1a8ac6; border-radius: 2px;}
        QPushButton:pressed { background-color: #158805; border-style: inset;}
        QPushButton:disabled { background-color: lightGray; border: 1px solid gray; }
        QPushButton:hover {background-color: #1588f5; }
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
        saveBtn = gui.button(None, self, "Save json", callback=self.saveJson)
        saveBtn.setStyleSheet(css)
        saveBtn.setFixedSize(80,20)
        createBtn = gui.button(None, self, "Create Widget", callback=self.saveWidget)
        createBtn.setStyleSheet(css)
        createBtn.setFixedSize(120,20)
        registerBtn = gui.button(None, self, "Register Widget", callback=self.registerWidget)
        registerBtn.setStyleSheet(css)
        registerBtn.setFixedSize(150,20)
        box=QtGui.QHBoxLayout()
        box.addWidget(saveBtn)
        box.addWidget(createBtn)
        box.addWidget(registerBtn)
        box.addStretch(1)
        layout.addLayout(box)
    
    def buildData(self):
        self.data={}
        for attr in ('name','description','category','docker_image_name','docker_image_tag',
            'priority','icon','inputs','outputs','volumes','parameters','command','autoMap'):
            self.data[attr]=deepcopy(getattr(self,attr))
        #separate multiple commands
        self.data['command']=self.data['command'].split('\n')
        
        #add persistance
        self.data['persistentSettings']='all'
        #add required elements
        reqList=[]
        if 'parameters' in self.data:
            for pname,pvalue in self.data['parameters'].items():
                if not pvalue['optional']:
                    reqList.append(pname)
        self.data['requiredParameters']=reqList
        #add volume mappings
        if 'volumes' in self.data and self.data['volumes']:
            myMappings=[]
            for attr, volume in self.data['volumes'].items():
                myMappings.append({'conVolume':volume['containerVolume'], 'attr' : attr})
            self.data['volumeMappings']=myMappings
            self.data.pop('volumes',None)
        #replace text str with type(str)
        for pname in ('inputs','outputs'):
            if pname in self.data:
                for key, myDict in self.data[pname].items():
                    if myDict['type'] == 'str':
                        myDict['type'] = type('str')

        #for now just keep default, flags, 
        #also make all defaults false for booleans - will fix these kluges after cleaning up BwBase code
        if 'parameters' in self.data:
            for key, myDict in self.data['parameters'].items():
                newDict={}
                if  'default' in myDict and myDict['default'] is not None:
                    #make sure that these are the correct type
                    if myDict['type'] == 'bool':
                        if myDict['default'] == 'False':
                            newDict['default']=False
                        elif myDict['default'] == 'True':
                            newDict['default']=True
                        else:
                            raise Exception ('{} is boolean - default values must be True or False not {}'.format(pname,myDict['default']))
                    elif 'type' in myDict and myDict['type'] == 'int':
                        newDict['default']=int(myDict['default'])
                    elif 'type' in myDict and (myDict['type'] == 'double' or myDict['type'] == 'float') :
                        newDict['default']=float(myDict['default'])                
                    elif  'type' in myDict and myDict['type'] == 'str' :
                        newDict['default']=str(myDict['default'])
                    else:
                        newDict['default']=json.loads(myDict['default'])
                
                if 'flag' in myDict:
                    newDict['flag']=myDict['flag']
                #arguments are the same as having a null flag value
                #parameters that are not arguments e.g. environment variables will not have a flag key    
                if 'argument' in myDict and myDict['argument'] == True:
                    newDict['flag']=None
                if 'label' in myDict:
                    newDict['label']=myDict['label']
                else:
                    newDict['label']=None
                if 'type' in myDict and myDict['type']:
                    newDict['type']=myDict['type']
                if 'env' in myDict and myDict['env']:
                    newDict['env']=myDict['env']
                self.data['parameters'][key]=newDict 
    
              
    def saveJson(self):
        self.buildData()
        #get json filename
        dataJ=jsonpickle.encode(self.data)
        jSaveFile=QtWidgets.QFileDialog.getSaveFileName(self,"QFileDialog.getSaveFileName()","","All Files (*);;json Files (*.json)")[0]
        if jSaveFile:
            sys.stderr.write('save file is {}\n'.format(jSaveFile))
            with open(jSaveFile,"w") as f:
                f.write(dataJ)
            f.close()
            
    def saveWidget(self):
        self.buildData()
        #get json filename
        outputWidget=QtWidgets.QFileDialog.getSaveFileName(self,"QFileDialog.getSaveFileName()","","Python Files (*.py);;All Files (*)")[0]
        if outputWidget:
            createWidget(None,outputWidget,inputData=self.data)
            
    def registerWidget(self):
        self.buildData()
        jsonFile=QtWidgets.QFileDialog.getOpenFileName(self,"QFileDialog.getOpenFileName()","Choose json file","json Files (*.json);;All Files (*)")[0]
        if jsonFile:
            widgetFile=QtWidgets.QFileDialog.getOpenFileName(self,"QFileDialog.getOpenFileName()","Choose widget file","Python Files (*.py);;All Files (*)")[0]
            if widgetFile:
                register(jsonFile,widgetFile,self.data['icon'])

                    
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
        
        #initialize the attr list
        for attr in ('name','description','category','docker_image_name','docker_image_tag',
                     'priority','icon','inputs','outputs','volumes','parameters','command','autoMap'):
            if attr in self.allAttrs:
                setattr(self,attr,self.allAttrs[attr])
            else:
                setattr(self,attr,None)
                self.allAttrs[attr]=None
    
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
        requiredBox = gui.widgetBox(self.bigBox, "Widget entries")
        #draw Ledits for the frequired elements
        leditRequiredLayout=QtGui.QGridLayout()
        leditRequiredLayout.setSpacing(5)
        setattr(leditRequiredLayout,'nextRow',1)        
        requiredBox.layout().addLayout(leditRequiredLayout)
        for pname in ['name','description','category','docker_image_name','docker_image_tag']:
            self.drawLedit(pname,layout=leditRequiredLayout)
        self.drawLedit('priority',layout=leditRequiredLayout)
        #file entry for icon
        self.drawLedit('icon',layout=leditRequiredLayout,addBrowseButton=True, fileType=None)

        #listwidgets for inputs and outputs 
        #define widgetItems for the different widgetLists
        #top level widgets are drawXXX - these can have multiple substituents
        #lower level widgets are makeXXX - these can also have multiple substituents
        
        self.drawIListWidget('inputs',layout=leditRequiredLayout)
        self.drawOListWidget('outputs',layout=leditRequiredLayout)
        self.drawVolumeListWidget('volumes',layout=leditRequiredLayout)
        self.drawParamsListWidget('parameters',layout=leditRequiredLayout)
        self.drawCommand('command',layout=leditRequiredLayout)
        self.drawExec(self.controlArea.layout())
    
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
            setattr(self,pname,OrderedDict(widgetList))
        else:
            self.allStates[pname]=None
            setattr(self,pname,qWidgetItem.getStateValue(None))
            
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
            setattr(self,pname,OrderedDict(widgetList))
        else:
            setattr(self,pname,None)
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
            
             
    def drawCommand(self,pname,layout=None):
        labelTextBox=self.makeTextBox(pname,label='Enter command:')
        self.initAllStates(pname,labelTextBox)
        labelTextBox.textBox.textChanged.connect(lambda: self.updateAllStates(pname,labelTextBox,labelTextBox.getState()))
        layout.addWidget(labelTextBox.label,layout.nextRow,0)
        layout.addWidget(labelTextBox.textBox,layout.nextRow,1,1,4)
        layout.nextRow = layout.nextRow + 1
    
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
        if not hasattr(self,attr):
            setattr(self,attr,default)
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
        if not hasattr(self,pname):
            setattr(self,pname,None); 
        boxEdit=QtGui.QListWidget(self)
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
        #buttons
        addBtn=gui.button(None, self, "", callback=lambda: self.addListWidget(addBtn,boxEdit,lineItem), autoDefault=False)
        removeBtn=gui.button(None, self, "", callback=lambda: self.removeListWidget(removeBtn,boxEdit,lineItem), autoDefault=False)
        removeBtn.setEnabled(bool(boxEdit.selectedItems()))
        boxEdit.itemSelectionChanged.connect(lambda: removeBtn.setEnabled(bool(boxEdit.selectedItems())))
        boxEdit.itemSelectionChanged.connect(lambda: self.onListWidgetSelect(boxEdit,addBtn,removeBtn,lineItem))
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
        label=QtGui.QLabel(pname+':')
        label.setAlignment(Qt.AlignTop)
        layout.addWidget(label,layout.nextRow,0)
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

    def drawLedit(self,pname,layout=None,addBrowseButton=False,fileType=None):
        #make labeledLedit combo layout
        labelLedit=self.makeLedit(pname,'Enter '+ pname, label=pname+':')
        self.initAllStates(pname,labelLedit)
        labelLedit.ledit.textChanged.connect(lambda: self.updateAllStates(pname,labelLedit,labelLedit.getState()))
        layout.addWidget(labelLedit.label,layout.nextRow,0)
        if addBrowseButton:
            layout.addWidget(labelLedit.ledit,layout.nextRow,1,1,1)
            button=gui.button(None, self, "", callback= lambda : self.browseFileDir(pname,ledit=labelLedit.ledit,fileType=fileType),autoDefault=True, width=19, height=19)
            button.setIcon(self.browseIcon)
            layout.addWidget(button,layout.nextRow,2)
        else:
            layout.addWidget(labelLedit.ledit,layout.nextRow,1,1,2)
        layout.nextRow+=1
        return labelLedit
    
    def updateAllStates(self,pname,widget,state):
        self.allStates[pname]=state
        setattr(self,pname,widget.getStateValue(state))
        sys.stderr.write('updating all - pname {} state {}\n'.format(pname,state))
        
    def initAllStates(self, attr ,widget):
        if attr in self.allStates:
            widget.setState(self.allStates[attr])
        else:
            widget.setState(None)
            self.allStates[attr]=widget.getState()
        setattr(self,attr,widget.getStateValue(self.allStates[attr]))
            
    def browseFileDir(self, attr,ledit=None,fileType=None):
        defaultDir = '/root'
        if os.path.exists('/data'):
            defaultDir = '/data'
        if fileType == 'Directory':
            myFileDir=QtWidgets.QFileDialog.getExistingDirectory(self, caption="Locate directory", directory=defaultDir)
        else:
            myFileDir=QtWidgets.QFileDialog.getOpenFileName(self, "Locate file", defaultDir)[0]
        if myFileDir:
            setattr(self,attr,myFileDir)
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
