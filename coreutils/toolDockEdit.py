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
from PyQt5 import QtWidgets, QtGui, QtCore
from PyQt5.QtWidgets import QInputDialog, QLineEdit, QMainWindow, QApplication, QPushButton, QWidget, QAction, QTabWidget, QVBoxLayout, QMessageBox

from AnyQt.QtWidgets import (
    QWidget, QButtonGroup, QGroupBox, QRadioButton, QSlider,
    QDoubleSpinBox, QComboBox, QSpinBox, QListView, QLabel,
    QScrollArea, QVBoxLayout, QHBoxLayout, QFormLayout,
    QSizePolicy, QApplication, QCheckBox
)
defaultIconFile='/biodepot/Bwb_core/icons/default.png'     
class ToolDockEdit(widget.OWWidget):
    name = "ToolDockEditor"
    description = "Edits contents of tool dock"
    category = "Utility"
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
    
    def __init__(self,widgetID=None):
        super().__init__()
        self.css = '''
        QPushButton {background-color: #1588c5; color: white; height: 20px; border: 1px solid black; border-radius: 2px;}
        QPushButton:hover {background-color: #1555f5; }
        QPushButton:hover:pressed { background-color: #1588c5; color: black; border-style: inset; border: 1px solid white} 
        QPushButton:disabled { background-color: lightGray; border: 1px solid gray; }
        '''
        self.setStyleSheet(self.css)
        self.browseIcon=QtGui.QIcon('/biodepot/orangebiodepot/icons/bluefile.png')
        self.addIcon=QtGui.QIcon('/biodepot/orangebiodepot/icons/add.png')
        self.removeIcon=QtGui.QIcon('/biodepot/orangebiodepot/icons/remove.png')
        self.submitIcon=QtGui.QIcon('/biodepot/orangebiodepot/icons/submit.png')
        self.reloadIcon=QtGui.QIcon('/biodepot/orangebiodepot/icons/reload.png')
        self.controlArea.setMinimumWidth(500)
        self.controlArea.setMinimumHeight(60)
        self.startWidget()

    
    def clearLayout(self,layout):
        if layout != None:
            while layout.count():
                child = layout.takeAt(0)
                if child.widget() is not None:
                    child.widget().deleteLater()
                elif child.layout() is not None:
                    self.clearLayout(child.layout())
    
    def startWidget(self):
        self.setWindowTitle('WidgetToolDock Editor')
        #self.setStyleSheet(":disabled { color: #282828}")
        
    def addWidget(self):
        self.categories=(str(os.popen('''grep -oP 'name="\K[^"]+' /biodepot/setup.py''').read())).split()
        #directories are not same as categories because Python/Linux unfriendly characters are changed
        directoryList=(str(os.popen('''grep -oP 'packages=\["\K[^"]+' /biodepot/setup.py''').read())).split()
        self.categoryToDirectory={}
        for index, category in enumerate(self.categories):
            self.categoryToDirectory[category]=directoryList[index]
        self.clearLayout(self.controlArea.layout())
        if not hasattr(self,'grid'):
            self.grid=QtGui.QGridLayout()
        self.clearLayout(self.grid)
        ledit=self.makeLedit(self.grid,'Enter widget name','Choose widget ',startRow=1,startColumn=1,browse=True)
        cbox=self.makeComboBox(self.grid,'from category:',self.categories,startRow=1,startColumn=4)
        widgetAddBtn = gui.button(None, self, "Add", callback= lambda: self.widgetAdd(ledit,cbox))
        widgetAddBtn.setFixedSize(30,20)
        widgetAddBtn.setStyleSheet(self.css)
        self.grid.addWidget(widgetAddBtn,1,6)        
        self.controlArea.layout().addLayout(self.grid)

    def widgetAdd(self,ledit,cbox):
        category=self.getComboValue(cbox)
        inputDir=self.getLeditValue(ledit)
        if not inputDir:
            return
        inputName=os.path.basename(os.path.normpath(inputDir))
        if not os.path.exists(inputDir):
            if os.path.exists('/widgets/{}'.format(inputName)):
            #ask if we want to overwrite
                qm = QtGui.QMessageBox
                ret=qm.question(self,'', "{} exists - OverWrite ?".format(inputName), qm.Yes | qm.No)
                if ret == qm.No:
                    return
            #safer way of removing the widget
                os.system("cd /widgets && rm -rf {}".format(inputName))    
            os.system('cp -r {} /widgets/{}'.format(inputDir,inputName))
        try:
            directory=self.categoryToDirectory[category]
            os.system ("ln -sf  /widgets/{}/{}.py /biodepot/{}/OW{}.py".format(inputName,inputName,directory,inputName))
            title='Add {}'.format(inputName)
            message='Added {} to {} in ToolDock'.format(inputName,category)
            QtGui.QMessageBox.information(self, title,message,QtGui.QMessageBox.Ok)
        except:
            pass
            
    def removeWidget(self):
        self.categories=(str(os.popen('''grep -oP 'name="\K[^"]+' /biodepot/setup.py''').read())).split()
        #directories are not same as categories because Python/Linux unfriendly characters are changed
        directoryList=(str(os.popen('''grep -oP 'packages=\["\K[^"]+' /biodepot/setup.py''').read())).split()
        self.categoryToDirectory={}
        for index, category in enumerate(self.categories):
            self.categoryToDirectory[category]=directoryList[index]
        widgetList=(str(os.popen('cd /widgets && ls -d *').read())).split()
        if not hasattr(self,'grid'):
            self.grid=QtGui.QGridLayout()
        self.clearLayout(self.controlArea.layout())
        self.clearLayout(self.grid)
        self.wbox=self.makeComboBox(self.grid,'Remove widget ',widgetList,startRow=1,callback=self.__onWidgetChange)
        categoryList=self.getCategoryList(self.getComboValue(self.wbox))
        self.cbox=self.makeComboBox(self.grid,' from category ',categoryList,startRow=1,startColumn=4)
        widgetRemoveBtn = gui.button(None, self, "Remove", callback=self.widgetRemove)
        widgetRemoveBtn.setFixedSize(60,20)
        widgetRemoveBtn.setStyleSheet(self.css)
        self.grid.addWidget(widgetRemoveBtn,1,6)
        self.controlArea.layout().addLayout(self.grid)
        
    def widgetRemove(self):
        qm = QtGui.QMessageBox
        widgetName=self.getComboValue(self.wbox)
        if not os.path.exists('/widgets/{}'.format(widgetName)):
            return
        category=self.getComboValue(self.cbox)
        ret=qm.No
        if category == '_ALL_':
            ret=qm.question(self,'', "Remove {} completely ?".format(widgetName), qm.Yes | qm.No)
            if ret == qm.Yes:
                os.system('cd /widgets && rm {} -rf'.format(widgetName))
                #cleanup links
                os.system('find -L /biodepot -type l -delete')
                qm.information(self,'Removed widget','Completely removed widget {}'.format(widgetName),QtGui.QMessageBox.Ok)
            #remove from widgetList and set text to new top
        else:
            ret=qm.question(self,'', "Remove {} from category {} ?".format(widgetName,category), qm.Yes | qm.No)
            if ret == qm.Yes:
                directory=self.categoryToDirectory[category]
                os.system('cd /biodepot/{} && rm OW{}.py '.format(directory,widgetName))
                qm.information(self,'Removed widget','Removed widget {} from {}'.format(widgetName,category),QtGui.QMessageBox.Ok)
    
    def addCategory(self):
        if not hasattr(self,'grid'):
            self.grid=QtGui.QGridLayout()
        self.clearLayout(self.controlArea.layout())
        self.clearLayout(self.grid)
        nameLedit=self.makeLedit(self.grid,'Enter category name','Choose category:',startRow=1)
        iconLedit=self.makeLedit(self.grid,'Enter icon file','Choose icon:',startRow=2,browse=True)
        widgetAddBtn = gui.button(None, self, "Add", callback= lambda: self.categoryAdd(nameLedit,iconLedit))
        widgetAddBtn.setFixedSize(30,20)
        widgetAddBtn.setStyleSheet(self.css)
        self.grid.addWidget(widgetAddBtn,2,4)        
        self.controlArea.layout().addLayout(self.grid) 
        
    def categoryAdd(self,nameLedit,iconLedit):
        pass
    def getCategoryList(self,widgetName):
        #categories may not be directories
        if not widgetName:
            return self.categories
        returnList=[]
        for category in self.categories:
            directory=self.categoryToDirectory[category]
            sys.stderr.write('/biodepot/{}/OW{}.py\n'.format(directory,widgetName))
            if os.path.exists('/biodepot/{}/OW{}.py'.format(directory,widgetName)):
                returnList.append(category)
        returnList.insert(0,'_ALL_')
        return returnList

    def makeLedit(self,layout,text=None,label=None,startRow=1,startColumn=1,browse=False):
        leditLabel=None
        if(label):
            leditLabel=QtGui.QLabel(label)
        ledit = QtGui.QLineEdit(self)
        ledit.setClearButtonEnabled(True)
        ledit.setPlaceholderText(text)
        ledit.setStyleSheet(":disabled { color: #282828}")           
        layout.addWidget(leditLabel,startRow,startColumn)
        layout.addWidget(ledit,startRow,startColumn+1)
        if browse:
            button=gui.button(None, self, "",callback= lambda: self.browseWidget(ledit),autoDefault=True, width=19, height=19)
            button.setIcon(self.browseIcon)
            layout.addWidget(button,startRow,startColumn+2)
        return ledit
        
    def browseWidget(self, ledit):
        myFileDir=QtWidgets.QFileDialog.getExistingDirectory(self, caption="Locate widget", directory='/widgets')
        ledit.setText(myFileDir)
            
    def makeComboBox (self,layout,label, elements,startRow=1,startColumn=1,callback=None):
        comboBoxLabel=QtGui.QLabel(label)
        comboBox=QtGui.QComboBox()
        if elements:
            comboBox.addItems(elements)
        comboBox.currentIndex=0
        if callback:
            comboBox.currentIndexChanged.connect(lambda: callback(comboBox)) 
        layout.addWidget(comboBoxLabel,startRow,startColumn)
        layout.addWidget(comboBox,startRow,startColumn+1)
        return comboBox
    
    def __onWidgetChange(self,box):
        self.cbox.clear()
        self.cbox.addItems(self.getCategoryList(box.currentText()))
        
        
    def getComboValue(self,comboBox):
        if comboBox.isEnabled():
            return comboBox.currentText()  
        return None 
        
    def getLeditValue(self,ledit):
        if ledit.isEnabled():
            return ledit.text()
        return None
