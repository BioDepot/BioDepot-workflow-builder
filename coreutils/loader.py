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
class ToolDockEdit(widget.OWWidget):
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
    
    def __init__(self,widgetID=None):
        super().__init__()
        css = '''
        QPushButton {background-color: #1588c5; color: white; height: 20px; border: 1px solid black; border-radius: 2px;}
        QPushButton:hover {background-color: #1555f5; }
        QPushButton:hover:pressed { background-color: #1588c5; color: black; border-style: inset; 1px solid white} 
        QPushButton:disabled { background-color: lightGray; border: 1px solid gray; } 
         
        '''  
        self.setStyleSheet(css)
        self.browseIcon=QtGui.QIcon('/biodepot/orangebiodepot/icons/bluefile.png')
        self.addIcon=QtGui.QIcon('/biodepot/orangebiodepot/icons/add.png')
        self.removeIcon=QtGui.QIcon('/biodepot/orangebiodepot/icons/remove.png')
        self.submitIcon=QtGui.QIcon('/biodepot/orangebiodepot/icons/submit.png')
        self.reloadIcon=QtGui.QIcon('/biodepot/orangebiodepot/icons/reload.png')
        self.startWidget()

    
    def clearLayout(self,layout):
        while layout.count():
            child = layout.takeAt(0)
            if child.widget():
                child.widget().deleteLater()
    
    def startWidget(self):
        self.setWindowTitle('Add Widget to Dock')
        #self.setStyleSheet(":disabled { color: #282828}")
        
    def addWidget(self):
        categories=(str(os.popen('''grep -oP 'name="\K[^"]+' /biodepot/setup.py''').read())).split()
        grid=QtGui.QGridLayout()
        self.categoryValue=self.makeComboBox(grid,'Choose category:',categories)
        self.widgetDir=self.makeLedit(grid,'Enter widget name','Choose widget:')        
        self.controlArea.layout().addLayout(grid)

    def widgetAdd(self):
        inputDir=self.widgetDir()
        category=self.categoryValue()
        sys.stderr.write('input widget is {} category is {}\n'.format(inputDir,category))
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
            os.system ("ln -sf  /widgets/{}/{}.py /biodepot/{}/OW{}.py".format(inputName,inputName,category,inputName))
            title='Add {}'.format(inputName)
            message='Added {} to {} in ToolDock'.format(inputName,category)
            QtGui.QMessageBox.information(self, title,message,QtGui.QMessageBox.Ok)
        except:
            pass
    def removeWidget(self):
        pass
    def makeLedit(self,layout,text=None,label=None):
        leditLabel=None
        if(label):
            leditLabel=QtGui.QLabel(label)
        ledit = QtGui.QLineEdit(self)
        ledit.setClearButtonEnabled(True)
        ledit.setPlaceholderText(text)
        ledit.setStyleSheet(":disabled { color: #282828}")           
        button=gui.button(None, self, "",callback= lambda: self.browseWidget(ledit),autoDefault=True, width=19, height=19)
        button.setIcon(self.browseIcon)
        widgetAddBtn = gui.button(None, self, "Add", callback=self.widgetAdd)
        widgetAddBtn.setFixedSize(30,20)
        layout.addWidget(leditLabel,1,1)
        layout.addWidget(ledit,1,2)
        layout.addWidget(button,1,3)
        layout.addWidget(widgetAddBtn,1,4)
        return lambda : self.getLeditValue(ledit)
        
    def browseWidget(self, ledit):
        myFileDir=QtWidgets.QFileDialog.getExistingDirectory(self, caption="Locate widget", directory='/widgets')
        ledit.setText(myFileDir)
            
    def makeComboBox (self,layout,label, elements):
        comboBoxLabel=QtGui.QLabel(label)
        comboBox=QtGui.QComboBox()
        comboBox.addItems(elements)
        comboBox.currentIndex=0
        layout.addWidget(comboBoxLabel,2,1)
        layout.addWidget(comboBox,2,2)
        return lambda : self.getComboValue(comboBox)
        
    def getComboValue(self,comboBox):
        if comboBox.isEnabled():
            return comboBox.currentText()  
        return None 
        
    def getLeditValue(self,ledit):
        if ledit.isEnabled():
            return ledit.text()
        return None
