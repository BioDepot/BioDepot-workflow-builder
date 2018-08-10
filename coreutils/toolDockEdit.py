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
        QPushButton:hover:pressed { background-color: #1588c5; color: black; border-style: inset;} 
        QPushButton:disabled { background-color: lightGray; border: 1px solid gray; } 
         
        '''  
        self.setStyleSheet(css)
        self.browseIcon=QtGui.QIcon('/biodepot/Bwb_core/icons/bluefile.png')
        self.addIcon=QtGui.QIcon('/biodepot/Bwb_core/icons/add.png')
        self.removeIcon=QtGui.QIcon('/biodepot/Bwb_core/icons/remove.png')
        self.submitIcon=QtGui.QIcon('/biodepot/Bwb_core/icons/submit.png')
        self.reloadIcon=QtGui.QIcon('/biodepot/Bwb_core/icons/reload.png')
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
        elements=('User','RNA_seq','Utilities','Miscellaneous')
        cbox=self.makeComboBox('Choose category',elements)
        ledit=self.makeLedit('Enter widget name','Choose widget')
        addWidgetBtn = gui.button(None, self, "Add", callback=self.add)
        addWidgetBtn.setFixedSize(30,20)

        self.controlArea.layout().addLayout(ledit)
        self.controlArea.layout().addLayout(cbox)
        self.controlArea.layout().addWidget(addWidgetBtn)
        
    def add(self):
        pass
  #workhorse widget
    def makeLedit(self,text=None,label=None):
        leditLabel=None
        if(label):
            leditLabel=QtGui.QLabel(label)
        ledit = QtGui.QLineEdit(self)
        ledit.setClearButtonEnabled(True)
        ledit.setPlaceholderText(text)
        ledit.setStyleSheet(":disabled { color: #282828}")           
        box=QHBoxLayout()
        if leditLabel:
            box.addWidget(leditLabel)
            box.addWidget(ledit)
        setattr(box,'getValue',lambda : self.getLeditValue(ledit))
        return box
        

    def makeComboBox (self,label, elements):
        comboBoxLabel=QtGui.QLabel(label)
        comboBox=QtGui.QComboBox()
        comboBox.addItems(elements)
        comboBox.currentIndex=0
        box=QHBoxLayout()
        box.addWidget(comboBoxLabel)
        box.addWidget(comboBox)
        setattr(box,'getValue',lambda : self.getComboValue(comboBox))
        return box
        
    def getComboValue(self,comboBox):
        if comboBox.isEnabled():
            return comboBox.currentText()  
        return None 
        
    
