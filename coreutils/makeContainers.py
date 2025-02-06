import os
import re
import sys
import json
import jsonpickle
import pickle
import csv,time
import tempfile, shutil
from xml.dom import minidom
from glob import glob
from pathlib import Path
from shutil import copyfile
from createWidget import mergeWidget, createWidget, findIconFile
from copy import deepcopy
from collections import OrderedDict
from functools import partial
from AnyQt.QtCore import QThread, pyqtSignal, Qt
from Orange.widgets import widget, gui, settings
from DockerClient import DockerClient, PullImageThread, ConsoleProcess
from PyQt5 import QtWidgets, QtGui, QtCore
from PyQt5.QtCore import QProcess
from PyQt5.QtWidgets import *
from PyQt5.QtGui import *

class PullBuildProcess:
    def __init__(self,console,workflow,forceLoad,preferDockerfile,maxAttempts,finishHandler=None):
        self.process = QProcess()
        self.console = console
        self.workflow=workflow
        self.forceLoad=forceLoad
        self.preferDockerfile=preferDockerfile
        self.maxAttempts=maxAttempts
        self.state = "stopped"
        self.finishHandler=finishHandler
        self.userStopped=False
        if console:
            self.process.readyReadStandardOutput.connect(
                lambda: self.writeConsole(
                    self.process, console, self.process.readAllStandardOutput, Qt.white
                )
            )
            self.process.readyReadStandardError.connect(
                lambda: self.writeConsole(
                    self.process, console, self.process.readAllStandardError, Qt.red
                )
            )        
        self.process.finished.connect(self.onFinish)
        for var in ["workflow","forceLoad","preferDockerfile","maxAttempts"]:
            self.writeMessage("{} {}".format(var,getattr(self,var)))

    def writeConsole(self, process, console, read, color):
        console.setTextColor(color)
        console.append(read().data().decode("utf-8", errors="ignore").rstrip())

    def writeMessage(self, message, color=Qt.green):
        self.console.setTextColor(color)
        self.console.append(message)

    def stop(self, message=None):
        self.state = "stopped"
        #the build process should clean itself up
        self.userStopped=True
        self.process.terminate()
        if message:
            self.writeMessage(message,color=Qt.red)
        if self.finishHandler:
            self.finishHandler()
            
    def start(self):
        self.writeMessage('Starting build of {}\n'.format(self.workflow))
        cmdStr='build_workflow_containers.sh -w {} -m {} '.format(self.workflow,self.maxAttempts)
        if self.forceLoad:
            cmdStr+='-f '
        if self.preferDockerfile:
            cmdStr+='-d '
        self.writeMessage(cmdStr)
        self.process.start(cmdStr)
    
        
    def onFinish(self,code,status):
        if not code and not self.userStopped:   
            self.writeMessage('Successfully finished build\n')
        if self.finishHandler:
            self.finishHandler()
        

def categoryWidgetPath(categoryDir):
    pyFiles = glob("/biodepot/{}/OW*.py".format(categoryDir))
    if not pyFiles:
        return None
    widgetPath = os.path.dirname((os.readlink(pyFiles[0])))
    # check if absolute path
    if widgetPath[0] == "/":
        absWidgetPath = widgetPath
    else:
        absWidgetPath = os.path.abspath(
            "/biodepot/{}/{}".format(categoryDir, widgetPath)
        )
    return os.path.dirname(absWidgetPath)

class BuildContainers(widget.OWWidget):
    name = "ContainerBuilder"
    description = "Pulls and builds containers"
    category = "Utility"
    icon = "icons/build.png"
    priority = 2
    inputs = []
    outputs = []
    pset = partial(settings.Setting, schema_only=True)
    want_main_area = False
    want_control_area = True
    allStates = pset({})
    allAttrs = pset({})
    data = {}

    def __init__(self, renameData=None, canvasMainWindow=None):
        super().__init__()
        self.canvas = canvasMainWindow
        if renameData:
            # run the rename routine and then finish
            self.updateCategories()
            self.widgetRename(
                widgetName=renameData["widgetName"],
                category=renameData["category"],
                newName=renameData["newName"],
            )
            return
        self.css = """
        QPushButton {background-color: #1588c5; color: white; height: 20px; border: 1px solid black; border-radius: 2px;}
        QPushButton:hover {background-color: #1555f5; }
        QPushButton:hover:pressed { background-color: #1588c5; color: black; border-style: inset; border: 1px solid white} 
        QPushButton:disabled { background-color: lightGray; border: 1px solid gray; }
        """
        self.baseToolPath = "/biodepot"
        self.setStyleSheet(self.css)
        self.browseIcon = QIcon("/icons/bluefile.png")
        self.addIcon = QIcon("/icons/add.png")
        self.removeIcon = QIcon("/icons/remove.png")
        self.submitIcon = QIcon("/icons/submit.png")
        self.reloadIcon = QIcon("/icons/reload.png")
        self.controlArea.setMinimumWidth(500)
        self.controlArea.setMinimumHeight(120)
        self.startWidget()

    def clearLayout(self, layout):
        if layout != None:
            while layout.count():
                child = layout.takeAt(0)
                if child.widget() is not None:
                    child.widget().deleteLater()
                elif child.layout() is not None:
                    self.clearLayout(child.layout())

    def startWidget(self):
        self.setWindowTitle("Load Containers")
        self.grid = QGridLayout()
        self.execute_box=QHBoxLayout()
        self.controlArea.layout().addLayout(self.grid)
        self.initCategories()
        self.drawLoadContainers()
        
    def updateCategories(self):
        temp_categories = (
            str(os.popen("""grep -oP 'name="\K[^"]+' /biodepot/setup.py""").read())
        ).split()
        # directories are not same as categories because Python/Linux unfriendly characters are changed
        temp_directoryList = (
            str(
                os.popen("""grep -oP 'packages=\["\K[^"]+' /biodepot/setup.py""").read()
            )
        ).split()
        self.categoryToDirectory = {}
        self.categoryToPath = {}
        self.categories = []
        # check for icon link for categories - it is a symLink for workflows and not built-ins
        for index, category in enumerate(temp_categories):
            iconLink = "/biodepot/{}/icon".format(temp_directoryList[index])
            if os.path.islink(iconLink):
                self.categoryToDirectory[category] = temp_directoryList[index]
                self.categories.append(category)
                iconPath = os.readlink(iconLink)
                self.categoryToPath[category] = os.path.dirname(
                    os.path.dirname(os.path.dirname(iconPath))
                )

    def initCategories(self):
        self.updateCategories()
        self.clearLayout(self.controlArea.layout())
      
    def drawLoadContainers(self):
        self.controlArea.layout().addLayout(self.grid)
        self.RWcbox = self.makeComboBox(
            self.grid,
            "Choose workflow:",
            self.categories,
            startRow=1,
            callback=self.updateRWcbox
        )
        self.workflow_path=self.categoryToPath[self.categories[0]]
        self.forceLoad=False
        forceLoadcb=gui.checkBox(None, self, "forceLoad", label="Reload even if container is present")
        self.grid.addWidget(forceLoadcb, 3, 1)
        self.preferFile=False
        preferFilecb=gui.checkBox(None, self, "preferFile", label="Try building from Dockerfile first")
        self.grid.addWidget(preferFilecb, 4, 1)
        self.maxAttempts=1
        maxAttemptsLabel = QLabel("Max attempts at pulling image: ")
        self.maxAttemptsSpin = QtWidgets.QSpinBox(self)
        self.maxAttemptsSpin.setSuffix(' attempts')
        self.maxAttemptsSpin.setRange(1,10)
        self.grid.addWidget(maxAttemptsLabel,2,1)    
        self.grid.addWidget(self.maxAttemptsSpin, 2,2)
        self.loadBtn = gui.button(None, self, "Load", callback=self.onLoadClicked)
        self.loadBtn.setFixedSize(40, 20)
        self.loadBtn.setStyleSheet(self.css)
        self.stopBtn = gui.button(None, self, "Stop", callback=self.onStopClicked)
        self.stopBtn.setFixedSize(40, 20)
        self.stopBtn.setStyleSheet(self.css)
        self.stopBtn.setEnabled(False)
        self.console = QTextEdit()
        self.console.setReadOnly(True)
        pal = QPalette()
        pal.setColor(QPalette.Base, Qt.black)
        pal.setColor(QPalette.Text, Qt.green)
        self.console.setPalette(pal)
        self.console.setAutoFillBackground(True)
        self.controlArea.layout().addWidget(self.console)
        self.controlArea.layout().addLayout(self.execute_box)
        self.execute_box.addWidget(self.loadBtn)
        self.execute_box.addWidget(self.stopBtn)
        self.execute_box.addStretch(2)
        
    def updateRWcbox(self):
        sys.stderr.write("current category is {}\n".format(self.RWcbox.currentText()))
        sys.stderr.write("current path is {}\n".format(self.categoryToPath[self.RWcbox.currentText()]))
        self.workflow_path=self.categoryToPath[self.RWcbox.currentText()]
        
    def onLoadClicked(self):
        self.loadBtn.setEnabled(False)
        self.stopBtn.setEnabled(True)
        self.start()

    def onStopClicked(self):
        self.loadBtn.setEnabled(True)
        self.stopBtn.setEnabled(False)
        self.pullBuildProc.stop("User stopped build")
    
    def onBuildFinish(self):
        self.loadBtn.setEnabled(True)
        self.stopBtn.setEnabled(False)
        
    def start(self):
        self.pullBuildProc = PullBuildProcess(self.console,self.workflow_path,self.forceLoad,self.preferFile,self.maxAttemptsSpin.value(),finishHandler=self.onBuildFinish)
        self.pullBuildProc.start()

    def makeLedit(
        self,
        layout,
        text=None,
        label=None,
        startRow=1,
        startColumn=1,
        browse=False,
        browseFileFlag=False,
    ):
        leditLabel = None
        if label:
            leditLabel = QLabel(label)
        ledit = QLineEdit(self)
        ledit.setClearButtonEnabled(True)
        ledit.setPlaceholderText(text)
        ledit.setStyleSheet(":disabled { color: #282828}")
        layout.addWidget(leditLabel, startRow, startColumn)
        layout.addWidget(ledit, startRow, startColumn + 1)
        if browse:
            button = gui.button(
                None,
                self,
                "",
                callback=lambda: self.browseWidget(ledit, browseFileFlag),
                autoDefault=True,
                width=19,
                height=19,
            )
            button.setIcon(self.browseIcon)
            layout.addWidget(button, startRow, startColumn + 2)
        return ledit

    def browseWidget(self, ledit, browseFileFlag=False):
        if browseFileFlag:
            ledit.setText(
                str(
                    QtWidgets.QFileDialog.getOpenFileName(
                        self, "Locate file", "/widgets"
                    )[0]
                )
            )
            return
        myFileDir = QtWidgets.QFileDialog.getExistingDirectory(
            self, caption="Locate widget", directory="/widgets"
        )
        ledit.setText(str(myFileDir))

    def makeComboBox(
        self, layout, label, elements, startRow=1, startColumn=1, callback=None
    ):
        comboBoxLabel = QLabel(label)
        comboBox = QComboBox()
        if elements:
            comboBox.addItems(elements)
        comboBox.currentIndex = 0
        if callback:
            comboBox.currentIndexChanged.connect(callback)
        layout.addWidget(comboBoxLabel, startRow, startColumn)
        layout.addWidget(comboBox, startRow, startColumn + 1)
        return comboBox


    def getComboValue(self, comboBox):
        if comboBox.isEnabled():
            return comboBox.currentText()
        return None

    def getLeditValue(self, ledit):
        if ledit.isEnabled():
            return ledit.text()
        return None
