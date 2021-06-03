import os
import re
import sys
import json
import jsonpickle
import pickle
import csv
import tempfile, shutil
import OWImageBuilder, toolDockEdit
from glob import glob
from pathlib import Path
from shutil import copyfile, copytree
from createWidget import mergeWidget, createWidget, findIconFile
from copy import deepcopy
from collections import OrderedDict
from functools import partial
from AnyQt.QtCore import QThread, pyqtSignal, Qt
from Orange.widgets import widget, gui, settings
from DockerClient import DockerClient, PullImageThread, ConsoleProcess
from PyQt5 import QtWidgets, QtGui, QtCore
from PyQt5.QtWidgets import *
from makeToolDockCategories import niceForm

# TODO
# fix the intialization step by checking what fields are necessary for each tab
# and checking if those fields are defined in the pickle
# if not then fill the empties with default values or None
from AnyQt.QtWidgets import (
    QWidget,
    QButtonGroup,
    QGroupBox,
    QRadioButton,
    QSlider,
    QDoubleSpinBox,
    QComboBox,
    QSpinBox,
    QListView,
    QLabel,
    QScrollArea,
    QVBoxLayout,
    QHBoxLayout,
    QFormLayout,
    QSizePolicy,
    QApplication,
    QCheckBox,
)

defaultIconFile = "/icons/default.png"
def breakpoint(title=None, message=None):
    QtGui.QMessageBox.warning(None,title,message)
    return


class SaveWorkflowForm(QDialog):
    def __init__(self, returnData, parent=None):
        self.browseCSS = """
        QPushButton {background-color: rgba(30,30,200,128); color: white; height: 20px; border: 1px solid black; border-radius: 2px;}
        QPushButton:hover {background-color: #1555f5; }
        QPushButton:hover:pressed { background-color: #1588c5; color: black; border-style: inset; border: 1px solid white} 
        QPushButton:disabled { background-color: lightGray; border: 1px solid gray; }
        """
        self.colorCSS = """
        QPushButton {background-color: white; color: white; height: 20px; border: 1px solid black; border-radius: 2px;}
        QPushButton:hover {background-color: gray; }
        QPushButton:hover:pressed { background-color: black; color: black; border-style: inset; border: 1px solid white} 
        QPushButton:disabled { background-color: lightGray; border: 1px solid gray; }
        """
        self.defaultDir = "/"
        self.returnData = returnData
        self.returnData["success"] = False
        self.initialData = returnData.copy()
        if "start_dir" in self.initialData and self.initialData["start_dir"]:
            self.defaultDir = self.initialData["start_dir"]
        self.browseIcon = QtGui.QIcon("/icons/bluefile.png")
        self.colorIcon = QtGui.QIcon("/icons/colorWheel.png")
        self.ledits = {}
        self.required = ["name", "dir"]
        super(SaveWorkflowForm, self).__init__(parent)
        self.setWindowTitle("Save workflow")
        name_ledit = QLineEdit()
        self.ledits["name"] = name_ledit
        name_ledit.setClearButtonEnabled(True)
        name_ledit.setPlaceholderText("Enter workflow name")
        if "name" in self.initialData and self.initialData["name"]:
            name_ledit.setText(self.initialData["name"])
        name_ledit.setStyleSheet(":disabled { color: #282828}")
        name_label = QtGui.QLabel("Workflow name:")
        name_box = QHBoxLayout()
        name_box.addWidget(name_label)
        name_box.addWidget(name_ledit)

        # ledit for directory
        dir_ledit = QLineEdit()
        self.ledits["dir"] = dir_ledit
        dir_ledit.setClearButtonEnabled(True)
        dir_ledit.setPlaceholderText("Enter workflow directory")
        if "dir" in self.initialData and self.initialData["dir"]:
            dir_ledit.setText(self.initialData["dir"])
        dir_ledit.setStyleSheet(":disabled { color: #282828}")
        dir_label = QtGui.QLabel("Workflow parent directory:")
        dir_button = gui.button(
            None,
            self,
            "",
            callback=lambda: self.browseFileDir(
                "dirname", ledit=dir_ledit, fileType="Directory"
            ),
            autoDefault=True,
            width=19,
            height=19,
        )
        dir_button.setIcon(self.browseIcon)
        dir_button.setStyleSheet(self.browseCSS)
        dir_box = QHBoxLayout()
        dir_box.addWidget(dir_label)
        dir_box.addWidget(dir_ledit)
        dir_box.addWidget(dir_button)

        # ledit for iconFile
        icon_ledit = QLineEdit()
        self.ledits["icon"] = icon_ledit
        icon_ledit.setClearButtonEnabled(True)
        icon_ledit.setPlaceholderText("Enter iconFile")
        if "icon" in self.initialData and self.initialData["icon"]:
            icon_ledit.setText(self.initialData["icon"])
        else:
            icon_ledit.setText("")
        icon_ledit.setStyleSheet(":disabled { color: #282828}")
        icon_label = QtGui.QLabel("Change workflow icon:")
        icon_button = gui.button(
            None,
            self,
            "",
            callback=lambda: self.browseFileDir("iconFile", ledit=icon_ledit),
            autoDefault=True,
            width=19,
            height=19,
        )
        icon_button.setIcon(self.browseIcon)
        icon_button.setStyleSheet(self.browseCSS)
        icon_box = QHBoxLayout()
        icon_box.addWidget(icon_label)
        icon_box.addWidget(icon_ledit)
        icon_box.addWidget(icon_button)

        # ledit for color
        color_ledit = QLineEdit()
        self.ledits["color"] = color_ledit
        color_ledit.setClearButtonEnabled(True)
        color_ledit.setPlaceholderText("Enter color")
        if "color" in self.initialData and self.initialData["color"]:
            color_ledit.setText(self.initialData["color"])
        color_ledit.setStyleSheet(":disabled { color: #282828}")
        color_label = QtGui.QLabel("Change workflow color:")
        color_button = gui.button(
            None,
            self,
            "",
            callback=lambda: self.browseColor("color", ledit=color_ledit),
            autoDefault=True,
            width=19,
            height=19,
        )
        color_button.setIcon(self.colorIcon)
        color_button.setStyleSheet(self.colorCSS)
        color_box = QHBoxLayout()
        color_box.addWidget(color_label)
        color_box.addWidget(color_ledit)
        color_box.addWidget(color_button)

        # checkboxes
        self.merge_checkBox = QtGui.QCheckBox("Merge all widget types", self)
        self.merge_checkBox.stateChanged.connect(
            lambda: self.getCheckBoxState(self.merge_checkBox, "merge")
        )
        if "merge" in self.initialData and self.initialData["merge"] is not None:
            self.merge_checkBox.setCheckState(self.booleanToQtChecked(self.initialData["merge"]))
        else:
            self.merge_checkBox.setCheckState(Qt.Checked)

        # OK cancel buttons
        buttonbox = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        buttonbox.accepted.connect(self.checkFields)
        buttonbox.rejected.connect(self.cancel)

        layout = QtGui.QGridLayout()
        layout.addLayout(name_box, 1, 1, 1, 2)
        layout.addLayout(dir_box, 2, 1, 1, 2)
        layout.addLayout(color_box, 3, 1, 1, 2)
        layout.addLayout(icon_box, 4, 1, 1, 2)
        layout.addWidget(self.merge_checkBox, 5, 1, 1, 1)
        layout.addWidget(buttonbox, 6, 1, 1, 1)
        self.setLayout(layout)

        # save/cancel buttons
    def booleanToQtChecked(self,state):
        if state:
            return Qt.Checked
        return Qt.Unchecked
        
    def cancel(self):
        self.returnData = self.initialData
        self.returnData["success"] = False
        self.close()

    def checkFields(self):
        # make sure that the ledits are populated with something if required
        for attr in ("name", "dir", "color", "icon"):
            ledit = self.ledits[attr]
            if not ledit.text():
                if attr in self.required:
                    warning = QMessageBox()
                    warning.setText("Must enter a value for {}".format(attr))
                    warning.setWindowTitle("Missing required entry")
                    warning.setStandardButtons(QMessageBox.Ok)
                    warning.exec_()
                    return
            else:
                self.returnData[attr] = ledit.text()
        self.returnData["merge"] = self.merge_checkBox.isChecked()
        # make sure that the name is nice Form (can keep dashes until it is converted to a ows file
        self.returnData["name"] = niceForm(self.returnData["name"], useDash=True)
        self.returnData["success"] = True
        self.close()

    def getCheckBoxState(self, checkbox, attr):
        self.returnData[attr] = checkbox.isChecked()
        return

    def browseColor(self, attr, ledit=None):
        color = QColorDialog.getColor()
        if color.isValid():
            if ledit:
                ledit.setText(color.name())
            return color.name()
        return None

    def browseFileDir(self, attr, ledit=None, fileType=None):
        if fileType == "Directory":
            myFileDir = QtWidgets.QFileDialog.getExistingDirectory(
                self, caption="Locate directory", directory=self.defaultDir
            )
        else:
            myFileDir = QtWidgets.QFileDialog.getOpenFileName(
                self, "Locate file", self.defaultDir
            )[0]
        if myFileDir:
            self.returnData[attr] = myFileDir
        if ledit:
            ledit.setText(myFileDir)


class tabbedWindow(QTabWidget):
    def __init__(self, parent=None):
        super(tabbedWindow, self).__init__(parent)

    def add(self, title, minHeight=240):
        tab = QWidget()
        tab.setMinimumHeight(minHeight)
        self.addTab(tab, title)
        box = gui.widgetBox(tab)
        tab.setLayout(QVBoxLayout())
        tab.layout().addWidget(self.getScrollArea(box))
        return self.getLeditLayout(box)

    def addBox(self, title, minHeight=240):
        tab = QWidget()
        tab.setMinimumHeight(minHeight)
        self.addTab(tab, title)
        box = gui.widgetBox(tab)
        tab.setLayout(QVBoxLayout())
        tab.layout().addWidget(self.getScrollArea(box))
        return box, self.getLeditLayout(box)

    def getScrollArea(self, box):
        scroll_area = QScrollArea(verticalScrollBarPolicy=Qt.ScrollBarAlwaysOn)
        scroll_area.setWidget(box)
        scroll_area.setWidgetResizable(True)
        return scroll_area

    def getLeditLayout(self, box):
        layout = QtGui.QGridLayout()
        layout.setSpacing(5)
        setattr(layout, "nextRow", 1)
        box.layout().addLayout(layout)
        return layout


class DragAndDropList(QtGui.QListWidget):
    # overloads the Drag and dropEvents to emit code
    itemMoved = pyqtSignal(int, int)  # Old index, new index, item

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


class WidgetItem:
    # each widget item has:
    # A set of gui widgets that correspond to the name and each parameter
    # A set of state vectors corresponding to the different elements that make up the item

    def __init__(self, guiElements=OrderedDict(), states={}, values={}):
        # guielements has parameter as key and qwidget as value
        self.guiElements = guiElements
        self.states = states
        self.values = values
        if self.states and self.guiElements:
            for (key, element) in self.guiElements.items():
                element.setState(self.states[key])
        elif self.guiElements:
            self.blankState()

    def printValue(self):
        output = "{} : {{".format(self.guiElements["name"].getValue())
        for key, element in self.guiElements.items():
            if key != "name" and element.getValue() is not None:
                output += " {} : {}, ".format(key, element.getValue())
        output = output[:-2] + " }"
        return output

    def getStateValue(self, states):
        name = None
        myDict = {}
        i = 0
        if states:
            for key, element in self.guiElements.items():
                if key == "name" or key == "Name":
                    name = element.getStateValue(states[i])
                else:
                    myDict[key] = element.getStateValue(states[i])
                i += 1
            return (name, myDict)
        return None

    def getState(self):
        serialStates = []
        for key, element in self.guiElements.items():
            self.states[key] = element.getState()
            serialStates.append(self.states[key])
        return serialStates

    def setState(self, serialStates):
        if serialStates is None:
            self.blankState()
        else:
            i = 0
            for key, element in self.guiElements.items():
                self.states[key] = serialStates[i]
                element.setState(self.states[key])
                i += 1

    def blankState(self):
        for key, element in self.guiElements.items():
            element.setState(None)
            self.states[key] = None


class OWWidgetBuilder(widget.OWWidget):
    name = "Widget Builder"
    description = (
        "Build a new widget from a set of bash commands and a Docker container"
    )
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

    def drawExec(self, layout=None):
        self.saveMode = QtGui.QComboBox()
        self.saveMode.addItem("Overwrite")
        self.saveMode.addItem("Merge")
        self.saveMode.addItem("Data")
        self.saveMode.setCurrentIndex(self.saveModeIndex)
        self.saveMode.currentIndexChanged.connect(
            lambda: self.onSaveModeChange(self.saveMode)
        )
        self.saveWidgetBtn = gui.button(None, self, "Save", callback=self.saveWidget)
        self.saveWidgetBtn.setStyleSheet(self.css)
        self.saveWidgetBtn.setFixedSize(50, 20)
        self.saveWidgetAsBtn = gui.button(
            None, self, "Save as", callback=self.saveWidgetAs
        )
        self.saveWidgetAsBtn.setStyleSheet(self.css)
        self.saveWidgetAsBtn.setFixedSize(70, 20)
        self.loadWidgetBtn = gui.button(None, self, "Load", callback=self.loadWidget)
        self.loadWidgetBtn.setStyleSheet(self.css)
        self.loadWidgetBtn.setFixedSize(70, 20)
        self.renameBtn = gui.button(None, self, "Rename", callback=self.renameWidget)
        self.renameBtn.setStyleSheet(self.css)
        self.renameBtn.setFixedSize(100, 20)
        # choose save mode

        saveLabel = QtGui.QLabel("Save mode:")
        saveLabel.setAlignment(Qt.AlignRight | Qt.AlignVCenter)
        self.saveMode.setCurrentIndex(self.saveModeIndex)
        box = QtGui.QHBoxLayout()
        box.addWidget(saveLabel)
        box.addWidget(self.saveMode)
        box.addWidget(self.saveWidgetBtn)
        box.addWidget(self.saveWidgetAsBtn)
        box.addWidget(self.loadWidgetBtn)
        if not self.newFlag:
            box.addWidget(self.renameBtn)
        box.addStretch(1)
        layout.addLayout(box)

    def onSaveModeChange(self, widget):
        self.saveModeIndex = widget.currentIndex()
        self.allAttrs["saveModeIndex"] = self.saveModeIndex

    def loadWidget(self, loadWidgetDir=None, loadNameCheck=False, startDir=None, newLoad=False):
        if not startDir:
            if self.widgetDir:
                startDir = self.widgetDir
            else:
                startDir = "/templates/Generic"
        if startDir == "__New":
            self.widgetName = self.getWidgetName()
            if not self.widgetName:
                return
            allAttrsFile = "/templates/Generic/Generic.attrs"
            allStatesFile = "/templates/Generic/Generic.states"
            self.allAttrs = self.unPickleData(allAttrsFile)
            self.allStates = self.unPickleData(allStatesFile)
            self.addGroupToStates()
            self.makeDefaultFiles()
            self.startWidget()
            return
        if not loadWidgetDir:
            loadWidgetDir = QtWidgets.QFileDialog.getExistingDirectory(
                self, caption="Choose widget to load", directory=startDir
            )
        if loadWidgetDir:
            widgetName = os.path.split(loadWidgetDir)[-1]
            allAttrsFile = "{}/{}.attrs".format(loadWidgetDir, widgetName)
            allStatesFile = "{}/{}.states".format(loadWidgetDir, widgetName)
            if os.path.exists(allAttrsFile):
                self.allAttrs = self.unPickleData(allAttrsFile)
            if os.path.exists(allStatesFile):
                self.allStates = self.unPickleData(allStatesFile)
                self.addGroupToStates()
            #copy the location of the icon file to the icon ledit
            if not (self.widgetDir) or (os.path.normpath(loadWidgetDir) != os.path.normpath(self.widgetDir)):
                if  not (self.widgetDir):
                    self.widgetDir='/templates/Generic'
                if hasattr(self,"widgetIconLedit"):
                    iconFiles=os.listdir(loadWidgetDir + "/icon")
                    if iconFiles:
                        self.widgetIconLedit.setState([None,None,[True,loadWidgetDir+"/icon/"+iconFiles[0]]])
                newDockerFiles="{}/Dockerfiles".format(loadWidgetDir)
                self.replaceDockerDir(newDockerFiles)
            # we need to make default files in case the default files change from versions
            if loadNameCheck or not self.widgetName:
                self.widgetName = self.getWidgetName()
            if not self.widgetName:
                return
            # new widget 
            self.makeDefaultFiles()
            self.checkIconFile(loadWidgetDir)
            if self.isDrawn:
                self.updateWidget()
            else:
                self.startWidget()
        return
    def getWidgetName(self):
        niceName = self.widgetName
        if not niceName:
            if "name" in self.allAttrs and self.allAttrs["name"]:
                niceName = self.allAttrs["name"]
            else:
                niceName = "generic{}".format(os.getpid())
            niceName = re.sub(r"[^\w\s]", "", niceName)
            niceName = re.sub(r"\s+", "_", niceName)
        nameDir, okPressed = QInputDialog.getText(
            self, "Widget Name ", "Enter widget name:", QLineEdit.Normal, niceName
        )
        if okPressed and nameDir:
            niceName = re.sub(r"[^\w\s]", "", nameDir)
            niceName = re.sub(r"\s+", "_", nameDir)
        else:
            # return what we started with
            return self.widgetName
        return niceName

    def renameWidget(self):
        renameData = {}
        newName, okPressed = QInputDialog.getText(
            self, "Rename widget ", "Enter new name:", QLineEdit.Normal, ""
        )
        if okPressed and newName:
            niceName = niceForm(newName, useDash=False)
            renameData["newName"] = niceName
            renameData["widgetName"] = self.widgetName
            renameData["category"] = niceForm(
                os.path.basename(os.path.dirname(self.widgetDir)), useDash=True
            )
            td = toolDockEdit.ToolDockEdit(renameData, canvasMainWindow=self.canvas)
            del td
        else:
            return
        # renames the widget definition, not the instance
        pass

    def buildData(self):
        myData = {}
        self.data["name"] = self.widgetName
        for attr in (
            "name",
            "description",
            "docker_image_name",
            "docker_image_tag",
            "priority",
            "icon",
            "inputs",
            "outputs",
            "volumes",
            "ports",
            "parameters",
            "command",
            "autoMap",
            "runDockerMap"
        ):
            if attr in self.data and self.data[attr]:
                myData[attr] = deepcopy(self.data[attr])
            else:
                myData[attr] = None
        # separate multiple commands and give default command
        if "command" in myData and myData["command"]:
            if "\n" in myData["command"]:
                myData["command"] = myData["command"].split("\n")
            else:
                myData["command"] = [myData["command"]]

        # add persistance -all for now - add option to for transient data later
        myData["persistentSettings"] = "all"

        # add required elements
        reqList = []
        if "parameters" in myData and myData["parameters"]:
            for pname, pvalue in myData["parameters"].items():
                if not pvalue["optional"]:
                    reqList.append(pname)
        myData["requiredParameters"] = reqList
        # add volume mappings
        if "volumes" in myData and myData["volumes"]:
            myMappings = []
            for attr, volume in myData["volumes"].items():
                myMappings.append(
                    {"conVolume": volume["containerVolume"], "attr": attr}
                )
            myData["volumeMappings"] = myMappings
            myData.pop("volumes", None)
        # add port mappings
        if "ports" in myData and myData["ports"]:
            myPortMappings = []
            for attr, pvalue in myData["ports"].items():
                myPortMappings.append(
                    {"containerPort": pvalue["containerPort"], "attr": attr}
                )
            myData["portMappings"] = myPortMappings
            myData.pop("ports", None)

        # replace text str with type(str)
        for pname in ("inputs", "outputs"):
            if pname in myData and myData[pname]:
                for key, myDict in myData[pname].items():
                    if myDict["type"] == "str":
                        myDict["type"] = type("str")

        # for now just keep default, flags,
        # also make all defaults false for booleans - will fix these kluges after cleaning up BwBase code
        if "parameters" in myData and myData["parameters"]:
            for key, myDict in myData["parameters"].items():
                newDict = {}

                if "default" in myDict and myDict["default"] is not None:
                    sys.stderr.write(
                        "key is {} default is {}\n".format(key, myDict["default"])
                    )
                    # make sure that these are the correct type
                    if myDict["type"] == "bool":
                        if myDict["default"] == "False":
                            newDict["default"] = False
                        elif myDict["default"] == "True":
                            newDict["default"] = True
                        else:
                            raise Exception(
                                "{} is boolean - default values must be True or False not {}".format(
                                    pname, myDict["default"]
                                )
                            )
                    # check for lists - these are treated as column delimited
                    elif "type" in myDict and "list" in myDict["type"]:
                        reader = csv.reader([myDict["default"]], skipinitialspace=True)
                        for row in reader:
                            newDict["default"] = row
                            sys.stderr.write(
                                "default list is type{} {}\n".format(
                                    type(newDict["default"]), newDict["default"]
                                )
                            )
                    elif "type" in myDict and myDict["type"] == "int":
                        newDict["default"] = int(myDict["default"])
                        sys.stderr.write(
                            "key is {} new default is {}\n".format(
                                key, myDict["default"]
                            )
                        )
                    elif "type" in myDict and (
                        myDict["type"] == "double" or myDict["type"] == "float"
                    ):
                        newDict["default"] = float(myDict["default"])
                    elif "type" in myDict and (
                        myDict["type"] == "str"
                        or myDict["type"] == "file"
                        or myDict["type"] == "directory"
                    ):
                        newDict["default"] = str(myDict["default"])
                    else:
                        newDict["default"] = myDict["default"]
                if "flag" in myDict:
                    newDict["flag"] = myDict["flag"]
                # arguments are the same as having a null flag value
                if "argument" in myDict and myDict["argument"]:
                    newDict["argument"] = True
                if "label" in myDict and myDict["label"] is not None:
                    newDict["label"] = myDict["label"]
                if "type" in myDict and myDict["type"]:
                    newDict["type"] = myDict["type"]
                if "env" in myDict and myDict["env"]:
                    newDict["env"] = myDict["env"]
                if "group" in myDict and myDict["group"]:
                    newDict["group"] = myDict["group"]
                myData["parameters"][key] = newDict
        return myData

    def saveJson(self):
        myData = self.buildData()
        # get json filename
        dataJ = jsonpickle.encode(myData)
        jSaveFile = QtWidgets.QFileDialog.getSaveFileName(
            self,
            "QFileDialog.getSaveFileName()",
            "",
            "All Files (*);;json Files (*.json)",
        )[0]
        if jSaveFile:
            sys.stderr.write("save file is {}\n".format(jSaveFile))
            with open(jSaveFile, "w") as f:
                f.write(dataJ)
            f.close()

    def saveWidget(self):
        if not self.widgetDir or not self.widgetName:
            self.saveWidgetAs()
            return
        self.pickleWidget()
        qm = QtGui.QMessageBox
        title = "Save {}".format(self.widgetName)
        ret = qm.question(
            self,
            title,
            "Saved widget to {} - Reload?".format(self.widgetDir),
            qm.Yes | qm.No,
        )
        if ret == qm.No:
            return
        if self.canvas:
            self.canvas.reload_current()
        return

    def saveWidgetAs(self):
        qm = QtGui.QMessageBox
        self.widgetName = self.getWidgetName()
        if not self.widgetName:
            return
        outputDir = QtWidgets.QFileDialog.getExistingDirectory(
            self,
            caption="Choose directory to save the widget in",
            directory=self.defaultDir,
        )
        if outputDir:
            myWidgetDir = "{}/{}".format(outputDir, self.widgetName)
            if os.path.exists(myWidgetDir) and not (
                self.widgetDir and os.path.samefile(myWidgetDir, self.widgetDir)
            ):
                # same directory is occupied
                # ask permission to nuke it
                ret = qm.question(
                    self,
                    "",
                    "{}/{} exists - OverWrite ?".format(outputDir, self.widgetName),
                    qm.Yes | qm.No,
                )
                if ret == qm.No:
                    return
                os.system("cd {} && rm {}/* -rf ".format(outputDir, self.widgetName))
            self.widgetDir = outputDir + "/" + self.widgetName
            self.allAttrs["name"] = self.widgetName
            self.makeDefaultFiles()
            self.pickleWidget()
            self.nameLabel.setText(self.widgetName)
            self.setWindowTitle(self.widgetName + ":Definition")

            drawerName = os.path.basename(outputDir)
            if drawerName in self.directoryList:
                initPy = "/biodepot/{}/__init__.py".format(drawerName)
                workflowPath = ""
                if os.path.islink(initPy):
                    workflowPath = os.path.dirname(os.path.realpath(initPy))
                if not workflowPath or workflowPath == outputDir:
                    destLink = "/biodepot/{}/OW{}.py".format(
                        drawerName, self.widgetName
                    )
                    pythonFile = glob("{}/*.py".format(self.widgetDir))[0]
                    os.system("ln -sf {} {}".format(pythonFile, destLink))

                # make a symlink to the
            title = "Save {}".format(self.widgetName)
            message = "Saved widget to {}".format(self.widgetDir)
            ret = qm.question(
                self,
                title,
                "Saved widget to {} - Reload?".format(self.widgetDir),
                qm.Yes | qm.No,
            )
            if ret == qm.No:
                return
            if self.canvas:
                self.canvas.reload_current()
        return
        
    def checkIconFile(self,loadWidgetDir):
        if 'icon' in self.allStates and self.allStates['icon'][2][1] and os.path.isfile(self.allStates['icon'][2][1]):
            loadIconFile = self.allStates['icon'][2][1]
            if self.widgetDir:
                if os.path.realpath(self.widgetDir) != os.path.realpath(loadWidgetDir):
                    os.system("rm {}/icon/*".format(self.widgetDir))
                    os.system("cp {} {}/icon/.".format(loadIconFile,self.widgetDir))
                    self.resetIconFile()
                else:
                    self.resetIconFile()
        else:
            loadIconFile = os.listdir(loadWidgetDir + "/icon")[0]
            if self.widgetDir:
                if os.path.realpath(self.widgetDir) != os.path.realpath(loadWidgetDir):
                    breakpoint(message="cp {}/icon/{} {}/icon/.".format(loadWidgetDir,loadIconFile,self.widgetDir))
                    os.system("cp {}/icon/{} {}/icon/.".format(loadWidgetDir,loadIconFile,self.widgetDir)) 
                self.resetIconFile()
            else:
                self.allStates['icon'][2][1]="{}/icon/{}".format(loadWidgetDir,loadIconFile) 

    def resetIconFile(self):
        if "icon" in self.allAttrs:
            self.allAttrs["icon"]=None
        if "icon" in self.allStates:
            self.allStates["icon"][2][1]=''
        return

        
    def makeDefaultFiles(self):
        if not self.widgetDir:
            return
        #this routine makes the default files/directories that might not be present in earlier versions
        iconDir = "{}/icon".format(self.widgetDir)
        dockerFiles="{}/Dockerfiles".format(self.widgetDir)
        os.system("mkdir -p {} ".format(iconDir))
        os.system("mkdir -p {} ".format(dockerFiles))
        #add iconFile if none exists
        if not os.listdir("{}/icon".format(self.widgetDir)):
            copyfile(defaultIconFile, iconDir + "/" + os.path.basename(defaultIconFile))
        

    def pickleWidget(self):
        myData = self.buildData()
        Path(self.widgetDir).mkdir(parents=True, exist_ok=True)
        outputWidget = "{}/{}.py".format(self.widgetDir, self.widgetName)
        # check if data only
        if self.saveMode.currentIndex() < 2:
            if self.saveMode.currentIndex() == 1 and os.path.exists(self.outputWidget):
                mergeWidget(None, outputWidget, self.widgetName, inputData=myData)
            else:
                createWidget(None, outputWidget, self.widgetName, inputData=myData)
        allStatesFile = "{}/{}.states".format(self.widgetDir, self.widgetName)
        allAttrsFile = "{}/{}.attrs".format(self.widgetDir, self.widgetName)
        self.pickleData(self.allStates, allStatesFile)
        self.pickleData(self.allAttrs, allAttrsFile)

    def getDefaultDir(self):
        defaultDir = "/root"
        if os.path.exists("/data"):
            defaultDir = "/data"
        return defaultDir

    def pickleData(self, data, filename, jsonFlag=True):
        if jsonFlag:
            myJdata = jsonpickle.encode(data)
            with open(filename, "w") as f:
                f.write(myJdata)
            f.close()
        else:
            with open(filename, "wb") as f:
                pickle.dump(data, f, protocol=pickle.HIGHEST_PROTOCOL)
            f.close()

    def unPickleData(self, filename, jsonFlag=True):
        if jsonFlag:
            sys.stderr.write("opening file {}\n".format(filename))
            try:
                with open(filename, "r") as f:
                    data = jsonpickle.decode(f.read())
                f.close()
            except Exception as e:
                with open(filename, "rb") as f:
                    data = pickle.load(f)
                f.close()
        else:
            with open(filename, "rb") as f:
                data = pickle.load(f)
            f.close()
        return data

    def __init__(self, widgetID=None, canvasMainWindow=None):
        super().__init__()
        # Directories it toolDock
        self.directoryList = (
            str(
                os.popen("""grep -oP 'packages=\["\K[^"]+' /biodepot/setup.py""").read()
            )
        ).split()
        # minimum sizes
        self.canvas = canvasMainWindow
        self.controlArea.setMinimumWidth(600)
        self.controlArea.setMinimumHeight(400)

        self.css = """
        QPushButton {background-color: #1588c5; color: white; height: 20px; border: 1px solid black; border-radius: 2px;}
        QPushButton:hover {background-color: #1555f5; }
        QPushButton:hover:pressed { background-color: #1588c5; color: black; border-style: inset; border: 1px solid white} 
        QPushButton:disabled { background-color: lightGray; border: 1px solid gray; } 
        """
        self.addRemoveCSS = """
        QPushButton {background-color: lightBlue; color: white; height: 20px; border: 1px solid black; border-radius: 2px;}
        QPushButton:hover {background-color: blue; }
        QPushButton:hover:pressed { background-color: lightBlue; color: black; border-style: inset; border: 1px solid white} 
        QPushButton:disabled { background-color: white; border: 1px solid gray; } 
        """
        self.browseCSS = """
        QPushButton {background-color: rgba(30,30,200,128); color: white; height: 20px; border: 1px solid black; border-radius: 2px;}
        QPushButton:hover {background-color: #1555f5; }
        QPushButton:hover:pressed { background-color: #1588c5; color: black; border-style: inset; border: 1px solid white} 
        QPushButton:disabled { background-color: lightGray; border: 1px solid gray; } 
        """
        # self.setStyleSheet(self.css)
        self.browseIcon = QtGui.QIcon("/icons/bluefile.png")
        self.addIcon = QtGui.QIcon("/icons/add.png")
        self.removeIcon = QtGui.QIcon("/icons/remove.png")
        self.submitIcon = QtGui.QIcon("/icons/submit.png")
        self.reloadIcon = QtGui.QIcon("/icons/reload.png")
        self.outputWidget = ""
        self.defaultDir = self.getDefaultDir()
        self.widgetDir = None
        self.widgetName = None
        self.isDrawn = False
        self.containerID = None
        self.saveModeIndex = 0
        self.newFlag = False
        if widgetID == "__New":
            qm = QtGui.QMessageBox
            ret = qm.question(self, "", "Create new widget?", qm.Yes | qm.No)
            if ret == qm.No:
                raise ValueError("User cancelled")
            else:
                self.newFlag = True
                tmp = tempfile.mkdtemp()
                ret = qm.question(
                    self, "", "Use existing widget as template?", qm.Yes | qm.No
                )
                if ret == qm.Yes:
                    self.loadWidget(startDir="/widgets")
                else:
                    self.loadWidget(startDir="__New")
                if self.widgetName:
                    self.setWindowTitle(self.widgetName + ":Definition")
                else:
                    raise ValueError("no widget name given")
        else:
            widgetSplit = widgetID.split(".")
            self.widgetName = widgetSplit[-1][2:]
            widgetPy = os.path.realpath("/biodepot/{}.py".format("/".join(widgetSplit)))
            self.widgetDir = os.path.dirname(widgetPy)
            self.setWindowTitle(self.widgetName + ":Definition")
            sys.stderr.write(
                "widgetPy is {} widgetDir is {} widgetName is {}\n".format(
                    widgetPy, self.widgetDir, self.widgetName
                )
            )
            self.loadWidget(loadWidgetDir=self.widgetDir, loadNameCheck=False)
            self.defaultDir = os.path.dirname(self.widgetDir)

    def clearLayout(self, layout):
        if layout != None:
            while layout.count():
                child = layout.takeAt(0)
                if child.widget() is not None:
                    child.widget().deleteLater()
                elif child.layout() is not None:
                    self.clearLayout(child.layout())

    def startWidget(self):
        self.isDrawn = True
        self.setWindowTitle(self.widgetName + ":Definition")
        for attr in (
            "name",
            "description",
            "docker_image_name",
            "docker_image_tag",
            "priority",
            "icon",
            "inputs",
            "outputs",
            "volumes",
            "parameters",
            "command",
            "autoMap",
            "runDockerMap",
            "buildCommand",
        ):
            if attr in self.allAttrs:
                self.data[attr] = self.allAttrs[attr]
            else:
                self.data[attr] = None
                self.allAttrs[attr] = None
        if "saveModeIndex" in self.allAttrs:
            if self.allAttrs["saveModeIndex"] is None:
                self.allAttrs["saveModeIndex"] = 0
            self.saveModeIndex = self.allAttrs["saveModeIndex"]
        self.tabs = tabbedWindow()
        self.controlArea.layout().addWidget(self.tabs)
        self.setStyleSheet(":disabled { color: #282828}")
        # self.clearLayout(self.controlArea.layout())
        # self.controlArea.layout().addWidget(self.scroll_area)

        # controlBox = gui.vBox(self.generalBox)
        # self.drawExec(box=controlBox)ne)
        # requiredBox = gui.widgetBox(self.generalBox, "Widget entries")
        # draw Ledits for the frequired elements
        leditGeneralLayout = self.tabs.add("General")
        self.nameLabel = QtGui.QLabel("Name: " + self.widgetName)
        leditGeneralLayout.addWidget(self.nameLabel, leditGeneralLayout.nextRow, 0)
        leditGeneralLayout.nextRow = leditGeneralLayout.nextRow + 1
        for pname in ["description", "docker_image_name", "docker_image_tag"]:
            self.drawLedit(pname, layout=leditGeneralLayout)
        self.drawLedit("priority", layout=leditGeneralLayout)
        # file entry for icon
        self.widgetIconLedit=self.drawLedit(
            "icon", layout=leditGeneralLayout, addBrowseButton=True, fileType=None
        )

        # listwidgets for inputs and outputs
        # define widgetItems for the different widgetLists
        # top level widgets are drawXXX - these can have multiple substituents
        # lower level widgets are makeXXX - these can also have multiple substituents

        self.drawIListWidget("inputs", layout=self.tabs.add("Inputs"))
        self.drawOListWidget("outputs", layout=self.tabs.add("Outputs"))
        self.drawVolumeListWidget("volumes", layout=self.tabs.add("Volumes"))
        self.drawPortsListWidget("ports", layout=self.tabs.add("Ports"))
        self.drawParamsListWidget("parameters", layout=self.tabs.add("Parameters"))
        self.drawCommand("command", layout=self.tabs.add("Command"))
        self.drawDocker("buildCommand", layout=self.tabs.add("Docker"))
        self.drawExec(self.controlArea.layout())

    def updateWidget(self):
        self.clearLayout(self.controlArea.layout())
        self.startWidget()

    def updateCheckBox(self, checkBox, widget=None):
        if checkBox.isEnabled():
            widget.setEnabled(checkBox.isChecked())

    def addListWidget(self, addBtn, qWidgetList, qWidgetItem):
        if qWidgetList.selectedItems():
            index = qWidgetList.row(qWidgetList.selectedItems()[0])
            qWidgetList.states[index] = qWidgetItem.getState()
            qWidgetList.values[index] = qWidgetList.selectedItems()[0].setText(qWidgetItem.printValue())
            qWidgetList.selectedItems()[0].setText(qWidgetItem.printValue())
        else:
            item = qWidgetItem.printValue()
            qWidgetList.addItem(item)
            qWidgetList.states.append(qWidgetItem.getState())
            qWidgetList.values.append(item)
            sys.stderr.write("adding state {}\n".format(qWidgetItem.getState()))
            qWidgetItem.blankState()
        qWidgetList.updateAllStates()

    def removeListWidget(self, removeBtn, qWidgetList, qWidgetItem):
        if qWidgetList.selectedItems():
            for item in qWidgetList.selectedItems():
                del qWidgetList.states[qWidgetList.row(item)]
                del qWidgetList.values[qWidgetList.row(item)]
                qWidgetList.takeItem(qWidgetList.row(item))
        if not qWidgetList.count():
            removeBtn.setEnabled(False)
        qWidgetItem.blankState()
        qWidgetList.updateAllStates()

    def qwUpdateAllStates(self, pname, qWidgetList, qWidgetItem):
        serialState = []
        widgetList = []
        for state in qWidgetList.states:
            serialState.append(state)
            widgetList.append(qWidgetItem.getStateValue(state))
            sys.stderr.write(
                "append state {} value {}\n".format(
                    state, qWidgetItem.getStateValue(state)
                )
            )
        if serialState:
            self.allStates[pname] = serialState
            self.data[pname] = OrderedDict(widgetList)
        else:
            self.allStates[pname] = None
            self.data[pname] = qWidgetItem.getStateValue(None)

    def addGroupToStates(self):
        # for older files
        if "parameters" in self.allStates and  self.allStates["parameters"]:
            for serialState in self.allStates["parameters"]:
                if len(serialState) < 9:
                    serialState.append(serialState[7])
                    serialState[7] = [[True, False], [True, "group"], [False, ""]]

    def qwInitAllStates(self, pname, qWidgetList, qWidgetItem):
        widgetList = []
        if pname in self.allStates and self.allStates[pname]:
            for serialState in self.allStates[pname]:
                sys.stderr.write("init pname {} state {}\n".format(pname, serialState))
                qWidgetItem.blankState()
                qWidgetItem.setState(serialState)
                item = qWidgetItem.printValue()
                qWidgetList.addItem(item)
                qWidgetList.states.append(qWidgetItem.getState())
                qWidgetList.values.append(item)
                widgetList.append(qWidgetItem.getStateValue(serialState))
        if widgetList:
            self.data[pname] = OrderedDict(widgetList)
        else:
            self.data[pname] = None
        qWidgetItem.blankState()

    def onListWidgetSelect(self, qWidgetList, addBtn, removeBtn, qWidgetItem):
        if qWidgetList.selectedItems():
            if len(qWidgetList.selectedItems()) > 1:
                addBtn.setEnabled(False)
            else:
                item = qWidgetList.selectedItems()[0]
                myState = qWidgetList.states[qWidgetList.row(item)]
                sys.stderr.write("loading {}\n".format(myState))
                qWidgetItem.setState(myState)
                addBtn.setEnabled(True)
        else:
            removeBtn.setEnabled(False)
            addBtn.setEnabled(True)


    def synchronizeStates(self,qlist):
        #the states and values in a list will be out of sync after a drag and drop event
        #could do shuffle with temp but simpler just to mark and make changes
        indices=[]
        states=[]        
        for i in range (qlist.count()):
            itemText=qlist.item(i).text()
            if itemText != qlist.values[i]:
                indices.append(i)
                states.append(qlist.states[qlist.values.index(itemText)])
        for i in range (len(indices)):
            index=indices[i]
            qlist.states[index]=states[i]
            qlist.values[index]=qlist.item(index).text()
            
    def onItemMoved(self, oldRow, newRow, qWidgetList):
        sys.stderr.write("oldRow is {} newRow is {} len {}\n".format(oldRow, newRow,len(qWidgetList.selectedItems())))
        #algorithm works if newRow > oldRow
        sliceLen=len(qWidgetList.selectedItems())
        if newRow != oldRow:
            self.synchronizeStates(qWidgetList)
        qWidgetList.updateAllStates()


    def drawCommand(self, pname, layout=None):
        labelTextBox = self.makeTextBox(pname, label="Enter command:")
        self.initAllStates(pname, labelTextBox)
        labelTextBox.textBox.textChanged.connect(
            lambda: self.updateAllStates(pname, labelTextBox, labelTextBox.getState())
        )
        # layout.addWidget(labelTextBox.label,layout.nextRow,0)
        layout.addWidget(labelTextBox.textBox, layout.nextRow, 1, 1, 4)
        layout.nextRow = layout.nextRow + 1

    def drawDocker(self, pname, layout=None):
        dockerFileLedit = self.makeLedit(
            "dockerFileToAdd", "Enter file to add", label="Add Dockerfile:"
        )
        layout.addWidget(dockerFileLedit.label, layout.nextRow, 0)
        layout.addWidget(dockerFileLedit.ledit, layout.nextRow, 1, 1, 1)
        callback = self.browseFileDir
        browseFileBtn = gui.button(
            None,
            self,
            "",
            callback=lambda: callback(
                "dockerFileToAdd", ledit=dockerFileLedit.ledit, fileType=None
            ),
            autoDefault=True,
            width=20,
            height=20,
        )
        browseFileBtn.setIcon(self.browseIcon)
        browseFileBtn.setStyleSheet(self.browseCSS)
        layout.addWidget(browseFileBtn, layout.nextRow, 2)
        addFileBtn = gui.button(
            None,
            self,
            "",
            callback=lambda: self.addDockerFile(dockerFileLedit.ledit),
            autoDefault=False,
        )
        addFileBtn.setIcon(self.addIcon)
        addFileBtn.setStyleSheet(self.addRemoveCSS)
        layout.addWidget(addFileBtn, layout.nextRow, 3)
        layout.nextRow += 1

        dockerDirLedit = self.makeLedit(
            "dockerDirToAdd", "Enter directory to add", label="Add directory:"
        )
        layout.addWidget(dockerDirLedit.label, layout.nextRow, 0)
        layout.addWidget(dockerDirLedit.ledit, layout.nextRow, 1, 1, 1)
        callback = self.browseFileDir
        browseDirBtn = gui.button(
            None,
            self,
            "",
            callback=lambda: callback(
                "dockerDirToAdd", ledit=dockerDirLedit.ledit, fileType="Directory"
            ),
            autoDefault=True,
            width=20,
            height=20,
        )
        browseDirBtn.setIcon(self.browseIcon)
        browseDirBtn.setStyleSheet(self.browseCSS)
        layout.addWidget(browseDirBtn, layout.nextRow, 2)
        addDirBtn = gui.button(
            None,
            self,
            "",
            callback=lambda: self.addDockerDir(dockerDirLedit.ledit),
            autoDefault=False,
        )
        addDirBtn.setIcon(self.addIcon)
        addDirBtn.setStyleSheet(self.addRemoveCSS)
        layout.addWidget(addDirBtn, layout.nextRow, 3)
        layout.nextRow += 1
        # self.drawLedit('Add file to Dockerfiles',layout=layout,addBrowseButton=True, fileType=None)
        # self.drawLedit('Add directory to Dockerfiles',layout=layout,addBrowseButton=True, fileType='Directory')

        # addDateCb=self.makeCheckBox ('addBuildDate','Add date to docker tag',default=True,persist=True,track=True)
        # containerIDLabel=QtGui.QLabel('Container ID: {}'.format(self.containerID))
        imageBuilderLabel = QtGui.QLabel("Launch Image Builder")
        imageBuilderBtn = gui.button(
            None, self, "Launch", callback=self.startImageBuilder
        )
        imageBuilderBtn.setStyleSheet(self.css)
        imageBuilderBtn.setFixedSize(80, 20)
        updateLabel = QtGui.QLabel("Clear Dockerfiles")
        updateBtn = gui.button(None, self, "Clear", callback=self.clearDockerfiles)
        updateBtn.setStyleSheet(self.css)
        updateBtn.setFixedSize(80, 20)
        buildCommandBox = self.makeTextBox(pname, label="Docker build command:")
        self.initAllStates(pname, buildCommandBox)
        buildCommandBox.textBox.textChanged.connect(
            lambda: self.updateAllStates(
                pname, buildCommandBox, buildCommandBox.getState()
            )
        )
        layout.addWidget(updateLabel, layout.nextRow, 0)
        layout.addWidget(updateBtn, layout.nextRow, 1)
        layout.addWidget(imageBuilderLabel, layout.nextRow + 1, 0)
        layout.addWidget(imageBuilderBtn, layout.nextRow + 1, 1)
        # layout.addWidget(addDateCb,layout.nextRow+2,0)
        # layout.addWidget(containerIDLabel,layout.nextRow+3,0)
        # layout.addWidget(addDateCb,layout.nextRow+4,0)
        layout.addWidget(buildCommandBox.label, layout.nextRow + 2, 0)
        layout.addWidget(buildCommandBox.textBox, layout.nextRow + 2, 1, 1, 4)

    def addDockerFile(self, ledit):
        if not hasattr(ledit, "text"):
            return
        src = ledit.text().strip()
        if src and os.path.isfile(src):
            dest = self.widgetDir + "/" + "Dockerfiles/" + os.path.basename(src)
            copyfile(src, dest)
            qm = QtGui.QMessageBox
            title = "Added Dockerfile"
            message = "Added {}".format(src)
            qm.information(self, title, message, QtGui.QMessageBox.Ok)
            
    def replaceDockerDir(self,newDockerfiles):
        if newDockerfiles and os.path.isdir(newDockerfiles):
            currentDockerfiles = self.widgetDir + "/" + "Dockerfiles"
            if os.path.exists(currentDockerfiles):
                qm = QtGui.QMessageBox
                ret = qm.question(
                    self,
                    "Confirm delete",
                    "The directory {} exits Delete?".format(currentDockerfiles),
                    qm.Yes | qm.No,
                )
                if ret == qm.No:
                    return False
                shutil.rmtree(currentDockerfiles)
            shutil.copytree(newDockerfiles, currentDockerfiles)
            qm = QtGui.QMessageBox
            title = "Copied Dockerfiles"
            message = "Added Dockerfiles {}".format(newDockerfiles)
            qm.information(self, title, message, QtGui.QMessageBox.Ok)
        else:
            breakpoint(message="{} not found".format(newDockerfiles))
            
    def addDockerDir(self, ledit):
        if not hasattr(ledit, "text"):
            return
        src = ledit.text().strip()
        if src and os.path.isdir(src):
            dest = self.widgetDir + "/" + "Dockerfiles/" + os.path.basename(src)
            if os.path.isdir(dest):
                qm = QtGui.QMessageBox
                ret = qm.question(
                    self,
                    "Confirm delete",
                    "The directory {} exits Delete?".format(os.path.basename(dest)),
                    qm.Yes | qm.No,
                )
                if ret == qm.No:
                    return
                shutil.rmtree(dest)
            shutil.copytree(src, dest)
            qm = QtGui.QMessageBox
            title = "Added Directory"
            message = "Added {}".format(src)
            qm.information(self, title, message, QtGui.QMessageBox.Ok)

    def startImageBuilder(self):
        dockerDir = self.widgetDir + "/" + "Dockerfiles"
        widget = OWImageBuilder.OWImageBuilder(dockerDir)
        widget.showNormal()
        widget.raise_()
        widget.activateWindow()

    def clearDockerfiles(self):
        dockerDir = self.widgetDir + "/" + "Dockerfiles"
        if not os.listdir(dockerDir):
            return
        qm = QtGui.QMessageBox
        ret = qm.question(
            self,
            "Confirm delete",
            "Are you sure you want to remove all files in {} ?".format(dockerDir),
            qm.Yes | qm.No,
        )
        if ret == qm.No:
            return False
        os.system("cd {} && rm -rf * ".format(dockerDir))
        title = "Cleared docker files"
        message = "Cleared docker files"
        qm.information(self, title, message, QtGui.QMessageBox.Ok)
        return True

    def makeTextBox(self, attr, label):
        box = QHBoxLayout()
        textLabel = None
        if label:
            textLabel = QtGui.QLabel(label)
            textLabel.setAlignment(Qt.AlignTop)
        textBox = QtGui.QPlainTextEdit()
        textBox.setStyleSheet(":disabled { color: #282828}")
        setattr(box, "textBox", textBox)
        setattr(box, "label", textLabel)
        setattr(box, "getValue", textBox.toPlainText())
        setattr(
            box, "getStateValue", lambda state: self.getTextBoxStateValue(box, state)
        )
        setattr(box, "getState", lambda: self.getTextBoxState(box))
        setattr(box, "setState", lambda state: self.setTextBoxState(box, state))
        return box

    def getTextBoxState(self, widget):
        return [widget.textBox.isEnabled(), widget.textBox.toPlainText()]

    def getTextBoxStateValue(self, widget, state):
        if state and state[0] and state[1]:
            return state[1]
        return None

    def setTextBoxState(self, widget, state):
        if state is None:
            widget.label.setEnabled(True)
            widget.textBox.setEnabled(True)
            widget.textBox.clear()
        else:
            widget.textBox.setEnabled(state[0])
            widget.textBox.setPlainText(state[1])

    # workhorse widget
    def makeLedit(
        self, leditAttr, text=None, label=None, addCheckBox=False, initialValue=None
    ):
        leditLabel = None
        checkBox = None
        if label:
            leditLabel = QtGui.QLabel(label)
        ledit = QtGui.QLineEdit(self)
        ledit.setClearButtonEnabled(True)
        ledit.setPlaceholderText(text)
        ledit.setStyleSheet(":disabled { color: #282828}")
        if not initialValue:
            ledit.clear()
        if addCheckBox:
            checkBox = QtGui.QCheckBox(None, self)
            checkBox.stateChanged.connect(lambda: self.updateCheckBox(checkBox, ledit))
            ledit.setEnabled(checkBox.isChecked())
        box = QHBoxLayout()
        if checkBox:
            box.addWidget(checkBox)
        if leditLabel:
            box.addWidget(leditLabel)
            box.addWidget(ledit)
        setattr(box, "ledit", ledit)
        setattr(box, "label", leditLabel)
        setattr(box, "checkBox", checkBox)
        setattr(box, "getValue", lambda: self.getLeditValue(ledit))
        setattr(
            box,
            "getStateValue",
            lambda state: self.getLeditStateValue(checkBox, leditLabel, ledit, state),
        )
        setattr(
            box, "getState", lambda: self.getLeditState(checkBox, leditLabel, ledit)
        )
        setattr(
            box,
            "setState",
            lambda state: self.setLeditState(checkBox, leditLabel, ledit, state),
        )
        return box

    def getLeditValue(self, ledit):
        if ledit.isEnabled():
            if ledit.text() is None:
                return ""
            return ledit.text()
        return None

    def getLeditStateValue(self, checkBox, label, ledit, state):
        sys.stderr.write("getledit state vector is {}\n".format(state))
        if state is None:
            return None
        if checkBox:
            if not state[0][1]:
                return None
        return state[2][1]

    def getLeditState(self, checkBox, label, ledit):
        # can't serialize qt objects
        # so we just concatenate the data as a list instead of a dict
        state = []
        if checkBox:
            state.append([checkBox.isEnabled(), checkBox.isChecked()])
        else:
            state.append([None, None])
        if label:
            state.append([label.isEnabled(), label.text()])
        else:
            state.append([None, None])
        if ledit:
            state.append([ledit.isEnabled(), ledit.text()])
        else:
            state.append([None, None])
        return state

    def setLeditState(self, checkBox, label, ledit, state):
        if state is None:
            # initialize it
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
        if checkBox is not None and state[0] is not None:
            checkBox.setEnabled(state[0][0])
            checkBox.setChecked(state[0][1])

        if label is not None and state[1] is not None:
            sys.stderr.write(
                "labelEnable is {} labelState is {}\n".format(state[1][0], state[1][1])
            )
            label.setEnabled(state[1][0])
            label.setText(state[1][1])

        if ledit is not None and state[2] is not None:
            ledit.setEnabled(state[2][0])
            if state[2][1]:
                ledit.setText(state[2][1])
            else:
                ledit.clear()

    def makeComboBox(self, pname, label, elements):
        comboBoxLabel = QtGui.QLabel(label)
        comboBox = QtGui.QComboBox()
        comboBox.addItems(elements)
        comboBox.currentIndex = 0
        box = QHBoxLayout()
        box.addWidget(comboBoxLabel)
        box.addWidget(comboBox)
        setattr(box, "getValue", lambda: self.getComboValue(comboBox))
        setattr(
            box, "getStateValue", lambda state: self.getComboStateValue(comboBox, state)
        )
        setattr(box, "getState", lambda: self.getComboState(comboBox, comboBoxLabel))
        setattr(
            box,
            "setState",
            lambda state: self.setComboState(comboBox, comboBoxLabel, state),
        )
        return box

    def getComboValue(self, comboBox):
        if comboBox.isEnabled():
            return comboBox.currentText()
        return None

    def getComboStateValue(self, comboBox, state):
        if not state or not state[0]:
            return None
        elif state[0][0]:  # check if enabled
            index = state[0][1]
            if index >= 0:
                allItems = state[0][2]
                return allItems[index]
            else:
                return None
        return None

    def getComboState(self, comboBox, label):
        state = []
        if comboBox:
            allItems = [comboBox.itemText(i) for i in range(comboBox.count())]
            enabled = comboBox.isEnabled()
            index = comboBox.findText(comboBox.currentText())
            state.append([enabled, index, allItems])
        else:
            state.append(None)
        if label:
            state.append([label.isEnabled(), label.text()])
        else:
            state.append(None)
        return state

    def setComboState(self, comboBox, label, state):
        if state is None:
            # intialize
            comboBox.setEnabled(True)
            comboBox.setCurrentIndex(0)
            return None
        if state[0] is None:
            comboBox = None
        else:
            comboBox.setEnabled(state[0][0])
            # can have different states stored if the comboBox values are changed
            self.adjustComboBoxState(state, comboBox)
            comboBox.setCurrentIndex(state[0][1])
        if state[1] is None:
            label = None
        else:
            label.setEnabled(state[1][0])
            label.setText(state[1][1])

    def adjustComboBoxState(self, state, comboBox):
        allItems = [comboBox.itemText(i) for i in range(comboBox.count())]
        if allItems == state[0][2]:
            return
        # find the value of the current element
        currentStateElement = state[0][2][state[0][1]]
        for i in range(len(allItems)):
            if allItems[i] == currentStateElement:
                state[0][2] = allItems
                state[0][1] = i
                return
        # no match then set current element to 0
        state[0][2] = allItems
        state[0][1] = 0

    def makeCheckBox(self, attr, label, default=False, persist=False, track=False):
        # checkbox for binary options
        # not used as part of other elements
        if attr not in self.data:
            self.data[attr] = default
        checkBox = QtGui.QCheckBox(label, self)
        setattr(
            checkBox, "getValue", lambda: checkBox.isEnabled() and checkBox.isChecked()
        )
        setattr(
            checkBox,
            "getStateValue",
            lambda state: self.getCheckBoxStateValue(checkBox, state, default),
        )
        setattr(checkBox, "getState", lambda: self.getCheckBoxState(checkBox))
        setattr(
            checkBox,
            "setState",
            lambda state: self.setCheckBoxState(checkBox, state, default),
        )
        if persist:
            self.initAllStates(attr, checkBox)
            checkBox.stateChanged.connect(
                lambda: self.updateAllStates(attr, checkBox, checkBox.getState())
            )
        if track:
            if not persist:
                checkBox.stateChanged.connect(
                    lambda: self.updateAllStates(attr, checkBox, checkBox.getState())
                )
        return checkBox

    def getCheckBoxStateValue(self, checkBox, state, default):
        if not state:
            return default
        if state[0] and state[1]:
            return True
        return False

    def getCheckBoxState(self, checkBox):
        state = []
        if checkBox:
            state = [checkBox.isEnabled(), checkBox.isChecked()]
        else:
            state = None
        return state

    def setCheckBoxState(self, checkBox, state, default=False):
        if state is None:
            # initialize
            checkBox.setEnabled(True)
            checkBox.setChecked(default)
            return state
        if state[0] is None:
            checkBox = None
        else:
            checkBox.setEnabled(state[0])
            checkBox.setChecked(state[1])

    def makeListWidget(self, pname, boxEdit):
        # setup boxEdit
        # logic is handled by add remove buttons
        if pname not in self.data:
            self.data[pname] = None
        boxEdit = DragAndDropList(self)
        boxEdit.setDragDropMode(QtWidgets.QAbstractItemView.InternalMove)
        boxEdit.setSelectionMode(QtWidgets.QAbstractItemView.ExtendedSelection)
        boxEdit.setStyleSheet(":disabled { color: #282828}")
        boxEdit.setMinimumHeight(60)
        return boxEdit
    
    def makeListWidgetUnit(
        self, pname, layout=None, lineWidgets=None, otherWidgets=None, elementsPerLine=5
    ):
        # widgets is ordered dict of widgets managed by listWidget
        boxEdit = self.makeListWidget(pname, None)
        setattr(boxEdit, "states", [])
        setattr(boxEdit, "values", [])
        lineItem = WidgetItem(OrderedDict(lineWidgets), {})
        setattr(
            boxEdit,
            "updateAllStates",
            lambda: self.qwUpdateAllStates(pname, boxEdit, lineItem),
        )
        self.qwInitAllStates(pname, boxEdit, lineItem)
        setattr(boxEdit, "oldRow", boxEdit.currentRow())
        # buttons
        addBtn = gui.button(
            None,
            self,
            "",
            callback=lambda: self.addListWidget(addBtn, boxEdit, lineItem),
            autoDefault=False,
        )
        removeBtn = gui.button(
            None,
            self,
            "",
            callback=lambda: self.removeListWidget(removeBtn, boxEdit, lineItem),
            autoDefault=False,
        )
        removeBtn.setEnabled(bool(boxEdit.selectedItems()))
        boxEdit.itemSelectionChanged.connect(
            lambda: removeBtn.setEnabled(bool(boxEdit.selectedItems()))
        )
        boxEdit.itemSelectionChanged.connect(
            lambda: self.onListWidgetSelect(boxEdit, addBtn, removeBtn, lineItem)
        )
        boxEdit.itemMoved.connect(
            lambda oldRow, newRow: self.onItemMoved(oldRow, newRow, boxEdit)
        )

        # make button box for buttons so that they don't expand in linelayout
        # setmaximumwidth seems to be ignored
        # need blank field and high columns to push the buttons to the right and keep their size
        buttonBoxLabel = QtGui.QLabel("")
        buttonBox = QtGui.QGridLayout()
        buttonBox.addWidget(buttonBoxLabel, 0, 0, 1, 20)
        buttonBox.addWidget(addBtn, 0, 25)
        buttonBox.addWidget(removeBtn, 0, 26)
        addBtn.setIcon(self.addIcon)
        addBtn.setStyleSheet(self.addRemoveCSS)
        removeBtn.setIcon(self.removeIcon)
        removeBtn.setStyleSheet(self.addRemoveCSS)
        # layout
        # basic layout is other widgets on top with filesBox below and the lineItem underneath

        filesBoxLeditLayout = QtGui.QVBoxLayout()
        # add to the main parameters box
        myBox = gui.vBox(None)
        myBox.layout().addLayout(filesBoxLeditLayout)
        # label=QtGui.QLabel(pname+':')
        # label.setAlignment(Qt.AlignTop)
        # layout.addWidget(label,layout.nextRow,0)
        # otherwidgets here
        if otherWidgets:
            for widget in otherWidgets:
                layout.addWidget(widget, layout.nextRow, 1, 1, 2)
                layout.nextRow = layout.nextRow + 1
        # box with list here
        layout.addWidget(myBox, layout.nextRow, 1, 1, 2)
        layout.nextRow = layout.nextRow + 1
        # input with line layout here
        lineLayout = QtGui.QGridLayout()
        col = 0
        row = 1
        for item in lineWidgets:
            layout = item[1]
            # The Widgets are actually layouts except for the checkboxes
            if col >= elementsPerLine:
                row += 1
                col = 0
            if isinstance(layout, QtGui.QCheckBox):
                lineLayout.addWidget(layout, row, col)
            else:
                lineLayout.addLayout(layout, row, col)
            col += 1
        lineLayout.addLayout(buttonBox, row, col)
        # lineLayout.addWidget(addBtn,row,col+1)
        # lineLayout.addWidget(removeBtn,row,col+2)

        # now add the two layouts to the bigBox layout
        filesBoxLeditLayout.addWidget(boxEdit)
        filesBoxLeditLayout.addLayout(lineLayout)

    def drawIListWidget(self, pname, layout=None):
        nameBox = self.makeLedit(pname + "nameLedit", "Enter name", "Name")
        callbackBox = self.makeLedit(
            pname + "callbackLedit", "Enter callback", "callback", addCheckBox=True
        )
        comboBox = self.makeComboBox(
            pname, "Type:", ["str", "dict", "Orange.data.Table"]
        )
        widgetList = [("name", nameBox), ("callback", callbackBox), ("type", comboBox)]
        self.makeListWidgetUnit(pname, layout=layout, lineWidgets=widgetList)

    def drawOListWidget(self, pname, layout=None):
        nameBox = self.makeLedit(pname + "nameLedit", "Enter name", "Name")
        defaultBox = self.makeLedit(
            pname + "defaultLedit", "Enter default", "Default value", addCheckBox=True
        )
        comboBox = self.makeComboBox(
            pname, "Type:", ["str", "dict", "Orange.data.Table"]
        )
        widgetList = [("name", nameBox), ("default", defaultBox), ("type", comboBox)]
        self.makeListWidgetUnit(pname, layout=layout, lineWidgets=widgetList)

    def drawVolumeListWidget(self, pname, layout=None):
        nameBox = self.makeLedit(pname + "Name", "Enter variable", "Local volume variable")
        volumeBox = self.makeLedit(
            pname + "volumeLedit", "Enter volume", "Container volume"
        )
        widgetList = [("name", nameBox), ("containerVolume", volumeBox)]
        autoMapCb = self.makeCheckBox(
            "autoMap",
            "Pass local to Bwb mappings to container",
            default=True,
            persist=True,
            track=True,
        )
        runDockerMapCb = self.makeCheckBox(
            "runDockerMap",
            "Pass mappings for launching containers",
            default=False,
            persist=True,
            track=True,
        )
        self.makeListWidgetUnit(
            pname, layout=layout, lineWidgets=widgetList, otherWidgets=[autoMapCb, runDockerMapCb]
        )

    def drawPortsListWidget(self, pname, layout=None):
        nameBox = self.makeLedit(
            pname + "Name", "Enter variable name", "Host port variable"
        )
        containerPortBox = self.makeLedit(
            pname + "containerPort", "Enter container port", "container port"
        )
        widgetList = [("name", nameBox), ("containerPort", containerPortBox)]
        self.makeListWidgetUnit(pname, layout=layout, lineWidgets=widgetList)

    def drawParamsListWidget(self, pname, layout=None):
        nameBox = self.makeLedit(pname + "nameLedit", "Enter name", "Name")
        flagBox = self.makeLedit(
            pname + "flagLedit", "Enter flag", "flag", addCheckBox=True
        )
        labelBox = self.makeLedit(
            pname + "labelLedit", "Enter label", "label", addCheckBox=True
        )
        envBox = self.makeLedit(
            pname + "envLedit", "Enter ENV variable", "env", addCheckBox=True
        )
        defaultBox = self.makeLedit(
            pname + "defaultLedit", "Enter default", "default", addCheckBox=True
        )
        groupBox = self.makeLedit(
            pname + "groupLedit", "Enter group", "group", addCheckBox=True
        )
        setattr(self, "optional", False)
        setattr(self, "argument", False)
        optionalCb = self.makeCheckBox("optional", "Optional")
        argumentCb = self.makeCheckBox("argument", "Argument")
        # connect argument checkbox to disabling the flag checkbox
        flagBox.checkBox.stateChanged.connect(
            lambda: (not flagBox.checkBox.isChecked()) or (argumentCb.setChecked(False))
        )
        argumentCb.stateChanged.connect(
            lambda: (not argumentCb.isChecked())
            or (flagBox.checkBox.setChecked(False) or flagBox.ledit.clear())
        )
        comboBox = self.makeComboBox(
            pname,
            "Type:",
            [
                "str",
                "file",
                "file list",
                "directory",
                "directory list",
                "bool",
                "bool list",
                "text list",
                "int",
                "int list",
                "double",
                "double list",
                "patternQuery",
            ],
        )
        widgetList = [
            ("name", nameBox),
            ("type", comboBox),
            ("flag", flagBox),
            ("argument", argumentCb),
            ("env", envBox),
            ("label", labelBox),
            ("default", defaultBox),
            ("group", groupBox),
            ("optional", optionalCb),
        ]
        self.makeListWidgetUnit(pname, layout=layout, lineWidgets=widgetList)

    def drawLedit(
        self, pname, layout=None, addBrowseButton=False, fileType=None, callback=None
    ):
        # make labeledLedit combo layout
        labelLedit = self.makeLedit(pname, "Enter " + pname, label=pname + ":")
        self.initAllStates(pname, labelLedit)
        labelLedit.ledit.textChanged.connect(
            lambda: self.updateAllStates(pname, labelLedit, labelLedit.getState())
        )
        layout.addWidget(labelLedit.label, layout.nextRow, 0)
        if addBrowseButton:
            layout.addWidget(labelLedit.ledit, layout.nextRow, 1, 1, 1)
            if callback is None:
                callback = self.browseFileDir
            button = gui.button(
                None,
                self,
                "",
                callback=lambda: callback(
                    pname, ledit=labelLedit.ledit, fileType=fileType
                ),
                autoDefault=True,
                width=19,
                height=19,
            )
            button.setIcon(self.browseIcon)
            button.setStyleSheet(self.browseCSS)
            layout.addWidget(button, layout.nextRow, 2)
        else:
            layout.addWidget(labelLedit.ledit, layout.nextRow, 1, 1, 2)
        layout.nextRow += 1
        return labelLedit

    def updateAllStates(self, pname, widget, state):
        self.allStates[pname] = state
        self.data[pname] = widget.getStateValue(state)
        sys.stderr.write("updating all - pname {} state {}\n".format(pname, state))

    def initAllStates(self, attr, widget):
        if attr in self.allStates:
            widget.setState(self.allStates[attr])
        else:
            widget.setState(None)
            self.allStates[attr] = widget.getState()
        self.data[attr] = widget.getStateValue(self.allStates[attr])

    def browseFileDir(self, attr, ledit=None, fileType=None):
        self.defaultDir = self.getDefaultDir()
        if fileType == "Directory":
            myFileDir = QtWidgets.QFileDialog.getExistingDirectory(
                self, caption="Locate directory", directory=self.defaultDir
            )
        else:
            myFileDir = QtWidgets.QFileDialog.getOpenFileName(
                self, "Locate file", self.defaultDir
            )[0]
        if myFileDir:
            self.data[attr] = myFileDir
            self.defaultDir = myFileDir
        if ledit:
            ledit.setText(myFileDir)

    def updateCheckBoxLayout(self, cb, layout):
        if cb.isChecked():
            for i in reversed(range(layout.count())):
                widget = layout.itemAt(i).widget()
                if isinstance(widget, QCheckBox):
                    widget.setChecked(False)
                    widget.setEnabled(False)
                elif isinstance(widget, QLabel):
                    continue
                elif isinstance(widget, QtWidgets.QLineEdit):
                    widget.clear()
                    widget.setEnabled(False)
                else:
                    widget.setEnabled(False)
        else:
            for i in reversed(range(layout.count())):
                widget = layout.itemAt(i).widget()
                if isinstance(widget, QCheckBox):
                    widget.setChecked(False)
                    widget.setEnabled(True)
                elif isinstance(widget, QLabel):
                    continue
                elif isinstance(widget, QtWidgets.QLineEdit):
                    widget.setEnabled(False)
                else:
                    widget.setEnabled(True)
