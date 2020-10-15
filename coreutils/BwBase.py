import os
import re
import sys
import logging
import jsonpickle
import functools
import multiprocessing
import glob
from itertools import zip_longest
from OWWidgetBuilder import tabbedWindow
from ServerUtils import IterateDialog
from functools import partial
from pathlib import Path
from AnyQt.QtCore import QThread, pyqtSignal, Qt
from Orange.widgets import widget, gui, settings
from DockerClient import DockerClient, PullImageThread, ConsoleProcess
from PyQt5 import QtWidgets, QtGui, QtCore
from PyQt5.QtCore import *
from PyQt5.QtGui import *
from PyQt5.QtWidgets import *
from time import sleep


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


def breakpoint(title=None, message=None):
    return
    QtGui.QMessageBox.warning(title,'',message)


def getJsonName(filename, widgetName):
    widgetPy = os.path.realpath(filename)
    widgetDir = os.path.dirname(widgetPy)
    return "{}/{}.json".format(widgetDir, widgetName)


def getIconName(filename, iconFile):
    widgetPy = os.path.realpath(filename)
    widgetDir = os.path.dirname(widgetPy)
    return "{}/icon/{}".format(widgetDir, iconFile)


class ScrollMessageBox(QMessageBox):
    def __init__(self, l, *args, **kwargs):
        QMessageBox.__init__(self, *args, **kwargs)
        scroll = QScrollArea(self)
        scroll.setWidgetResizable(True)
        self.content = QWidget()
        scroll.setWidget(self.content)
        lay = QVBoxLayout(self.content)
        if l:
            for item in l:
                lay.addWidget(QLabel(item, self))
                self.layout().addWidget(scroll, 0, 0, 1, self.layout().columnCount())
            self.setStyleSheet("QScrollArea{min-width:300 px; min-height: 400px}")
        else:
            lay.addWidget(QLabel("No Match", self))
            self.layout().addWidget(scroll, 0, 0, 1, self.layout().columnCount())
            self.setStyleSheet("QScrollArea{min-width:50 px;}")


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


# for selecting multiple directories
class getExistingDirectories(QtWidgets.QFileDialog):
    def __init__(self, *args):
        super(getExistingDirectories, self).__init__(*args)
        self.setOption(self.DontUseNativeDialog, True)
        self.setFileMode(self.Directory)
        self.setOption(self.ShowDirsOnly, True)
        self.findChildren(QtWidgets.QListView)[0].setSelectionMode(
            QtWidgets.QAbstractItemView.ExtendedSelection
        )
        self.findChildren(QtWidgets.QTreeView)[0].setSelectionMode(
            QtWidgets.QAbstractItemView.ExtendedSelection
        )
        self.findChildren(QtWidgets.QTreeView)[0].setSelectionMode(
            QtWidgets.QAbstractItemView.ExtendedSelection
        )


class BwbHelperFunctions:
    def __init__(self):
        self.initialize = {}
        self.disableGUI = {}
        self.generate = {}
        self.receive = {}

    def add(self, attr, functionType, function):
        try:
            myDict = getattr(self, functionType)
            myDict[attr] = function
        except Exception as e:
            sys.stderr.write(
                "unable to add function {} for attribute {} to helper functions\n".format(
                    functionType, attr
                )
            )
            sys.exit()

    def function(self, attr, functionType):
        if hasattr(self, functionType):
            myDict = getattr(self, functionType)
            if attr in myDict:
                return myDict[attr]
        else:
            sys.stderr.write(
                "no functionType {} for attr {}".format(functionType, attr)
            )
        return None


class BwbGridLayout:
    # adds methods to keep track of rows and columns
    def __init__(self, spacing=5, startCol=0, startRow=0):
        self.css = """
        QPushButton {background-color: #1588c5; color: white; height: 20px; border: 1px solid black; border-radius: 2px;}
        QPushButton:hover {background-color: #1555f5; }
        QPushButton:hover:pressed { background-color: #1588c5; color: black; border-style: inset; border: 1px solid white} 
        QPushButton:disabled { background-color: lightGray; border: 1px solid gray; } 
        """
        self._layout = QtGui.QGridLayout()
        self._layout.setSpacing(spacing)
        self.nextCol = startCol
        self.nextRow = startRow

    def addWidget(self, widget, space=None, width=None, height=1, linefeed=None):
        if space is not None:
            self.nextCol += space
        if width is not None:
            self._layout.addWidget(widget, self.nextRow, self.nextCol, width, height)
            self.nextCol += width
        else:
            self._layout.addWidget(widget, self.nextRow, self.nextCol)
            self.nextCol += 1
        if linefeed is not None:
            self.nextRow += linefeed
            self.nextCol = 0

    def addSpace(self, space=1):
        self.nextCol += space

    def addLinefeed(self, lf=1, startCol=0):
        self.nextCol = startCol
        self.nextRow += lf

    def layout(self):
        return self._layout


class BwbGuiElements:
    # keep track of List of Gui elements
    # add support for initialization callbacks
    # keep track of values

    # need to check for tuple elements from Orange gui
    def disableElement(self, element):
        if isinstance(element, tuple):
            for g in element:
                g.setDisabled(True)
        else:
            element.setDisabled(True)

    def enableElement(self, element):
        if isinstance(element, tuple):
            for g in element:
                g.setEnabled(True)
        else:
            element.setEnabled(True)

    def __init__(self, required=None, active=None):
        self._enableCallbacks = {}
        self._dict = {}
        self.required = required
        self.active = {}

    def add(self, attr, guiElement, enableCallback=None, updateCallback=None):
        if attr not in self._dict:
            self._dict[attr] = []
        self._dict[attr].append(guiElement)
        if enableCallback is not None:
            self._enableCallbacks[attr] = enableCallback

    def addList(self, attr, guiElements, enableCallback=None, updateCallback=None):
        if attr not in self._dict:
            self._dict[attr] = guiElements
        else:
            self._dict[attr].extend(guiElements)
        if enableCallback is not None:
            self._enableCallbacks[attr] = enableCallback

    def disable(
        self, attr, ignoreCheckbox=False, disableFtn=None
    ):  # gray out to disable
        # if the disabling is due to an input - there may be a custom partial disable function
        if attr in self._dict:
            if disableFtn:
                disableFtn()
                return
            if ignoreCheckbox:
                for g in self._dict[attr]:
                    if not isinstance(g, QtWidgets.QCheckBox):
                        sys.stderr.write("ignore cb, {} disabling {}\n".format(attr, g))
                        self.disableElement(g)
            else:
                for g in self._dict[attr]:
                    sys.stderr.write("{} disabling {}\n".format(attr, g)),
                    self.disableElement(g)
            return True
        return False

    def clear(self, attr, activate=False):
        # will activate first cb when activate is True
        # otherwise
        if not attr in self._dict or not self._dict[attr]:
            return
        i = 0
        if isinstance(self._dict[attr][0], QtWidgets.QCheckBox):
            self._dict[attr][0].setCheckState(activate)
            i = 1
        # check for all ledits which are the only elements with a clear function
        # some of the orange gui elements are tuples
        while i < len(self._dict[attr]):
            g = self._dict[attr][i]
            if isinstance(g, tuple):
                for element in g:
                    if isinstance(element, QtWidgets.QLineEdit):
                        element.clear()
            else:
                if isinstance(g, QtWidgets.QLineEdit):
                    sys.stderr.write("clearing line edit {}\n".format(g))
                    g.clear()
            i = i + 1

    def enable(self, attr, value):
        sys.stderr.write("checking attr {}\n".format(attr))
        clearLedit = False
        if value is None or value is "":
            clearLedit = True
        if attr in self._dict:
            sys.stderr.write("found attr in dict {}\n".format(attr))
            if attr in self._enableCallbacks:
                sys.stderr.write("enable callback for {}\n".format(attr))
                self._enableCallbacks[attr](value, clearLedit)
            else:
                sys.stderr.write("enable for {}\n".format(attr))
                for g in self._dict[attr]:
                    sys.stderr.write("enable element {}\n".format(g))
                    self.enableElement(g)
        return True

    def disableAll(self):
        for attr in self._dict.keys():
            self.disable(attr)

    def reenableAll(self, OWself):
        for attr in self._dict.keys():
            if not OWself.inputConnections.isConnected(attr):
                self.enable(attr, getattr(OWself, attr))


class ConnectionDict:
    def __init__(self, inputDict):
        self._dict = {}
        self._dict = (
            inputDict
        )  # we do want the name not a copy - i.e. the inputDict should change

    def add(self, slot, connectionId=None):
        if connectionId is None:
            return
        if slot in self._dict and connectionId not in self._dict[slot]:
            self._dict[slot].append(connectionId)
        else:
            self._dict[slot] = [connectionId]
        sys.stderr.write("Adding {} to connection {}\n".format(slot, self._dict[slot]))

    def remove(self, slot, connectionId=None):
        sys.stderr.write("removing {} to connection {}\n".format(slot, connectionId))
        if slot in self._dict:
            if connectionId is None:
                del self._dict[slot]
            else:
                try:
                    self._dict[slot].remove(connectionId)
                except ValueError:
                    # this is probably due to an addition of a link that did not register
                    # so on error - delete the slot
                    sys.stderr.write(
                        "exception raise on remove {} - delete slot instead".format(
                            connectionId
                        )
                    )
                    del self._dict[slot]

        # need to update the volumes
        else:
            sys.stderr.write(
                "no matching slot for connection {}\n".format(connectionId)
            )

    def isConnected(self, slot, connectionId=None):
        if self._dict is None:
            return False
        if slot in self._dict:
            sys.stderr.write(
                "slot found in connection check with contents {} and ID {}\n".format(
                    self._dict[slot], connectionId
                )
            )
            if connectionId:
                if connectionId in self._dict[slot]:
                    return True
                else:
                    return False
            elif not self._dict[slot]:
                del self._dict[slot]
                sys.stderr.write("deleted slot\n")
                return False
            else:
                return True
        return False

    def isSet(self, slot):  # same as is connected but checks for non null value
        if self._dict is None:
            return False
        if (slot in self._dict) and (self._dict[slot] is not None):
            return True
        return False


class OWBwBWidget(widget.OWWidget):
    serversFile = "/biodepot/serverSettings.json"
    dockerClient = DockerClient("unix:///var/run/docker.sock", "local")
    defaultFileIcon = QtGui.QIcon("/icons/bluefile.png")
    browseIcon = QtGui.QIcon("/icons/bluefile.png")
    addIcon = QtGui.QIcon("/icons/add.png")
    removeIcon = QtGui.QIcon("/icons/remove.png")
    submitIcon = QtGui.QIcon("/icons/submit.png")
    reloadIcon = QtGui.QIcon("icons/reload.png")
    useScheduler = settings.Setting(False, schema_only=True)
    pset = functools.partial(settings.Setting, schema_only=True)
    nWorkers = pset(1)
    iterateSettings = pset({})
    iterate = pset(False)

    # Initialization
    def __init__(self, image_name, image_tag):
        super().__init__()
        self.helpers = BwbHelperFunctions()

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
        self.useTestMode = False
        self.jobRunning = False
        self.saveBashFile = None
        if not hasattr(self, "iterateSettings"):
            setattr(self, "iterateSettings", {})
            self.iterateSettings["iteratedAttrs"] = []
            self.iterateSettings["widgetThreads"] = 1
            self.iterateSettings["data"] = {}

        self.inputConnections = ConnectionDict(self.inputConnectionsStore)
        self._dockerImageName = image_name
        self._dockerImageTag = image_tag
        # setting the qbutton groups

        # drawing layouts for gui
        # file directory
        self.filesBoxLayout = QtGui.QVBoxLayout()
        self.fileDirRequiredLayout = BwbGridLayout()
        self.fileDirOptionalLayout = BwbGridLayout()
        self.fileDirScheduleLayout = QtGui.QGridLayout()
        # lineEdits
        self.leditRequiredLayout = BwbGridLayout()
        self.leditOptionalLayout = BwbGridLayout()
        self.leditOptionalLayout = BwbGridLayout()

        # keep track of gui elements associated with each attribute
        self.bgui = BwbGuiElements()

        # For compatibility if triggers are not being kept
        if not hasattr(self, "triggerReady"):
            self.triggerReady = {}

    def getQBGroups(self):
        QBList = []
        self.QButtonHash = {}
        if not hasattr(self, "QButtons"):
            self.QButtons = []
        if not hasattr(self, "QBGroups"):
            self.QBGroups = {}
        if not ("parameters" in self.data):
            return
        if  self.data["parameters"] is None:
            return
        for pname in self.data["parameters"]:
            pvalue = self.data["parameters"][pname]
            if "group" in pvalue and pvalue["group"]:
                myGroup = pvalue["group"]
                if pvalue["group"] not in QBList:
                    self.QBList.append(myGroup)
                    self.QButtons.append(QButtonGroup())
                self.QButtonHash[pname] = self.QButtons[self.QBList.index(myGroup)]

    def initVolumes(self):
        # initializes container volumes
        # also initializes button groups - to avoid changing old widgets we do that here
        self.getQBGroups()
        if "volumeMappings" in self.data and self.data["volumeMappings"]:
            for mapping in self.data["volumeMappings"]:
                bwbVolAttr = mapping["attr"]
                if not hasattr(self, bwbVolAttr):
                    setattr(self, bwbVolAttr, None)

    # Override the send function to handle the test flag
    def send(self, attr, *args, **kwargs):
        if self.useTestMode:
            super().send(attr, *args, test=self.useTestMode, **kwargs)
        else:
            super().send(attr, *args, **kwargs)

    # Drawing the GUI
    def drawGUI(self):
        self.tabs = tabbedWindow()
        self.controlArea.layout().addWidget(self.tabs)
        self.setStyleSheet(":disabled { color: #282828}")

        if "requiredParameters" in self.data and self.data["requiredParameters"]:
            self.requiredBox, requiredLayout = self.tabs.addBox(
                "Required entries", minHeight=160
            )
            self.requiredBox.layout().addLayout(self.fileDirRequiredLayout.layout())
            setattr(self.fileDirRequiredLayout.layout(), "added", True)
            self.requiredBox.layout().addLayout(self.leditRequiredLayout.layout())
            setattr(self.leditRequiredLayout.layout(), "added", True)
            self.drawRequiredElements()

        if self.findOptionalElements():
            self.optionalBox, optionalLayout = self.tabs.addBox(
                "Optional entries", minHeight=160
            )
            # self.optionalBox = gui.widgetBox(self.topBox, "Optional parameters")
            self.optionalBox.layout().addLayout(self.fileDirOptionalLayout.layout())
            setattr(self.fileDirOptionalLayout.layout(), "added", True)
            self.optionalBox.layout().addLayout(self.leditOptionalLayout.layout())
            setattr(self.leditOptionalLayout.layout(), "added", True)
            self.drawOptionalElements()

            self.scheduleBox,scheduleLayout=self.tabs.addBox('Scheduler',minHeight=160)
            self.scheduleBox.layout().addLayout(self.fileDirScheduleLayout.layout())
            self.drawScheduleElements()

        # disable connected elements
        for i in self.inputs:
            attr = i.name
            if self.inputConnections.isConnected(attr):
                self.bgui.disable(attr)

        # make a box for the console and console control - not abs necessary but it might help with future org and with the size hinting system
        self.consoleBox, consoleLayout = self.tabs.addBox("Console", minHeight=160)
        
        # self.consoleBox = gui.widgetBox(self.controlArea)
        consoleControlLayout = BwbGridLayout()
        self.consoleBox.layout().addLayout(consoleControlLayout.layout())
        setattr(consoleControlLayout, "added", True)
        self.drawConsoleControl(box=self.consoleBox, layout=consoleControlLayout)
        self.console = QtGui.QTextEdit()
        self.console.setReadOnly(True)
        pal = QtGui.QPalette()
        pal.setColor(QtGui.QPalette.Base, Qt.black)
        pal.setColor(QtGui.QPalette.Text, Qt.green)
        self.console.setPalette(pal)
        self.console.setAutoFillBackground(True)
        self.consoleBox.layout().addWidget(self.console)
        controlBox = QtGui.QVBoxLayout()
        self.pConsole = ConsoleProcess(
            console=self.console, finishHandler=self.onRunFinished
        )
        self.drawExec(box=self.controlArea.layout())
        self.checkTrigger()

    def drawConsoleControl(self, box=None, layout=None):

        layout.addWidget(QtGui.QLabel("Console: "))
        pname = "saveLog"
        pvalue = {"type": "directory", "label": "AutoLog"}
        if not hasattr(self, pname):
            setattr(self, pname, None)

        self.btnConsoleClear = gui.button(
            None, self, "Clear", callback=self.clearConsole
        )
        self.btnConsoleClear.setStyleSheet(self.css)
        self.btnConsoleClear.setFixedSize(60, 20)
        self.btnConsoleSave = gui.button(None, self, "Save", callback=self.saveConsole)
        self.btnConsoleSave.setStyleSheet(self.css)
        self.btnConsoleSave.setFixedSize(60, 20)
        self.cboThreadNumber = QtGui.QComboBox()
        self.cboThreadNumber.setMaximumWidth(100)
        self.cboThreadNumberPopulate()
        self.displayThread=self.cboThreadNumber.currentIndex()
        self.cboThreadNumber.currentIndexChanged.connect(self.cboThreadNumberUpdate)
        layout.addWidget(self.cboThreadNumber)
        layout.addWidget(self.btnConsoleClear)
        layout.addWidget(self.btnConsoleSave)
        # self.drawFileDirElements(pname, pvalue, box=box,layout=layout, addCheckbox=True)
        # keep track of buttons separately because the drawFileDir routine has its own callback for enabling the elements
        # should fix this at some point to check which elements have been dealt with using a dict or work out a recursive scheme with elements
        btnPname = pname + "Btns"
        if not hasattr(self, btnPname):
            setattr(self, btnPname, None)
        self.bgui.add(btnPname, self.btnConsoleClear)
        self.bgui.add(btnPname, self.btnConsoleSave)

    def cboThreadNumberUpdate(self):
        sys.stderr.write('changed cboThreadNumber index from {} to {}\n'.format(self.displayThread,self.cboThreadNumber.currentIndex()))
        self.displayThread= self.cboThreadNumber.currentIndex()
        self.pConsole.changeThreadNumber(self.displayThread)
        
    def cboThreadNumberPopulate(self):
        if hasattr(self,'cboThreadNumber'):
            self.cboThreadNumber.clear()            
            self.cboThreadNumber.addItem('All threads')
            if self.iterate:
                self.cboThreadNumber.setEnabled(True)
                for i in range (1, self.nWorkers+1):
                    self.cboThreadNumber.addItem('Thread {}'.format(i))
            else:
                self.cboThreadNumber.setEnabled(False)

    def clearConsole(self):
        self.console.clear()

    def saveConsole(self):
        defaultDir = "/root"
        if os.path.exists("/data"):
            defaultDir = "/data"
        try:
            fileName = QtWidgets.QFileDialog.getSaveFileName(
                self,
                "QFileDialog.getSaveFileName()",
                "",
                "Text files (*.txt);;All Files (*)",
            )[0]
            if fileName:
                with open(fileName, "w") as myFile:
                    myFile.write(str(self.console.toPlainText()))
                    myFile.close()
        except Exception as e:
            return

        # consoleControlLayout.addWidget(outputLabel,0,0)

    def drawElements(self, elementList, isOptional=False):
        for pname in elementList:
            if not ("parameters" in self.data) or not (
                pname in self.data["parameters"]
            ):
                continue
            pvalue = self.data["parameters"][pname]
            if (
                "gui" in pvalue
                and pvalue["gui"] != "file"
                and pvalue["gui"] != "directory"
            ) or (pvalue["type"] != "file" and pvalue["type"] != "directory"):
                continue
            if isOptional:
                self.drawFileDirElements(
                    pname,
                    pvalue,
                    box=self.optionalBox,
                    layout=self.fileDirOptionalLayout,
                    addCheckbox=True,
                )
            else:
                self.drawFileDirElements(
                    pname,
                    pvalue,
                    box=self.requiredBox,
                    layout=self.fileDirRequiredLayout,
                )

        for pname in elementList:
            if not ("parameters" in self.data) or not (
                pname in self.data["parameters"]
            ):
                continue
            pvalue = self.data["parameters"][pname]
            if ("gui" in pvalue and pvalue["gui"][-4:] != "list") or (
                pvalue["type"][-4:] != "list"
            ):
                continue
            sys.stderr.write(
                "drawing textBox for  pname {} pvalue {}\n".format(pname, pvalue)
            )
            if isOptional:
                self.drawTextBox(
                    pname,
                    pvalue,
                    box=self.optionalBox,
                    layout=self.fileDirOptionalLayout,
                    addCheckbox=True,
                )
            else:
                self.drawTextBox(
                    pname,
                    pvalue,
                    box=self.requiredBox,
                    layout=self.fileDirRequiredLayout,
                )

        for pname in elementList:
            if not ("parameters" in self.data) or not (
                pname in self.data["parameters"]
            ):
                continue
            pvalue = self.data["parameters"][pname]
            if ("gui" in pvalue and pvalue["gui"] != "Ledit") or (
                pvalue["type"] != "double"
                and pvalue["type"] != "str"
                and pvalue["type"] != type("text")
            ):
                sys.stderr.write("type is {} {}\n".format(pvalue["type"], type("text")))
                continue
            if isOptional:
                sys.stderr.write("type is {} {}\n".format(pvalue["type"], type("text")))
                self.drawLedit(
                    pname,
                    pvalue,
                    self.optionalBox,
                    layout=self.leditOptionalLayout,
                    addCheckbox=True,
                )
            else:
                self.drawLedit(
                    pname, pvalue, self.requiredBox, layout=self.leditRequiredLayout
                )

        for pname in elementList:
            if not ("parameters" in self.data) or not (
                pname in self.data["parameters"]
            ):
                continue
            pvalue = self.data["parameters"][pname]
            if ("gui" in pvalue and pvalue["gui"] != "Spin") or (
                pvalue["type"] != "int"
            ):
                continue
            if isOptional:
                self.drawSpin(pname, pvalue, self.optionalBox, addCheckbox=True)
            else:
                self.drawSpin(pname, pvalue, self.requiredBox)

        for pname in elementList:
            if not ("parameters" in self.data) or not (
                pname in self.data["parameters"]
            ):
                continue
            pvalue = self.data["parameters"][pname]
            if ("gui" in pvalue and pvalue["gui"] != "bool") or (
                pvalue["type"] != "bool"
            ):
                continue
            if isOptional:
                self.drawCheckbox(pname, pvalue, self.optionalBox)
            else:
                self.drawCheckbox(pname, pvalue, self.requiredBox)

        for pname in elementList:
            if not ("parameters" in self.data) or not (
                pname in self.data["parameters"]
            ):
                continue
            pvalue = self.data["parameters"][pname]
            if ("gui" in pvalue and pvalue["gui"] != "patternQuery") or (
                pvalue["type"] != "patternQuery"
            ):
                continue
            if isOptional:
                self.drawPatternQuery(
                    pname,
                    pvalue,
                    box=self.optionalBox,
                    layout=self.fileDirOptionalLayout,
                    addCheckbox=True,
                )
            else:
                self.drawPatternQuery(
                    pname,
                    pvalue,
                    box=self.requiredBox,
                    layout=self.fileDirRequiredLayout,
                )

    def drawRequiredElements(self):
        for pname in self.data["requiredParameters"]:
            if not ("parameters" in self.data) or not (
                pname in self.data["parameters"]
            ):
                continue
            pvalue = self.data["parameters"][pname]
            if not hasattr(self, pname):
                setattr(self, pname, None)
            if (getattr(self, pname) is None) and ("default" in pvalue):
                setattr(self, pname, pvalue["default"])
                sys.stderr.write(
                    "set {} to default {} with type {} - current value {} of type {}\n".format(
                        pname,
                        pvalue["default"],
                        type(pvalue["default"]),
                        getattr(self, pname),
                        type(getattr(self, pname)),
                    )
                )
            else:
                sys.stderr.write(
                    "{} has prevous value {} of type {}\n".format(
                        pname, getattr(self, pname), type(getattr(self, pname))
                    )
                )
        self.drawElements(self.data["requiredParameters"])

    def findOptionalElements(self):
        # checks that there are optional elements
        if not "parameters" in self.data or not self.data["parameters"]:
            return False
        for pname in self.data["parameters"]:
            if pname not in self.data["requiredParameters"]:
                return True
        return False

    def drawOptionalElements(self):
        optionalList = []
        for pname in self.data["parameters"]:
            if pname not in self.data["requiredParameters"]:
                pvalue = self.data["parameters"][pname]
                if not hasattr(self, pname):
                    setattr(self, pname, None)
                if (getattr(self, pname) is None) and ("default" in pvalue):
                    setattr(self, pname, pvalue["default"])
                    sys.stderr.write("default value is {}\n".format(pvalue["default"]))
                if not (pname in self.optionsChecked):
                    self.optionsChecked[pname] = False
                optionalList.append(pname)
        self.drawElements(optionalList, isOptional=True)

    def drawScheduleElements(self):
        with open(self.serversFile) as f:
            serverSettings = jsonpickle.decode(f.read())
        self.widgetThreads = 1
        self.widgetThreadsChecked = True
        self.IPs = []
        if "data" in serverSettings and serverSettings["data"]:
            for addr in serverSettings["data"]:
                self.IPs.append(addr)
        if not self.IPs:
            self.IPs = ["127.0.0.1"]
        self.schedulers = ["Default", "GKE", "AWS batch", "AWS lambda"]

        threadBox = QtGui.QHBoxLayout()
        scheduleBox = QtGui.QHBoxLayout()

        # have to wait until the OWxxx.py instance loads data from json before finding iterables - so don't move

        self.IPMenuItems = {}
        self.schedulerMenuItems = {}
        self.IPMenu = QtGui.QMenu(self)
        self.IPBtn = QtGui.QToolButton(self)
        self.IPBtn.setText("Servers")
        for attr in self.IPs:
            action = self.IPMenu.addAction(attr)
            action.setCheckable(True)
            action.setChecked(bool(attr in self.IPs))
            action.changed.connect(self.chooseIP)
            self.IPMenuItems[action] = attr
        self.IPBtn.setMenu(self.IPMenu)
        self.IPBtn.setPopupMode(QtGui.QToolButton.InstantPopup)

        schedulerBox = QHBoxLayout()
        self.schedulerLabel = QtGui.QLabel("Scheduler")
        self.schedulerComboBox = QtGui.QComboBox()
        self.schedulerComboBox.addItems(self.schedulers)
        self.schedulerComboBox.currentIndex = 0
        schedulerBox.addWidget(self.schedulerLabel)
        schedulerBox.addWidget(self.schedulerComboBox)

        cbLabel = QtGui.QLabel("Number of workers: ")
        if not hasattr(self, "nWorkers"):
            self.nWorkers = 1
        self.threadSpin = QSpinBox()
        self.threadSpin.setRange(1, 128)
        self.threadSpin.valueChanged.connect(self.updateThreadSpin)
        self.threadSpin.setValue(self.nWorkers)

        # only draw this if there are iterables

        if not hasattr(self, "iterate"):
            self.iterate = False

        self.iterateSettings["iterableAttrs"] = self.findIterables()
        if self.iterateSettings["iterableAttrs"]:
            self.iterateSettingsBtn = gui.button(
                None, self, "Settings", callback=self.setIteration
            )
            self.iterateSettingsBtn.setStyleSheet(self.css)
            self.iterateSettingsBtn.setFixedSize(80, 20)
            iterateCheckbox = gui.checkBox(None, self, "iterate", label="Iterate")
            self.iterateSettingsBtn.setEnabled(iterateCheckbox.isChecked())
            iterateCheckbox.stateChanged.connect(
                lambda: self.iterateSettingsBtn.setEnabled(iterateCheckbox.isChecked())
            )
            iterateCheckbox.stateChanged.connect(
                lambda: self.threadSpin.setEnabled(iterateCheckbox.isChecked())
            )
            iterateCheckbox.stateChanged.connect(
                self.cboThreadNumberPopulate
            )
            iterateBox = QtGui.QHBoxLayout()
            iterateBox.addWidget(iterateCheckbox)
            iterateBox.addWidget(self.iterateSettingsBtn)
            iterateBox.addStretch(1)

            self.threadSpin.setEnabled(iterateCheckbox.isChecked())
            self.fileDirScheduleLayout.addLayout(iterateBox, 0, 0)

        self.scheduleCheckbox = gui.checkBox(None, self, "useScheduler", label="")
        self.IPBtn.setEnabled(self.scheduleCheckbox.isChecked())
        self.schedulerLabel.setEnabled(self.scheduleCheckbox.isChecked())
        self.schedulerComboBox.setEnabled(self.scheduleCheckbox.isChecked())

        self.scheduleCheckbox.stateChanged.connect(
            lambda: self.IPBtn.setEnabled(self.scheduleCheckbox.isChecked())
        )
        self.scheduleCheckbox.stateChanged.connect(
            lambda: self.schedulerLabel.setEnabled(self.scheduleCheckbox.isChecked())
        )
        self.scheduleCheckbox.stateChanged.connect(
            lambda: self.schedulerComboBox.setEnabled(self.scheduleCheckbox.isChecked())
        )

        self.fileDirScheduleLayout.setAlignment(Qt.AlignTop)

        scheduleBox.addWidget(self.scheduleCheckbox)
        scheduleBox.addLayout(schedulerBox)
        scheduleBox.addWidget(self.IPBtn)
        scheduleBox.addStretch(1)

        threadBox.addWidget(cbLabel)
        threadBox.addWidget(self.threadSpin)
        threadBox.addStretch(1)

        self.fileDirScheduleLayout.addLayout(scheduleBox, 1, 0)
        self.fileDirScheduleLayout.addLayout(threadBox, 2, 0)
        
    def setIteration(self):
        iterateDialog = IterateDialog(self.iterateSettings)
        iterateDialog.exec_()
        self.iterateSettings = iterateDialog.iterateSettings

    def updateThreadSpin(self):
        self.nWorkers = self.threadSpin.value()
        self.cboThreadNumberPopulate()

    def drawCheckbox(self, pname, pvalue, box=None):
        # for booleans - their value is the same as the checkbox state
        sys.stderr.write(
            "drawCB pname {} pvalue {} label {}\n".format(
                pname, pvalue, pvalue["label"]
            )
        )
        if not hasattr(self, pname):
            if "default" in pvalue:
                if type(pvalue["default"]) is str:
                    if pvalue["default"] == "True":
                        setattr(self, pname, True)
                    if pvalue["default"] == "False":
                        setattr(self, pname, False)
                    else:
                        raise Exception(
                            "{} is boolean - default values must be True or False not {}".format(
                                pname, pvalue["default"]
                            )
                        )
                elif type(pvalue["default"]) is bool:
                    setattr(self, pname, pvalue["default"])
            else:
                setattr(self, pname, False)
        sys.stderr.write(
            "draw CB pname {} value {}\n".format(pname, getattr(self, pname))
        )
        cb = gui.checkBox(box, self, pname, pvalue["label"])
        if pname in self.QButtonHash:
            self.QButtonHash[pname].addButton(cb)
        checkAttr = pname + "Checked"
        setattr(self, checkAttr, getattr(self, pname))
        # check if inactive
        self.bgui.add(pname, cb)

    def drawSpin(self, pname, pvalue, box=None, addCheckbox=False):
        # for drawSpin - we use the origin version which already has a checkbox connected
        # TODO could change this to the same way we handle ledits with separate cbo
        # the gui spin box returns either (cb,sbox) or just sbox depending on whether there is a checkabox
        checkbox = None
        checkAttr = None
        if addCheckbox:
            if pname not in self.optionsChecked:
                self.optionsChecked[pname] = False
            checkAttr = pname + "Checked"
            setattr(self, checkAttr, self.optionsChecked[pname])
        default = 0
        if "default" in pvalue:
            default = pvalue["default"]
        if pvalue["type"] == "int":
            if not hasattr(self, pname):
                setattr(self, pname, int(default))
            elif getattr(self, pname):
                setattr(self, pname, int(getattr(self, pname)))
            else:
                setattr(self, pname, int(default))
        else:
            if not hasattr(self, pname):
                setattr(self, pname, float(default))
            elif getattr(self, pname):
                setattr(self, pname, int(getattr(self, pname)))
            else:
                setattr(self, pname, float(default))
        if addCheckbox:
            (checkbox, mySpin) = gui.spin(
                box,
                self,
                pname,
                minv=-2147483648,
                maxv=2147483647,
                label=pvalue["label"],
                checked=checkAttr,
                checkCallback=lambda: self.updateSpinCheckbox(pname),
            )
            if pname in self.QButtonHash:
                self.QButtonHash[pname].addButton(checkbox)
        else:
            mySpin = gui.spin(
                box,
                self,
                pname,
                minv=-2147483648,
                maxv=2147483647,
                label=pvalue["label"],
                checked=checkAttr,
                checkCallback=lambda: self.updateSpinCheckbox(pname),
            )

        if getattr(self, pname) is None:
            mySpin.clear()
        self.bgui.add(
            pname,
            mySpin,
            enableCallback=lambda value, clearLedit: self.enableSpin(
                value, clearLedit, checkbox, mySpin
            ),
        )

    def drawPatternQuery(self, pname, pvalue, box=None, layout=None, addCheckbox=False):
        # stores a vector consisting of root directory, pattern, search for file, search for directory, search for both
        # add mindepth maxdepth later as options...

        # At *runtime* it will use the vector to generate a set of files
        # This type of input is needed so that Bwb's iterator can be used

        # This checkbox is here when the input is optional
        # Should disable everything but itself when checked

        rootAttr = pname + "Root"
        findFileAttr = pname + "findFile"
        findDirAttr = pname + "findDir"
        patternAttr = pname + "Pattern"
        if not hasattr(self, pname) or not getattr(self, pname):
            setattr(self, pname, {})
        self.initializePatterQueryAttrs(
            pname, rootAttr, patternAttr, findFileAttr, findDirAttr
        )

        queryElements = []
        checkbox = None

        # ask for root directory

        rootLedit = gui.lineEdit(None, self, rootAttr, disabled=addCheckbox)

        # note that using lambda does not work in this function - truncates string variable so partial used instead
        browseBtn = gui.button(
            None,
            self,
            "",
            callback=partial(self.browseFileDir, attr=rootAttr, filetype="directory"),
            autoDefault=True,
            width=19,
            height=19,
            disabled=addCheckbox,
        )
        if getattr(self, rootAttr) is None:
            rootLedit.clear()
        if addCheckbox:
            checkAttr = pname + "Checked"
            if pname not in self.optionsChecked:
                self.optionsChecked[pname] = False
            setattr(self, checkAttr, self.optionsChecked[pname])
            checkbox = gui.checkBox(None, self, checkAttr, label=None)
            if pname in self.QButtonHash:
                self.QButtonHash[pname].addButton(checkbox)
            checkbox.stateChanged.connect(
                lambda: self.updateCheckbox(
                    pname, checkbox.isChecked(), getattr(self, pname)
                )
            )
            queryElements.append(checkbox)
        labelValue="Directory"
        if pvalue["label"]:
            labelValue = pvalue["label"] + " directory"
        sys.stderr.write(
            "adding filedir for pname {} using layout {}\n".format(pname, layout)
        )
        self.bwbFileEntry(
            box,
            browseBtn,
            rootLedit,
            layout=layout,
            label=labelValue + ":",
            entryType=pvalue["type"],
            checkbox=checkbox,
            placeHolderText="Enter root directory",
        )

        # query layout section
        # box for query layout

        findFileCB = gui.checkBox(None, self, findFileAttr, label="Find files")
        findDirCB = gui.checkBox(None, self, findDirAttr, label="Find directories")
        findFileCB.setDisabled(addCheckbox)
        findDirCB.setDisabled(addCheckbox)
        startCol = 0
        if addCheckbox:
            startCol += 1
            checkbox.stateChanged.connect(
                lambda: findDirCB.setEnabled(checkbox.isChecked())
            )
            checkbox.stateChanged.connect(
                lambda: findFileCB.setEnabled(checkbox.isChecked())
            )
        layout.addLinefeed(startCol=startCol)
        patternLabel = QtGui.QLabel("Pattern")
        if pvalue["label"]:
            patternLabel = QtGui.QLabel(pvalue["label"] + " pattern")
        patternLedit = gui.lineEdit(None, self, patternAttr, disabled=addCheckbox)
        patternLedit.setClearButtonEnabled(True)
        patternLedit.setPlaceholderText("Enter search pattern eg. **/*.fq")

        testBtn = gui.button(
            None, self, "Test", callback=lambda: self.testPatternQuery(pname)
        )
        testBtn.setStyleSheet(self.css)
        testBtn.setFixedSize(60, 20)

        layout.addWidget(patternLabel, width=1)
        layout.addWidget(patternLedit, width=1)
        layout.addWidget(testBtn)
        layout.addLinefeed(startCol=startCol + 1)
        layout.addWidget(findFileCB, width=1)
        layout.addLinefeed(startCol=startCol + 1)
        layout.addWidget(findDirCB, width=1)
        layout.addLinefeed(startCol=startCol + 1)

        queryElements.extend(
            [rootLedit, browseBtn, patternLedit, findFileCB, findDirCB]
        )

        self.initializePatternQueryHelpers(
            pname, rootLedit, browseBtn, patternLedit, findFileCB, findDirCB
        )
        self.bgui.addList(
            pname,
            queryElements,
            enableCallback=lambda value, clearLedit: self.enablePatternQuery(
                value,
                clearLedit,
                checkbox,
                rootLedit,
                browseBtn,
                patternLedit,
                findFileCB,
                findDirCB,
            ),
            updateCallback=lambda: self.updatePatternQuery(
                pname, rootAttr, patternAttr, findFileAttr, findDirAttr
            ),
        )
        # update the values on changes
        rootLedit.textChanged.connect(
            lambda: self.updatePatternQuery(
                pname,
                getattr(self, rootAttr),
                getattr(self, patternAttr),
                getattr(self, findFileAttr),
                getattr(self, findDirAttr),
            )
        )
        patternLedit.textChanged.connect(
            lambda: self.updatePatternQuery(
                pname,
                getattr(self, rootAttr),
                getattr(self, patternAttr),
                getattr(self, findFileAttr),
                getattr(self, findDirAttr),
            )
        )
        findFileCB.stateChanged.connect(
            lambda: self.updatePatternQuery(
                pname,
                getattr(self, rootAttr),
                getattr(self, patternAttr),
                getattr(self, findFileAttr),
                getattr(self, findDirAttr),
            )
        )
        findDirCB.stateChanged.connect(
            lambda: self.updatePatternQuery(
                pname,
                getattr(self, rootAttr),
                getattr(self, patternAttr),
                getattr(self, findFileAttr),
                getattr(self, findDirAttr),
            )
        )

    def testPatternQuery(self, pname):
        queryResults = self.generatePatternQuery(pname)
        qm = ScrollMessageBox(queryResults, None)
        qm.setWindowTitle("{} query values".format(pname))
        qm.exec_()

    def updatePatternQuery(
        self, attr, root, pattern, findFile, findDir
    ):  # updates for input - called before and after addition and deletion of input
        patternQuery = getattr(self, attr)
        patternQuery["root"] = root
        patternQuery["pattern"] = pattern
        patternQuery["findFile"] = findFile
        patternQuery["findDir"] = findDir

    def enablePatternQuery(
        self,
        value,
        clearLedit,
        checkbox,
        browseBtn,
        rootLedit,
        patternLedit,
        findFileCB,
        findDirCB,
    ):
        # first element is checkbox
        # last element is browseBtn if it exists
        if checkbox:
            checkbox.setEnabled(True)
        if not checkbox or checkbox.isChecked():
            for g in [browseBtn, rootLedit, patternLedit, findFileCB, findDirCB]:
                if g:
                    g.setEnabled(True)

    def receivePatternQuery(
        self, attr, value, rootLedit, patternLedit, findFileCB, findDirCB
    ):
        patternQuery = getattr(self, attr)
        if type(value) is str:
            patternQuery["root"] = value
            rootLedit.setText(value)
        elif type(value) is dict:
            if "root" in value:
                patternQuery["root"] = value["root"]
                rootLedit.setText(value["root"])
            if "pattern" in value:
                patternQuery["pattern"] = value["pattern"]
                patterLedit.setText(value["pattern"])
            if "findFile" in value:
                patternQuery["findFile"] = value["findFile"]
                findFileCb.setChecked(value["findFile"])
            if "findDir" in value:
                patternQuery["findDir"] = value["findDir"]
                findDirCb.setChecked(value["findDir"])

    def initializePatterQueryAttrs(
        self, attr, rootAttr, patternAttr, findFileAttr, findDirAttr
    ):
        patternQuery = getattr(self, attr)
        if patternQuery:
            setattr(self, rootAttr, patternQuery["root"])
            setattr(self, patternAttr, patternQuery["pattern"])
            setattr(self, findFileAttr, patternQuery["findFile"])
            setattr(self, findDirAttr, patternQuery["findDir"])
        else:
            setattr(self, rootAttr, "")
            setattr(self, patternAttr, "")
            setattr(self, findFileAttr, True)
            setattr(self, findDirAttr, False)

    def initializePatternQueryHelpers(
        self, pname, rootLedit, browseBtn, patternLedit, findFileCB, findDirCB
    ):
        self.helpers.add(pname, "generate", lambda: self.generatePatternQuery(pname))
        self.helpers.add(
            pname,
            "receive",
            lambda value: self.receivePatternQuery(
                pname, value, rootLedit, patternLedit, findFileCB, findDirCB
            ),
        )
        self.helpers.add(
            pname,
            "initialize",
            lambda value=None, inputType=None: self.initializePatternQuery(
                pname, rootLedit, patternLedit, findFileCB, findDirCB, value, inputType
            ),
        )
        self.helpers.add(
            pname,
            "disableGUI",
            lambda value=None, inputType=None: self.disableGUIPatternQuery(
                pname,
                rootLedit,
                browseBtn,
                patternLedit,
                findFileCB,
                findDirCB,
                value=value,
                inputType=inputType,
            ),
        )

    def initializePatternQuery(
        self, pname, rootLedit, patternLedit, findFileCB, findDirCB, value, inputType
    ):
        if value is None:
            rootLedit.clear()
        elif type(value) is str and value:
            rootLedit.setText(value)
        if type(value) is str or inputType == "str":
            return
        patternLedit.clear()
        findFileCB.setChecked(True)
        findDirCB.setChecked(False)

    def disableGUIPatternQuery(
        self,
        pname,
        rootLedit,
        browseBtn,
        patternLedit,
        findFileCB,
        findDirCB,
        value=None,
        inputType=None,
    ):
        if type(value) is dict or inputType == "dict":
            self.pname = {}
            if value:
                self.pname = value
            for element in [rootLedit, browseBtn, patternLedit, findFileCB, findDirCB]:
                element.setDisabled(True)
        elif type(value) is str or inputType == "str":
            if value is None or value is "":
                rootLedit.clear()
            else:
                rootLedit.setText(value)
            rootLedit.setDisabled(True)
            browseBtn.setDisabled(True)

    def generatePatternQuery(self, pname):
        patternQuery = getattr(self, pname)
        path = patternQuery["root"]
        pattern = patternQuery["pattern"]
        findFile = patternQuery["findFile"]
        findDir = patternQuery["findDir"]
        matches = []
        if findFile:
            matches.extend(self.getGlobFiles(path, pattern))
        if findDir:
            matches.extend(self.getGlobDirs(path, pattern))
        if matches:
            matches.sort()
        return matches

    def getGlobFiles(self, path, pattern):
        globOutput = []
        pwd = os.getcwd()
        os.chdir(path)
        for match in glob.glob(pattern, recursive=True):
            # make sure it is a file
            if os.path.isfile(match):
                globOutput.append(os.path.join(path, match))
        os.chdir(pwd)
        return globOutput

    def getGlobDir(self, path, pattern):
        globOutput = []
        pwd = os.getcwd()
        os.chdir(path)
        if pattern[-1] != "/":
            pattern = pattern + "/"
        for match in glob.glob(pattern, recursive=True):
            # make sure it is a file
            if os.path.isdir(match):
                globOutput.append(match)
        os.chdir(pwd)
        return globOutput

    def drawLedit(self, pname, pvalue, box=None, layout=None, addCheckbox=False):
        checkAttr = None
        checkbox = None
        ledit = gui.lineEdit(None, self, pname, disabled=addCheckbox)
        if addCheckbox:
            if pname not in self.optionsChecked:
                self.optionsChecked[pname] = False
            checkAttr = pname + "Checked"
            setattr(self, checkAttr, self.optionsChecked[pname])
            checkbox = gui.checkBox(None, self, checkAttr, label=None)
            if pname in self.QButtonHash:
                self.QButtonHash[pname].addButton(checkbox)
            checkbox.stateChanged.connect(
                lambda: ledit.setEnabled(checkbox.isChecked())
            )
            checkbox.stateChanged.connect(
                lambda: self.updateCheckbox(pname, checkbox.isChecked(), ledit.text())
            )
            self.bgui.add(pname, checkbox)

        self.bwbLedit(box, checkbox, ledit, layout=layout, label=pvalue["label"])
        # check if the value is none - then we clear it
        if getattr(self, pname) is None or getattr(self, pname) is "":
            ledit.clear()
        self.bgui.add(
            pname,
            ledit,
            enableCallback=lambda value, clearLedit: self.enableLedit(
                value, clearLedit, checkbox, ledit
            ),
        )
        if addCheckbox:
            self.updateCheckbox(pname, checkbox.isChecked(), ledit.text())
        return {"ledit": ledit, "checkbox": checkbox}

    def enableLedit(self, value, clearLedit, cb, ledit):
        sys.stderr.write("cb is {} ledit is {}\n".format(cb, ledit))
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

    def drawFileDirElements(
        self, pname, pvalue, box=None, layout=None, addCheckbox=False
    ):
        checkbox = None
        ledit = gui.lineEdit(None, self, pname, disabled=addCheckbox)
        # note that using lambda does not work in this function - truncates string variable so partial used instead
        button = gui.button(
            None,
            self,
            "",
            callback=partial(self.browseFileDir, attr=pname, filetype=pvalue["type"]),
            autoDefault=True,
            width=19,
            height=19,
            disabled=addCheckbox,
        )
        if getattr(self, pname) is None:
            ledit.clear()
        if addCheckbox:
            checkAttr = pname + "Checked"
            if pname not in self.optionsChecked:
                self.optionsChecked[pname] = False
            setattr(self, checkAttr, self.optionsChecked[pname])
            checkbox = gui.checkBox(None, self, checkAttr, label=None)
            if pname in self.QButtonHash:
                self.QButtonHash[pname].addButton(checkbox)
            sys.stderr.write("updating filedir {}\n".format(pname))
            checkbox.stateChanged.connect(
                lambda: self.updateCheckbox(
                    pname, checkbox.isChecked(), getattr(self, pname)
                )
            )
            self.bgui.add(pname, checkbox)
        labelValue = pvalue["label"]
        if labelValue is None:
            labelValue = ""
        sys.stderr.write(
            "adding filedir for pname {} using layout {}\n".format(pname, layout)
        )
        self.bwbFileEntry(
            box,
            button,
            ledit,
            layout=layout,
            label=labelValue + ":",
            entryType=pvalue["type"],
            checkbox=checkbox,
        )
        self.bgui.add(pname, ledit)
        self.bgui.add(
            pname,
            button,
            enableCallback=lambda value, clearLedit: self.enableFileDir(
                value, clearLedit, checkbox, ledit, button
            ),
        )
        if addCheckbox:
            self.updateCheckbox(pname, checkbox.isChecked(), ledit.text())

    def enableFileDir(self, value, clearLedit, cb, ledit, btn):
        sys.stderr.write("cb is {} ledit is {}\n".format(cb, ledit))
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

    def drawTextBoxBtnRules(self, boxEdit, ledit, addBtn, removeBtn):
        if not ledit.text():
            ledit.clear()
            addBtn.setEnabled(False)
        if not boxEdit.selectedItems():
            removeBtn.setEnabled(False)

    def enableTextBox(
        self, value, clearLedit, checkbox, browseBtn, boxEdit, ledit, addBtn, removeBtn
    ):
        # first element is checkbox
        # last element is browseBtn if it exists
        if checkbox:
            checkbox.setEnabled(True)

        if not checkbox or checkbox.isChecked():
            ledit.clear()
            boxEdit.setEnabled(True)
            ledit.setEnabled(True)
            for g in [browseBtn, addBtn, removeBtn]:
                if g:
                    g.setEnabled(True)
        self.drawTextBoxBtnRules(boxEdit, ledit, addBtn, removeBtn)

    def updateTextBox(
        self, attr, ledit, boxEdit
    ):  # updates for input - called before and after addition and deletion of input
        if hasattr(self, attr):
            value = getattr(self, attr)
            ledit.clear()
            if (
                value is None
            ):  # explicitly check for this to avoid text None from appearing
                boxEdit.clear()
            else:
                boxEdit.clear()
                boxEdit.addItems(value)
                self.updateBoxEditValue(attr, boxEdit)

    def updateBoxEditValue(self, attr, boxEdit):
        myItems = []
        sys.stderr.write("{} objects in {}\n".format(boxEdit.count(), boxEdit))
        for i in range(boxEdit.count()):
            myItems.append(boxEdit.item(i).text())
            sys.stderr.write(
                "add item number {} value {}\n".format(i, boxEdit.item(i).text())
            )
        setattr(self, attr, myItems)

    def drawTextBox(self, pname, pvalue, box=None, layout=None, addCheckbox=False):
        # for multiple files or directories draw a widget list with a line edit below and buttons to browse, add, delete, submit and reload
        # is used for multiple entry of any type

        sys.stderr.write(
            "adding text box for pname{} pvalue {} box {} layout {} addCheckbox {}\n".format(
                pname, pvalue, box, layout, addCheckbox
            )
        )
        disabledFlag = False
        checkbox = None
        if addCheckbox or self.inputConnections.isConnected(pname):
            disabledFlag = True
        elements = []
        # add checkbox if necessary
        if addCheckbox:
            value = None
            if pname not in self.optionsChecked:
                self.optionsChecked[pname] = False
            if hasattr(self, pname):
                value = getattr(self, pname)
            checkAttr = (
                pname + "Checked"
            )  # this is not actually used but needed for orange gui checkbox element
            setattr(self, checkAttr, self.optionsChecked[pname])
            checkbox = gui.checkBox(None, self, checkAttr, label=None)
            if pname in self.QButtonHash:
                self.QButtonHash[pname].addButton(checkbox)
            checkbox.stateChanged.connect(
                lambda: self.updateCheckbox(pname, checkbox.isChecked(), value=value)
            )
            elements.append(checkbox)
            if checkbox.isChecked() and not self.inputConnections.isConnected(pname):
                disabledFlag = False

        # line edit for manual file entry
        leditAttr = pname + "Ledit"
        if not hasattr(self, leditAttr):
            setattr(self, leditAttr, None)
        ledit = gui.lineEdit(None, self, leditAttr, disabled=disabledFlag)
        ledit.setClearButtonEnabled(True)
        if pvalue["type"] == "file list":
            ledit.setPlaceholderText("Enter file")
        elif pvalue["type"] == "directory list":
            ledit.setPlaceholderText("Enter directory")
        else:
            ledit.setPlaceholderText("Enter parameter")

        ledit.setStyleSheet(":disabled { color: #282828}")
        ledit.clear()
        elements.append(ledit)

        # add boxEdit layout
        layoutAttr = pname + "Layout"
        setattr(
            self,
            layoutAttr,
            self.addBoxEdit(
                pname,
                pvalue,
                layout,
                ledit,
                checkbox,
                elements=elements,
                disabledFlag=disabledFlag,
            ),
        )

    def addBoxEdit(
        self, pname, pvalue, layout, ledit, checkbox, elements=None, disabledFlag=False
    ):
        # setup boxEdit - boxEdit element values other than the list itself are not tracked and not saved in settings
        if not hasattr(self, pname):
            setattr(self, pname, None)
        boxEdit = DragAndDropList(self)
        boxEdit.setDragDropMode(QtWidgets.QAbstractItemView.InternalMove)
        boxEdit.setSelectionMode(QtWidgets.QAbstractItemView.ExtendedSelection)
        boxEdit.setStyleSheet(":disabled { color: #282828}")
        boxEdit.setMinimumHeight(60)

        # fill boxEdit - ONLY part that is tracked
        if hasattr(self, pname):
            entryList = getattr(self, pname)
            sys.stderr.write(
                "filling with {} of type {}\n".format(entryList, type(entryList))
            )
            if type(entryList) == list:
                boxEdit.addItems(entryList)
            else:
                boxEdit.addItems([entryList])
        boxEdit.setDisabled(disabledFlag)
        if elements:
            elements.append(boxEdit)
        # buttons
        browseBtn = None
        addBtn = gui.button(
            None,
            self,
            "",
            callback=lambda: self.addLedit(pname, ledit, boxEdit, addBtn),
            autoDefault=False,
            disabled=disabledFlag,
        )
        removeBtn = gui.button(
            None,
            self,
            "",
            callback=lambda: self.removeItem(pname, boxEdit, removeBtn),
            autoDefault=False,
            disabled=disabledFlag,
        )
        if pvalue["type"] == "file list":
            browseBtn = gui.button(
                None,
                self,
                "",
                callback=partial(self.browseFiles, boxEdit=boxEdit, attr=pname),
                autoDefault=False,
                disabled=disabledFlag,
            )
        elif pvalue["type"] == "directory list":
            browseBtn = gui.button(
                None,
                self,
                "",
                callback=partial(self.browseDirs, boxEdit=boxEdit, attr=pname),
                autoDefault=False,
                disabled=disabledFlag,
            )
        # set styles
        if browseBtn:
            browseBtn.setStyleSheet(self.browseCSS)
            elements.append(browseBtn)
        for btn in (addBtn, removeBtn):
            if btn:
                btn.setStyleSheet(self.addRemoveCSS)
                elements.append(btn)
        # set icons
        self.bgui.addList(
            pname,
            elements,
            enableCallback=lambda value, clearLedit: self.enableTextBox(
                value,
                clearLedit,
                checkbox,
                browseBtn,
                boxEdit,
                ledit,
                addBtn,
                removeBtn,
            ),
            updateCallback=lambda: self.updateTextBox(pname, ledit, boxEdit),
        )
        if browseBtn:
            browseBtn.setIcon(self.browseIcon)
        addBtn.setIcon(self.addIcon)
        removeBtn.setIcon(self.removeIcon)

        # check rules for buttons
        self.drawTextBoxBtnRules(boxEdit, ledit, addBtn, removeBtn)

        # connects from ledit and boxEdit to buttons
        ledit.textChanged.connect(lambda: addBtn.setEnabled(bool(ledit.text())))
        boxEdit.itemSelectionChanged.connect(
            lambda: removeBtn.setEnabled(bool(boxEdit.selectedItems()))
        )

        # connect to changes in drag and drop
        boxEdit.itemMoved.connect(lambda: self.updateBoxEditValue(pname, boxEdit))

        # layout section
        filesBoxLeditLayout = QtGui.QVBoxLayout()
        # add to the main parameters box
        myBox = gui.vBox(None)
        myBox.layout().addLayout(filesBoxLeditLayout)
        startCol = 0
        if "label" in pvalue and pvalue["label"]:
            label = QtGui.QLabel(pvalue["label"] + ":")
        else:
            label = QtGui.QLabel(" ")
        label.setAlignment(Qt.AlignTop)
        if checkbox:
            layout.addWidget(checkbox)

        layout.addWidget(label)
        layout.addWidget(myBox, width=2, linefeed=2)
        # line layout
        lineLayout = BwbGridLayout()

        lineLayout.addWidget(ledit)
        if browseBtn:
            lineLayout.addWidget(browseBtn)
        lineLayout.addWidget(addBtn)
        lineLayout.addWidget(removeBtn)

        filesBoxLeditLayout.addWidget(boxEdit)
        filesBoxLeditLayout.addLayout(lineLayout.layout())

        # just in case the state has changed while drawing this - do an update
        self.updateBoxEditValue(pname, boxEdit)
        return filesBoxLeditLayout

    def addLedit(self, attr, ledit, boxEdit, addBtn):
        # adds text in ledit to items in boxEdit
        if ledit.text():
            boxEdit.addItem(ledit.text())
            sys.stderr.write("adding {} to {}\n".format(ledit.text(), attr))
            ledit.clear()
            addBtn.setEnabled(False)
            self.updateBoxEditValue(attr, boxEdit)
            sys.stderr.write("boxEdit values is {}\n".format(getattr(self, attr)))

    def removeItem(self, attr, boxEdit, removeBtn):
        if boxEdit.selectedItems():
            for item in boxEdit.selectedItems():
                boxEdit.takeItem(boxEdit.row(item))
            self.updateBoxEditValue(attr, boxEdit)
        if not boxEdit.count():
            removeBtn.setEnabled(False)

    def findIterables(self):
        retList = []
        if self.data["parameters"] is None:
            return retList
        for pname, pvalue in self.data["parameters"].items():
            if "list" in pvalue["type"] or pvalue["type"] == "patternQuery":
                retList.append(pname)
        return retList

    def drawExec(self, box=None):
        if not hasattr(self, "useDockerfile"):
            setattr(self, "useDockerfile", False)
        # find out if there are triggers
        self.candidateTriggers = []
        if self.data["inputs"] is not None:
            for pname in self.data["inputs"]:
                self.candidateTriggers.append(pname)
        # make sure that runTriggers are a subset of candidateTriggers
        self.runTriggers=[trigger for trigger in self.runTriggers if trigger in self.data["inputs"]]
        #initilialize any runTriggers
        for attr in self.runTriggers:
            self.triggerReady[attr] = False
        
        # initialize the exec state
        self.execLayout = QtGui.QGridLayout()
        self.execLayout.setSpacing(5)

        self.graphicsMode = QtGui.QCheckBox("Export graphics", self)
        self.graphicsMode.setChecked(self.exportGraphics)
        self.graphicsMode.stateChanged.connect(self.graphicsModeChange)

        self.testMode = QtGui.QCheckBox("Test mode", self)
        self.testMode.setChecked(self.useTestMode)
        self.testMode.stateChanged.connect(self.testModeChange)

        # self.dockerMode=QtGui.QCheckBox('Build container',self)
        # self.dockerMode.setChecked(self.useDockerfile)
        # self.dockerMode.stateChanged.connect(self.dockerModeChange)

        self.cboRunMode = QtGui.QComboBox()
        self.cboRunMode.addItem("Manual")
        self.cboRunMode.addItem("Automatic")
        if self.candidateTriggers:
            self.cboRunMode.addItem("Triggered")
        elif self.runMode == 2:  # reset in case there is a change in an older workflow
            self.runMode = 0
        self.cboRunMode.setCurrentIndex(self.runMode)
        if self.candidateTriggers:
            self.execBtn = QtGui.QToolButton(self)
            self.execBtn.setText("Select Triggers")
            self.execMenu = QtGui.QMenu(self)
            self.triggerMenu = {}
            for attr in self.candidateTriggers:
                action = self.execMenu.addAction(attr)
                action.setCheckable(True)
                action.setChecked(bool(attr in self.runTriggers))
                action.changed.connect(self.chooseTrigger)
                self.triggerMenu[action] = attr
            self.execBtn.setMenu(self.execMenu)
            self.execBtn.setPopupMode(QtGui.QToolButton.InstantPopup)
            if self.runMode == 2:
                self.execBtn.setEnabled(True)
            else:
                self.execBtn.setEnabled(False)
        self.cboRunMode.currentIndexChanged.connect(self.runModeChange)
        myLabel = QtGui.QLabel("RunMode:")
        myLabel.setAlignment(Qt.AlignRight | Qt.AlignVCenter)
        myLabel.setAlignment(Qt.AlignRight | Qt.AlignVCenter)
        self.btnRun = gui.button(
            None, self, "Start", callback=lambda: self.onRunClicked(button=self.btnRun)
        )
        self.btnRun.setStyleSheet(self.css)
        self.btnRun.setFixedSize(60, 20)
        self.btnStop = gui.button(None, self, "Stop", callback=self.onStopClicked)
        self.btnStop.setStyleSheet(self.css)
        self.btnStop.setFixedSize(60, 20)
        self.btnStop.setEnabled(False)
        self.execLayout.addWidget(self.btnRun, 1, 0)
        self.execLayout.addWidget(self.btnStop, 1, 1)
        # self.execLayout.addWidget(self.iterateBtn,1,2)
        self.execLayout.addWidget(self.graphicsMode, 1, 2)
        self.execLayout.addWidget(self.testMode, 1, 3)
        # self.execLayout.addWidget(self.dockerMode,1,3)
        self.execLayout.addWidget(myLabel, 1, 4)
        self.execLayout.addWidget(self.cboRunMode, 1, 5)
        if self.candidateTriggers:
            self.execLayout.addWidget(self.execBtn, 1, 6)
        box.layout().addLayout(self.execLayout)

    def testModeChange(self):
        self.useTestMode = self.testMode.isChecked()

    def dockerModeChange(self):
        self.useDockerfile = self.dockerMode.isChecked()

    def graphicsModeChange(self):
        self.exportGraphics = self.graphicsMode.isChecked()

    def runModeChange(self):
        self.runMode = self.cboRunMode.currentIndex()
        if self.candidateTriggers:
            if self.runMode == 2:
                self.execBtn.setEnabled(True)
            else:
                self.execBtn.setEnabled(False)
        self.checkTrigger()

    def chooseTrigger(self):
        action = self.execMenu.sender()
        attr = self.triggerMenu[action]
        checked = action.isChecked()
        if attr is None or checked is None:
            return
        if checked and attr not in self.runTriggers:
            self.runTriggers.append(attr)
            self.triggerReady[attr] = False
        elif not checked and attr in self.runTriggers:
            (self.runTriggers).remove(attr)
            self.triggerReady[attr] = False

    def chooseIterable(self):
        action = self.iterablesMenu.sender()
        attr = self.iterablesMenuItems[action]
        checked = action.isChecked()
        if attr is None or checked is None:
            self.iterateSettings["iteratedAttrs"].append(attr)
        if checked and attr not in self.iterateSettings["iteratedAttrs"]:
            self.iterateSettings["iteratedAttrs"].append(attr)
        if not checked and attr in self.iterateSettings["iteratedAttrs"]:
            del self.iterateSettings["iteratedAttrs"][attr]

    def chooseIP(self):
        action = self.IPMenu.sender()
        attr = self.IPMenuItems[action]
        checked = action.isChecked()
        if attr is None or checked is None:
            return

    def chooseScheduler(self):
        action = self.schedulerMenu.sender()
        attr = self.schedulerMenuItems[action]
        checked = action.isChecked()
        if attr is None or checked is None:
            return

    def checkTrigger(self, inputReceived=False):
        # this should be checked any time there is a change
        # but only triggers when an input is received
        if self.runMode == 0:  # manual - only go on start button
            return
        elif self.runMode == 1:  # automatic same as pushing start button
            self.onRunClicked()
            return
        elif self.runTriggers:
            if not inputReceived:
                return
            # check if the input triggers are set
            sys.stderr.write(
                "Checking triggers with runTriggers{}\n".format(self.runTriggers)
            )
            for trigger in self.runTriggers:
                if not self.inputConnections.isSet(trigger):
                    return
                if not self.triggerReady[trigger]:
                    return
            self.onRunClicked()

    def bwbFileEntry(
        self,
        widget,
        button,
        ledit,
        icon=browseIcon,
        layout=None,
        label=None,
        entryType="file",
        checkbox=None,
        placeHolderText=None,
    ):
        button.setIcon(icon)
        button.setStyleSheet(self.browseCSS)
        ledit.setClearButtonEnabled(True)
        if placeHolderText is None:
            ledit.setPlaceholderText("Enter {}".format(entryType))
        else:
            ledit.setPlaceholderText(placeHolderText)
        if checkbox:
            layout.addWidget(checkbox)
        if label:
            # myLabel=QtGui.QLabel(label)
            # myLabel.setAlignment(Qt.AlignLeft | Qt.AlignVCenter)
            layout.addWidget(QtGui.QLabel(label))
        layout.addWidget(ledit)
        layout.addWidget(button)
        if not hasattr(layout, "added") or not getattr(layout, "added"):
            widget.layout().addLayout(layout.layout())
            setattr(layout, "added", True)
        layout.addLinefeed()

    def bwbLedit(self, widget, checkbox, ledit, layout=None, label=None):
        ledit.setClearButtonEnabled(True)
        ledit.setPlaceholderText("Enter parameter")
        if checkbox:
            layout.addWidget(checkbox)
            layout.addWidget(QtGui.QLabel(label))
        else:
            layout.addWidget(QtGui.QLabel(label), space=1)
        layout.addWidget(ledit)
        if not hasattr(layout, "added") or not getattr(layout, "added"):
            widget.layout().addLayout(layout.layout())
            setattr(layout, "added", True)
        layout.addLinefeed()

    def browseFiles(self, boxEdit, attr=None):
        defaultDir = "/root"
        if os.path.exists("/data"):
            defaultDir = "/data"
        boxEdit.addItems(
            QtWidgets.QFileDialog.getOpenFileNames(
                self,
                caption="Choose files",
                directory=defaultDir,
                filter="Any file (*.*)",
            )[0]
        )
        if attr:
            self.updateBoxEditValue(attr, boxEdit)

    def browseDirs(self, boxEdit, attr=None):
        dlg = getExistingDirectories()
        if dlg.exec_() == QtWidgets.QDialog.Accepted:
            boxEdit.addItems(dlg.selectedFiles())
            if attr:
                self.updateBoxEditValue(attr, boxEdit)

    def browseFileDir(self, attr, filetype=None):
        defaultDir = "/root"
        if os.path.exists("/data"):
            defaultDir = "/data"
        if filetype == "file":
            myFile = QtWidgets.QFileDialog.getOpenFileName(
                self, "Locate file", defaultDir
            )[0]
            if myFile:
                setattr(self, attr, myFile)
                dirAttr = attr + "Dir"
                setattr(self, dirAttr, os.path.dirname(myFile))
        else:
            myDir = QtWidgets.QFileDialog.getExistingDirectory(
                self, caption="Locate directory", directory=defaultDir
            )
            if myDir:
                setattr(self, attr, myDir)

    def updateCheckbox(self, pname, state, value=None):
        self.optionsChecked[pname] = state
        sys.stderr.write(
            "updating checkbox pname {} connect {} isChecked {}\n".format(
                pname, self.inputConnections.isConnected(pname), state
            )
        )
        if not (self.inputConnections.isConnected(pname)) and state:
            self.bgui.enable(pname, value)
            sys.stderr.write("enabled\n")
        elif not (self.inputConnections.isConnected(pname)) and not state:
            self.bgui.disable(pname, ignoreCheckbox=True)
            sys.stderr.write("disabled\n")

    def enableSpin(self, value, clearLedit, cb, spin):
        if cb is None:
            spin.setEnabled(True)
        else:
            if cb.isChecked():
                spin.setEnabled(True)
            else:
                spin.setEnabled(False)
            cb.setEnabled(True)

    def updateSpinCheckbox(self, pname):
        checkAttr = pname + "Checked"
        self.optionsChecked[pname] = getattr(self, checkAttr)

    def removeRequiredFromOptionsChecked(self):
        for pname in self.data["requiredParameters"]:
            if pname in self.optionsChecked:
                del(self.optionsChecked[pname])
                
    def convertOrangeTypes(self, inputType):
        inputTypeArray = inputType.split(".")
        return inputTypeArray[-1]

    # Handle inputs
    def handleInputs(self, attr, value, sourceId, test=False):
        inputType = None
        sys.stderr.write(
            "INPUT RECEIVED value is {} sourceId is {}\n".format(value, sourceId)
        )
        # check for test mode pre-signal
        if type(value) is str and value[0:8] == "__purge ":
            # remove signal - this used to be None which was also passed when the value actually was None
            self.inputConnections.remove(attr, sourceId)
            sys.stderr.write(
                "sig handler removing {} disabled {}\n".format(
                    attr, self.inputConnections.isConnected(attr)
                )
            )
            inputType = self.convertOrangeTypes(value.split()[-1])
            self.initializeAttr(attr, inputType=inputType)
            if attr in self.runTriggers:
                self.triggerReady[attr] = False
            value = None
        elif type(value) is str and value[0:6] == "__add ":
            # this is a new connection - and not just a signal - does not activate
            inputType = self.convertOrangeTypes(value.split()[1])
            self.initializeAttr(attr, inputType=inputType)
            disableFtn = self.helpers.function(attr, "disableGUI")
            sys.stderr.write(
                "attr is {} disableFtn is {} inputType is {}\n".format(
                    attr, disableFtn, inputType
                )
            )
            if disableFtn:
                self.bgui.disable(
                    attr, disableFtn=lambda: disableFtn(inputType=inputType)
                )
            else:
                self.bgui.clear(attr, activate=True)
                self.bgui.disable(attr)
            self.inputConnections.add(attr, sourceId)
            sys.stderr.write(
                "sig handler adding node with no signal: attr {} sourceId {} value {}\n".format(
                    attr, sourceId, value
                )
            )
            sys.stderr.write(
                "connection dict value for {} is {}\n".format(
                    attr, self.inputConnections._dict[attr]
                )
            )
            self.checkTrigger(inputReceived=False)
            # don't go through the default update since we have already disabled the new connection
            return
        else:
            if test:
                self.saveBashFile = None
            self.testMode.setChecked(test)
            self.inputConnections.add(attr, sourceId)
            sys.stderr.write(
                "sig handler adding input: attr {} sourceId {} value {}\n".format(
                    attr, sourceId, value
                )
            )
            sys.stderr.write(
                "connection dict value for {} is {}\n".format(
                    attr, self.inputConnections._dict[attr]
                )
            )
            # put transfer value here
            receiveFtn = self.helpers.function(attr, "receive")
            if receiveFtn:
                receiveFtn(value)
            else:
                setattr(self, attr, value)
            if attr in self.runTriggers:
                sys.stderr.write(
                    "trigger value for attr {} is {}\n".format(
                        attr, self.triggerReady[attr]
                    )
                )
                self.triggerReady[attr] = True
            self.checkTrigger(inputReceived=True)
        if value is not None:
            inputType = type(value)
        self.updateGui(attr, value, inputType=inputType)

    def initializeAttr(self, attr, inputType=None):
        if not hasattr(self, attr):
            setattr(self, attr, None)
        value = getattr(self, attr)
        initFunction = self.helpers.function(attr, "initialize")
        if initFunction:
            initFunction(value=value, inputType=inputType)
        else:
            if type(getattr(self, attr)) is str:
                setattr(self, attr, "")
            else:
                setattr(self, attr, None)

    def updateGui(self, attr, value, removeFlag=False, disableFtn=None, inputType=None):
        # disables manual input when the value has been given by an input connection
        if self.inputConnections.isConnected(attr):
            if disableFtn:
                if value is not None:
                    inputType = type(value)
                self.bgui.disable(attr, disableFtn=lambda: disableFtn(value, inputType))
            else:
                self.bgui.disable(attr)
        else:
            sys.stderr.write("enabling {} with {}\n".format(attr, value))
            self.bgui.enable(attr, value)

    # Generate commands and run job

    def startJob(self):
        if self.jobRunning:
            return
        self.hostVolumes = {}
        # check for missing parameters and volumes
        missingParms = self.checkRequiredParms()

        if missingParms:
            self.console.append("missing required parameters: {}".format(missingParms))
            return
        missingVols = self.getRequiredVols()
        if missingVols:
            self.console.append(
                "missing or incorrect volume mappings to: {}".format(missingVols)
            )
            return

        # get ready to start
        #make sure that requiredParms are not in options checked
        self.removeRequiredFromOptionsChecked()
        attrList = self.__dict__.keys()
        self.bgui.disableAll()
        self.disableExec()
        self.jobRunning = True

        # generate any needed data

        cmd = self.generateCmdFromData()
        self.envVars = {}
        self.getEnvironmentVariables()
        if hasattr(self, "nWorkers") and self.nWorkers is not None:
            if self.nWorkers:
                self.iterateSettings["nWorkers"] = self.nWorkers
            else:
                self.iterateSettings["nWorkers"] = multiprocessing.cpu_count()
        #        try:
        imageName = "{}:{}".format(self._dockerImageName, self._dockerImageTag)
        self.pConsole.writeMessage(
            "Generating Docker command from image {}\nVolumes {}\nCommands {}\nEnvironment {}\n".format(
                imageName, self.hostVolumes, cmd, self.envVars
            )
        )
        if self.iterate:
            breakpoint(message="cmd is {}".format(cmd))
            cmds = self.findIteratedFlags(cmd)
            breakpoint(message="cmds are {}".format(cmds))
        else:
            cmds = [cmd]
        # generate cmds here
        self.status = "running"
        self.setStatusMessage("Running...")
        sys.stderr.write("cmds are {}\n".format(cmds))
        if hasattr(self, "useScheduler") and self.useScheduler:
            self.dockerClient.create_container_external(
                imageName,
                hostVolumes=self.hostVolumes,
                cmds=cmds,
                environment=self.envVars,
                consoleProc=self.pConsole,
                exportGraphics=self.exportGraphics,
                portMappings=self.portMappings(),
                testMode=self.useTestMode,
                logFile=self.saveBashFile,
                scheduleSettings=None,
                iterateSettings=self.iterateSettings,
                iterate=self.iterate,
            )
        else:
            self.dockerClient.create_container_iter(
                imageName,
                hostVolumes=self.hostVolumes,
                cmds=cmds,
                environment=self.envVars,
                consoleProc=self.pConsole,
                exportGraphics=self.exportGraphics,
                portMappings=self.portMappings(),
                testMode=self.useTestMode,
                logFile=self.saveBashFile,
                scheduleSettings=None,
                iterateSettings=self.iterateSettings,
                iterate=self.iterate,
            )
        # except BaseException as e:
        # self.bgui.reenableAll(self)
        # self.reenableExec()
        # self.jobRunning=False
        # self.pConsole.writeMessage("unable to start Docker command "+ str(e))
        # self.setStatusMessage('Error')

    def portMappings(self):
        if "portMappings" in self.data:
            mappings = []
            for mapping in self.data["portMappings"]:
                sys.stderr.write("mapping is {}\n".format(mapping))
                if mapping["attr"] and hasattr(self, mapping["attr"]):
                    hostPort = getattr(self, mapping["attr"])
                    sys.stderr.write("{} hostPort value\n".format(hostPort))
                    mappings.append("{}:{}".format(hostPort, mapping["containerPort"]))
            if mappings:
                return mappings
        return None

    def checkRequiredParms(self):
        for parm in self.data["requiredParameters"]:
            if hasattr(self, parm):
                if not getattr(self, parm) and getattr(self, parm) != 0:
                    return parm
        return None

    def getRequiredVols(self):
        # get all the autoMaps
        # the mountpoint is passed because it will be converted later into the global hostpath
        bwbDict = {}
        if "autoMap" in self.data and self.data["autoMap"]:
            for bwbVolume, containerVolume in self.dockerClient.bwbMounts.items():
                self.hostVolumes[containerVolume] = bwbVolume
                bwbDict[containerVolume] = bwbVolume
        if "volumeMappings" in self.data and self.data["volumeMappings"]:
            for mapping in self.data["volumeMappings"]:
                conVol = mapping["conVolume"]
                attr = mapping["attr"]
                if not hasattr(self, attr) or not getattr(self, attr):
                    if (
                        "autoMap" in self.data
                        and self.data["autoMap"]
                        and conVol in bwbDict
                    ):
                        setattr(self, attr, bwbDict[conVol])
                    else:
                        return conVol
                self.hostVolumes[conVol] = getattr(self, attr)
        return None

    def joinFlagValue(self, flag, value):
        if flag:
            if flag.strip():
                if flag.strip()[-1] == "=":
                    return flag.strip() + str(value).strip()
                else:
                    return flag.strip() + " " + str(value).strip()
            else:
                return flag + str(value).strip()
        return value

    def flagString(self, pname):
        if "parameters" in self.data and pname in self.data["parameters"]:
            pvalue = self.data["parameters"][pname]
            sys.stderr.write("checkign flag with pname {} pvalue {}\n".format(pname, pvalue))
            if "flag" in pvalue and hasattr(self, pname):
                flagName = pvalue["flag"]
                if flagName is None:
                    flagName = ""
                flagValue = self.getAttrValue(pname)
                if pvalue["type"] == "patternQuery":
                    if flagValue:
                        return " ".join([flagName] + flagValue)
                    else:
                        return flagName
                elif pvalue["type"] == "bool":
                    if flagValue:
                        return flagName
                    return None
                elif pvalue["type"] == "file":
                    sys.stderr.write("flagValue is {}\n".format(flagValue))
                    if not flagValue:
                        return None
                    filename = str(flagValue)
                    if filename.strip():
                        hostFilename = self.bwbPathToContainerPath(
                            filename, isFile=True, returnNone=False
                        )
                        sys.stderr.write(
                            "orig file {} convert to container {}\n".format(
                                filename, hostFilename
                            )
                        )
                        return self.joinFlagValue(flagName, hostFilename)
                    return None
                elif pvalue["type"] == "file list":
                    # check whether it is iterated
                    files = flagValue
                    sys.stderr.write(
                        "flagnName {} files are {}\n".format(flagName, files)
                    )
                    if files:
                        hostFiles = []
                        for f in files:
                            hostFiles.append(
                                self.bwbPathToContainerPath(
                                    f, isFile=True, returnNone=False
                                )
                            )
                            sys.stderr.write(
                                "flagnName {} adding file {} hostFiles {}\n".format(
                                    flagName, f, hostFiles
                                )
                            )
                        if flagName:
                            return " ".join([flagName] + hostFiles)
                        else:
                            return " ".join(hostFiles)
                    return None
                elif pvalue["type"] == "directory":
                    path = str(flagValue)
                    if path:
                        hostPath = self.bwbPathToContainerPath(path, returnNone=False)
                    else:
                        raise Exception(
                            "Exception -no path flagValue is {} path is {}\n".format(
                                flagValue, path
                            )
                        )
                    return self.joinFlagValue(flagName, str(hostPath))
                elif pvalue["type"][-4:] == "list":
                    if flagName:
                        return flagName + " " + " ".join(flagValue)
                    else:
                        return " ".join(flagValue)
                elif flagValue is not None:
                    return self.joinFlagValue(flagName, flagValue)
        return None

    def getAttrValue(self, attr):
        # wrapper - returns the value of self.attr if there is not generate function
        if hasattr(self, attr):
            value = getattr(self, attr)
            generateFunction = self.helpers.function(attr, "generate")
            if generateFunction:
                return generateFunction()
            return value
        return None

    def generateCmdFromData(self):
        flags = []
        args = []
        # map port variables if necessary
        if not hasattr(self, "portVars"):
            sys.stderr.write("Initializing portVars\n")
            self.portVars = []
            if "portMappings" in self.data:
                for mapping in self.data["portMappings"]:
                    sys.stderr.write("adding {} to portVars\n".format(mapping["attr"]))
                    self.portVars.append(mapping["attr"])
        if self.data["parameters"] is None:
            sys.stderr.write("No parameters skipping parameter parsing\n")
            return self.generateCmdFromBash(
                self.data["command"], flags=flags, args=args
            )
        for pname, pvalue in self.data["parameters"].items():
            if pname in self.portVars:
                continue
            sys.stderr.write(
                "Parsing parameters: pname {} pvalue{}\n".format(pname, pvalue)
            )
            # possible to have an requirement or parameter that is not in the executable line
            # environment variables can have a Null value in the flags field
            # arguments are the only type that have no flag
            if "argument" in pvalue:
                if (
                    self.iterate
                    and hasattr(self, "iterateSettings")
                    and "iteratedAttrs" in self.iterateSettings
                    and pname in self.iterateSettings["iteratedAttrs"]
                ):
                    fStr = "_iterate{{{}}}".format(pname)
                else:
                    fStr = self.flagString(pname)
                if fStr and fStr is not None:
                    args.append(fStr)
                continue
            # do not add to arguments if it is an environment variable and there is no flag
            # if you really want it added put a space in the flag field
            if pvalue["flag"] is None or pvalue["flag"] == "" and "env" in pvalue:
                sys.stderr.write("Found env: pname {} pvalue{} which is added to arguments \n".format(pname, pvalue))
                continue
            # if required or checked then it is added to the flags
            addParms = False

            # checkattr is needed for the orange gui checkboxes but is not otherwise updated
            if pname in self.optionsChecked and self.optionsChecked[pname]:
                addParms = True

            # also need to check for booleans which are not tracked by optionsChecked
            if (
                pvalue["type"] == "bool"
                and hasattr(self, pname)
                and getattr(self, pname)
            ):
                addParms = True

            if self.iterate:
                if (
                    hasattr(self, "iterateSettings")
                    and "iteratedAttrs" in self.iterateSettings
                    and pname in self.data["requiredParameters"]
                    and hasattr(self, pname)
                ):
                    addParms = True

            elif pname in self.data["requiredParameters"] and hasattr(self, pname):
                addParms = True

            if addParms:
                if (
                    self.iterate
                    and hasattr(self, "iterateSettings")
                    and "iteratedAttrs" in self.iterateSettings
                    and pname in self.iterateSettings["iteratedAttrs"]
                ):
                    fStr = "_iterate{{{}}}".format(pname)
                else:
                    fStr = self.flagString(pname)
                sys.stderr.write("fStr is {}\n".format(fStr))
                sys.stderr.write("pvalue flag is {}\n".format(pvalue["flag"]))
                if fStr:
                    flags.append(fStr)

        return self.generateCmdFromBash(self.data["command"], flags=flags, args=args)

    def iteratedfString(self, pname):
        # flagstrings and findIteratedflags put the iterated names in the correct order place
        # this routine creates a list of values either flagstrings or other values to substitute into the command strings

        if "parameters" in self.data and pname in self.data["parameters"]:
            pvalue = self.data["parameters"][pname]
            # check if there is a flag
            if "flag" in pvalue and hasattr(self, pname):
                flagName = pvalue["flag"]
            else:
                flagName = ""

            # get flag values and groupSize
            groupSize = 1
            if (
                "data" in self.iterateSettings
                and pname in self.iterateSettings["data"]
                and "groupSize" in self.iterateSettings["data"][pname]
                and self.iterateSettings["data"][pname]["groupSize"]
            ):
                groupSize = int(self.iterateSettings["data"][pname]["groupSize"])

            flagValues = self.getAttrValue(pname)
            sys.stderr.write("original flagvalues are {}\n".format(flagValues))
            # make list of tuplets of groupSize
            flagValues = list(
                zip_longest(*[iter(flagValues)] * groupSize, fillvalue=flagValues[-1])
            )
            sys.stderr.write("chunked flagvalues are {}\n".format(flagValues))
            if (
                pvalue["type"] == "file list"
                or pvalue["type"] == "directory list"
                or pvalue["type"] == "patternQuery"
            ):
                files = flagValues
                # can pass isFile flag for patternQuery - even if there are directories because the root directory is always given and must be mappable
                isFileFlag = True
                if pvalue["type"] == "directory list":
                    isFileFlag = False
                if files:
                    flags = []
                    baseFlag = ""
                    for fgroup in files:
                        if flagName:
                            baseFlag = flagName
                        else:
                            baseFlag = ""
                        for f in fgroup:
                            hostFile = self.bwbPathToContainerPath(
                                f, isFile=True, returnNone=False
                            )
                            baseFlag += hostFile + " "
                        flags.append(baseFlag)
                    return flags
            elif pvalue["type"][-4:] == "list":
                flags = []

                if flagValues:
                    for fgroup in flagValues:
                        if flagName:
                            baseFlag = flagName
                        else:
                            baseFlag = ""
                        baseFlag += " ".join(fgroup)
                        flags.append(baseFlag + value)
                return flags
        return None

    def findIteratedFlags(self, cmd):
        # replace positional values
        cmd = self.replaceIteratedVars(cmd)
        pattern = r"\_iterate\{([^\}]+)\}"
        regex = re.compile(pattern)
        subs = []
        subFlags = {}
        cmds = []
        sys.stderr.write("command is {}\n".format(cmd))
        maxLen = 0
        # find matches
        for match in regex.finditer(cmd):
            sys.stderr.write("matched {}\n".format(match.group(1)))
            sub = match.group(1)
            if sub not in subs:
                subs.append(sub)
                subFlags[sub] = self.iteratedfString(sub)
                sys.stderr.write("sub is {} flags are {}".format(sub, subFlags[sub]))
                if len(subFlags[sub]) > maxLen:
                    maxLen = len(subFlags[sub])
        sys.stderr.write("subs are {}\n".format(subs))

        if not maxLen:
            return [cmd]
        for i in range(maxLen):
            cmds.append(cmd)
            for sub in subs:
                index = i % len(subFlags[sub])
                sys.stderr.write("sub {} i = {}\n".format(sub, i))
                replaceStr = subFlags[sub][index]
                cmds[i] = cmds[i].replace("_iterate{{{}}}".format(sub), replaceStr)

        breakpoint(message="end of findIterated cmds are {}".format(cmds))
        return cmds

    def replaceIteratedVars(self, cmd):
        # replace any _bwb with _iter if iterated
        pattern = r"\_bwb\{([^\}]+)\}"
        regex = re.compile(pattern)
        subs = []
        for match in regex.finditer(cmd):
            pname = match.group(1)
            if (
                pname
                and self.iterateSettings["iteratedAttrs"]
                and pname in self.iterateSettings["iteratedAttrs"]
            ):
                cmd = cmd.replace(
                    "_bwb{{{}}}".format(pname), "_iterate{{{}}}".format(pname)
                )
        return cmd

    def replaceVars(self, cmd, pnames, varSeen):
        pattern = r"\_bwb\{([^\}]+)\}"
        regex = re.compile(pattern)
        subs = []
        sys.stderr.write("command is {}\n".format(cmd))
        for match in regex.finditer(cmd):
            sys.stderr.write("matched {}\n".format(match.group(1)))
            sub = match.group(1)
            if sub not in subs:
                subs.append(sub)
        for sub in subs:
            # remove front and trailing spaces
            # create match
            matchStr = "_bwb{" + sub + "}"
            psub = sub.strip()
            sys.stderr.write("looking for {} in {}\n".format(psub, pnames))
            if psub in pnames and getattr(self, psub):
                fStr = self.flagString(psub)
                sys.stderr.write("psub {} is  in {}\n".format(psub, pnames))
                sys.stderr.write("fstr {}\n".format(fStr))
                if fStr:
                    sys.stderr.write(
                        "replace matchStr {} in cmd {} with fStr {}\n".format(
                            matchStr, cmd, fStr
                        )
                    )
                    cmd = cmd.replace(matchStr, fStr)
                    sys.stderr.write("to give new cmd {}\n".format(cmd))
                    varSeen[fStr] = True
                else:
                    cmd = cmd.replace(matchStr, "")
            else:
                cmd = cmd.replace(matchStr, "")
        return cmd

    def generateCmdFromBash(self, executables, flags=[], args=[]):
        # by default all flags and arguments are applied to final command
        # use positional _bwb{} variables so specify flags and arguments if there are multiple commands
        # unused args and flags are applied to the final command

        # multi executable commands need to use bash -c 'cmd1 && cmd2' type syntax - note this can cause problems when stopping container
        # can also have no executable in which case we retun nothing

        # will return a list of commands if there are iterables in the command string - it is possible to have iterated variables outside the command string

        # iterable flags is a dict of lists of flags
        iteratedFlags = []
        if not executables:
            return ""
        cmdStr = "bash -c '"
        if len(executables) == 1:
            cmdStr = ""
        varSeen = {}
        sys.stderr.write("execs {}\n".format(executables))
        sys.stderr.write("executables {}\n".format(executables))
        for executable in executables:
            lastExecutable = executable == executables[-1]
            sys.stderr.write(
                "Orig executable {} flags {} args {}\n".format(executable, flags, args)
            )
            pnames = []
            if "parameters" in self.data and self.data["parameters"] is not None:
                pnames = self.data["parameters"].keys()
            executable = self.replaceVars(executable, pnames, varSeen)
            sys.stderr.write(
                "New executable {} flags {} args {}\n".format(executable, flags, args)
            )
            cmdStr += executable + " "
            # extra args and flags go after last executable
            if lastExecutable:
                for flag in flags:
                    if flag not in varSeen:
                        cmdStr += str(flag) + " "
                for arg in args:
                    if arg not in varSeen:
                        cmdStr += str(arg) + " "
                if len(executables) > 1:
                    cmdStr += "'"
            else:
                cmdStr += " && "
        return cmdStr

    def getEnvironmentVariables(self):
        # dynamic environment variables
        if "parameters" in self.data and self.data["parameters"] is not None:
            for pname in self.data["parameters"]:           
                pvalue = self.data["parameters"][pname]
                if "env" in pvalue and getattr(self, pname) is not None:
                    checkAttr = pname + "Checked"
                    setenv = False
                    # check if boolean
                    if pvalue["type"] == "bool":
                        if getattr(self, pname) is True:
                            self.envVars[pvalue["env"]] = getattr(self, pname)
                    else:
                        if pname in self.optionsChecked:
                            if self.optionsChecked[pname]:
                                self.envVars[pvalue["env"]] = self.getAttrValue(pname)
                        else:
                            self.envVars[pvalue["env"]] = self.getAttrValue(pname)


        # now assign static environment variables
        if "env" in self.data:
            for e in self.data["env"]:
                if e not in self.envVars:
                    self.envVars[e] = self.data["env"][e]

    def exportWorkflow(self,widgetName="chosen widget"):
        qm = QtGui.QMessageBox
        ret = qm.question(
                self, "Export?", "Export workflow starting from {}?".format(widgetName), qm.Yes | qm.No
            )
        if ret == qm.No:
            return
        else:
            self.saveBashFile = ""
            retValue = QtWidgets.QFileDialog.getSaveFileName(
                self,
                "Export Docker commands",
                "myscript.sh",
                "Text files (*.sh);;All Files (*)",
            )[0]
            if retValue:
                self.saveBashFile = retValue
                with open(self.saveBashFile, "w") as f:
                    f.write("#!/bin/bash\n")
                self.useTestMode=True
                self.startJob()
                mb = QtGui.QMessageBox
                ret  = mb.information(self,"Export to workflow to script","Saved workflow to {}".format(self.saveBashFile))
                self.useTestMode=False
                sleep(10)
                self.replaceMappingsWithBashVariables()
                return
            else:
                return       
                     
    def replaceMappingsWithBashVariables(self):
        homeDirs=self.dockerClient.bwbMounts.keys()
        with open(self.saveBashFile) as f:
            s = f.read()
        i=1
        for homeDir in homeDirs:
            mapping="-v "+ homeDir
            if mapping in s:
                parm="-v ${}".format(str(i))
                s=s.replace(mapping,parm)
            i= int(i)+1
        with open(self.saveBashFile, 'w') as f:
            f.write(s)
        # change this to something nicer when we start having user permissions for written files
        os.system("chmod +777 {}".format(self.saveBashFile))
        
    # Event handlers
    def onRunClicked(self, button=None):
        if button and self.useTestMode:
            qm = QtGui.QMessageBox
            ret = qm.question(
                self, "Test Run?", "Run without generating results?", qm.Yes | qm.No
            )
            if ret == qm.No:
                return           
        self.startJob()

    def onStopClicked(self):
        self.pConsole.stop("Stopped by user")
        self.setStatusMessage("Stopped")
        self.status = "stopped"
        if self.jobRunning:
            self.bgui.reenableAll(self)
            self.reenableExec()
            self.jobRunning = False

    def onRunFinished(self, code=None, status=None):
        self.jobRunning = False
        self.pConsole.writeMessage("Finished")
        if code is not None:
            if code:
                self.pConsole.writeMessage("Exit code is {}".format(code), color=Qt.red)
            else:
                self.pConsole.writeMessage("Exit code is {}".format(code))
        if status is not None:
            if status:
                self.pConsole.writeMessage("Exit status is {}".format(status), color=Qt.red)
            else:
                self.pConsole.writeMessage("Exit status is {}".format(status))
        self.bgui.reenableAll(self)
        self.reenableExec()
        if self.status != "stopped" and self.status != "finished":
            if status:
                self.setStatusMessage("Error Code {}".format(code))
                self.status = "error"
            else:
                self.setStatusMessage("Finished")
                self.status = "finished"
                if hasattr(self, "handleOutputs"):
                    self.handleOutputs()

    def onRunError(self, error):
        self.bgui.reenableAll(self)
        self.reenableExec()
        self.jobRunning = False
        self.console.writeMessage("Error occurred {}\n".format(error), color=Qt.red)

    def onRunMessage(self, message):
        self.pConsole.writeMessage(message)

    # Utilities
    def bwbPathToContainerPath(self, path, isFile=False, returnNone=False):
        # converts the path entered relative to the bwb mountpoint
        # will return None if not found or the original path (default) depending on variable
        # first map it to the  path
        pathFile = None
        if isFile:
            dirPath = os.path.dirname(path)
            pathFile = os.path.basename(path)
            hostPath = os.path.normpath(
                self.dockerClient.to_best_host_directory(dirPath, returnNone=False)
            )
            sys.stderr.write(
                "dirPath {} pathFile {} hostPath {}\n".format(
                    dirPath, pathFile, hostPath
                )
            )
        else:
            hostPath = os.path.normpath(
                self.dockerClient.to_best_host_directory(path, returnNone=False)
            )
            sys.stderr.write("bwbPath {} hostPath {}\n".format(path, hostPath))
        conPath = None
        # now get all the possible submappings to volumeMappings by comparing the true hostmappings
        # if submap is found convert the common path to the container path
        # return shortest path
        for conVol, bwbVol in self.hostVolumes.items():
            hostVol = os.path.normpath(
                self.dockerClient.to_best_host_directory(bwbVol, returnNone=False)
            )
            sys.stderr.write(
                "checking conVol {} bwbVol {} hostVol {} hostPath {}\n".format(
                    conVol, bwbVol, hostVol, hostPath
                )
            )
            try:
                commonPath=os.path.commonpath([hostVol, hostPath])
                sys.stderr.write(
                    "common {}\n".format(commonPath)
                )
            except Exception as e:
                pass
            prefix = None
            if hostVol == hostPath:
                prefix = ""
            elif Path(hostVol) in Path(hostPath).parents:
                prefix = str(Path(hostPath).relative_to(hostVol))
                sys.stderr.write("prefix is {}\n".format(prefix))
                if prefix == ".":
                    prefix = ""
            if prefix is not None:
                cleanConVol = os.path.normpath(conVol)
                myConPath = os.path.normpath(str.join(os.sep, (cleanConVol, prefix)))
                if conPath is None or len(myConPath) < len(conPath):
                    conPath = myConPath

        if conPath is not None:
            if isFile:
                return os.path.normpath(str.join(os.sep, (conPath, pathFile)))
            return conPath
        else:
            if returnNone:
                return conPath
            return path

    def findTopDirectory(self, files):
        # return the top directory in a set of files
        # actually returns the shortest directory
        bestPath = None
        for f in files:
            bwbVol = os.path.dirname(os.path.normpath(f))
            if bestPath is None:
                bestPath = bwbVol
            elif len(bwbVol) < len(bestPath):
                bestPath = bwbVol
        return bestPath

    def disableExec(self):
        self.btnRun.setText("Running")
        self.btnRun.setEnabled(False)
        self.btnStop.setEnabled(True)
        self.cboRunMode.setEnabled(False)
        self.graphicsMode.setEnabled(False)

    def reenableExec(self):
        self.btnRun.setText("Start")
        self.btnRun.setEnabled(True)
        self.btnStop.setEnabled(False)
        self.cboRunMode.setEnabled(True)
        self.graphicsMode.setEnabled(True)
