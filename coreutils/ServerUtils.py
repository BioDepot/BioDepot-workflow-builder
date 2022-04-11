import os
import re
import sys
import json
import jsonpickle
import pickle
import csv
import tempfile, shutil
import OWImageBuilder
import workflowTools
import socket
import copy
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
from PyQt5.QtWidgets import *
from PyQt5.QtGui import QStandardItemModel, QTableView, QDropEvent

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


class IterateDialog(QDialog):
    def __init__(self, iterateSettings):
        nRows = len(iterateSettings["iterableAttrs"])
        if not nRows:
            return
        super().__init__()
        self.css = """
        QPushButton {background-color: #1588c5; color: white; height: 20px; border: 1px solid black; border-radius: 2px;}
        QPushButton:hover {background-color: #1555f5; }
        QPushButton:hover:pressed { background-color: #1588c5; color: black; border-style: inset; border: 1px solid white} 
        QPushButton:disabled { background-color: lightGray; border: 1px solid gray; } 
        """
        defaults={
            "groupSize": "1",
            "threads": "1",
            "ram": "0",
            "splitType": "None",
            "splitSize": "0",
            "splitUnit": "Auto",
            "splitCmd": "Auto",
            "mergeCmd": "None",
            "mergeOutput": "None"
        }
        nCols = len(defaults)+2
        self.setMinimumSize(860, 240)
        self.iterateSettings = iterateSettings
        self.setWindowTitle("Edit iterate settings")
        self.table = QTableWidget()
        self.table.setColumnCount(nCols)
        for col in range(nCols - 1):
            self.table.horizontalHeader().setResizeMode(
                col, QtGui.QHeaderView.ResizeToContents
            )
        self.table.horizontalHeader().setResizeMode(
            nCols - 1, QtGui.QHeaderView.Stretch
        )
        self.table.horizontalHeader().setStretchLastSection(True)
        self.table.setHorizontalHeaderLabels(
            ["", "Parameter", "Group size", "Threads", "RAM (MB)","Split type", "Split size","Split unit","Split cmd", "Merge cmd", "Merge output"]
        )
        self.settingsCopy = copy.deepcopy(self.iterateSettings)
        self.table.setRowCount(nRows)
        rowNum = 0
        self.splitTypeOptions=['None','Auto','Count','Max file size','Min file size']
        self.splitUnitOptionsCounts=['None','Auto','By lines','By slices','By reads','By bytes','By kB','By MB']
        self.splitUnitOptionsValues=['None','Auto','lines','slices','reads','bytes','kB','MB']
        self.splitCmdOptions=['None','Auto','For list','For image','For fastq','For SAM','For BAM']
        self.mergeCmdOptions=['None','Auto','Concatenate','For image','For fastq','For SAM','For BAM']
        self.mergeOutputOptions=['None','Auto']
        if self.iterateSettings["outputableAttrs"]:
            self.mergeOutputOptions.extend(self.iterateSettings["outputableAttrs"])
        for parm in self.settingsCopy["iterableAttrs"]:
            splittable=False 
            if parm  in self.iterateSettings["splittableAttrs"]:
                splittable=True
            # init values if they are not there
            if "data" not in self.settingsCopy:
                self.settingsCopy["data"] = {}
            if parm not in self.settingsCopy["data"]:
                self.settingsCopy["data"][parm] = defaults
            for key in defaults:
                if (
                    key not in self.settingsCopy["data"][parm].keys()
                    or not self.settingsCopy["data"][parm][key]
                ):
                    self.settingsCopy["data"][parm][key] = defaults[key]
             # make column items
             
            cb = QTableWidgetItem()
            parmItem = QTableWidgetItem(parm)
            groupSizeItem = QTableWidgetItem(
                self.settingsCopy["data"][parm]["groupSize"]
            )
            threadItem = QTableWidgetItem(self.settingsCopy["data"][parm]["threads"])
            ramItem = QTableWidgetItem(self.settingsCopy["data"][parm]["ram"])
            splitTypeItem = QTableWidgetItem(self.settingsCopy["data"][parm]["splitType"])
            splitSizeItem = QTableWidgetItem(self.settingsCopy["data"][parm]["splitSize"])
            splitUnitItem = QTableWidgetItem(self.settingsCopy["data"][parm]["splitUnit"])
            splitCmdItem = QTableWidgetItem(self.settingsCopy["data"][parm]["splitCmd"])
            mergeCmdItem = QTableWidgetItem(self.settingsCopy["data"][parm]["mergeCmd"])
            mergeOutputItem = QTableWidgetItem(self.settingsCopy["data"][parm]["mergeOutput"])

            self.setSelect(parmItem, False)
            parmItem.setFlags(parmItem.flags() ^ Qt.ItemIsEditable)
            parmItem.setFlags(parmItem.flags() ^ Qt.ItemIsSelectable)
            cb.setFlags(cb.flags() ^ Qt.ItemIsSelectable)
            if (
                "iteratedAttrs" in self.settingsCopy
                and parm in self.settingsCopy["iteratedAttrs"]
            ):
                if not splittable:
                    self.setEnableSelect(splitTypeItem, False)
                    self.setEnableSelect(splitSizeItem, False)
                    self.setEnableSelect(splitUnitItem, False)
                    self.setEnableSelect(splitCmdItem, False)
                    self.setEnableSelect(mergeCmdItem, False)
                    self.setEnableSelect(mergeOutputItem, False)
                else:
                    self.enableComboBox(rowNum,5,self.splitTypeOptions,splitTypeItem.text())
                    self.setEnableSelect(splitSizeItem, True)
                    self.enableComboBox(rowNum,7,self.splitUnitOptionsValues,splitUnitItem.text())
                    self.enableComboBox(rowNum,8,self.splitCmdOptions,splitCmdItem.text())
                    self.enableComboBox(rowNum,9,self.mergeCmdOptions,mergeCmdItem.text())
                    self.enableComboBox(rowNum,10,self.mergeOutputOptions,mergeOutputItem.text())
                cb.setCheckState(QtCore.Qt.Checked)
            else:
                cb.setCheckState(QtCore.Qt.Unchecked)
                self.setEnable(parmItem, False)
                self.setEnableSelect(groupSizeItem, False)
                self.setEnableSelect(threadItem, False)
                self.setEnableSelect(ramItem, False)
                self.setEnableSelect(splitTypeItem, False)
                self.setEnableSelect(splitSizeItem, False)
                self.setEnableSelect(splitUnitItem, False)
                self.setEnableSelect(splitCmdItem, False)
                self.setEnableSelect(mergeCmdItem, False)
                self.setEnableSelect(mergeOutputItem, False)
            
            self.table.setItem(rowNum, 0, cb)
            self.table.setItem(rowNum, 1, parmItem)
            self.table.setItem(rowNum, 2, groupSizeItem)
            self.table.setItem(rowNum, 3, threadItem)
            self.table.setItem(rowNum, 4, ramItem)
            self.table.setItem(rowNum, 5, splitTypeItem)
            self.table.setItem(rowNum, 6, splitSizeItem)
            self.table.setItem(rowNum, 7, splitUnitItem)
            self.table.setItem(rowNum, 8, splitCmdItem)
            self.table.setItem(rowNum, 9, mergeCmdItem)
            self.table.setItem(rowNum, 10, mergeOutputItem)
            rowNum = rowNum + 1
        self.table.cellChanged.connect(self.onCellChange)            
        
        # buttons for save and load

        saveBtn = gui.button(None, self, "Save", callback=self.save)
        saveBtn.setStyleSheet(self.css)
        saveBtn.setFixedSize(70, 20)

        # table
        tableBox = QGroupBox()
        scroll_area = QScrollArea(verticalScrollBarPolicy=Qt.ScrollBarAlwaysOn)
        scroll_area.setWidget(tableBox)
        scroll_area.setWidgetResizable(True)
        tableLayout = QHBoxLayout()
        tableLayout.addWidget(self.table)
        tableBox.setLayout(tableLayout)

        # buttons
        buttonLayout = QHBoxLayout()
        buttonLayout.setAlignment(Qt.AlignTop)
        buttonLayout.addWidget(saveBtn)

        iterateLayout = QVBoxLayout()
        iterateLayout.addWidget(tableBox)
        iterateLayout.addLayout(buttonLayout)
        self.setLayout(iterateLayout)

    def setEnable(self, item, state):
        if state:
            item.setFlags(item.flags() | Qt.ItemIsEnabled)
        else:
            item.setFlags(item.flags() ^ Qt.ItemIsEnabled)

    def setSelect(self, item, state):
        if state:
            item.setFlags(item.flags() | Qt.ItemIsSelectable)
        else:
            item.setFlags(item.flags() ^ Qt.ItemIsSelectable)

    def setEnableSelect(self, item, state):
        if not item:
            return
        if state:
            item.setFlags(item.flags() | Qt.ItemIsEnabled)
            item.setFlags(item.flags() | Qt.ItemIsSelectable)
        else:
            item.setFlags(item.flags() & ~Qt.ItemIsEnabled)
            item.setFlags(item.flags() & ~Qt.ItemIsSelectable)

    def onCellChange(self, row, column):
        if column == 0 :
            return self.onCheckBoxChange(row);
        if column == 5:
            return self.onSplitTypeChange(row);
        if column > 6 :
            if self.table.cellWidget(row,column):
                if self.table.item(row,column).text() != self.table.cellWidget(row,column).currentText():
                    self.table.item(row,column).setText(self.table.cellWidget(row,column).currentText())
    def onSplitTypeChange(self, row):
        parmItem = self.table.item(row, 1)
        parm = parmItem.text()
        if self.table.cellWidget(row,5):
            self.table.item(row,5).setText(self.table.cellWidget(row,5).currentText())
            if self.table.cellWidget(row,5).currentText() == "None":
                sys.stderr.write("inside if \n");
                self.setEnableSelect(self.table.item(row, 6), False)
                for col in range(7, 9):
                    if parm  in self.iterateSettings["splittableAttrs"]:
                        item = self.table.item(row, col)
                        self.setEnableSelect(item, False)
                        self.disableBox(row,col)
            else:
                self.setEnableSelect(self.table.item(row, 6), True)
                for col in range(7, 9):
                    if col == 7:
                        #find the value
                        value=self.table.item(row, col).text()
                        currentIndex=0
                        currentCount=0
                        currentTable=self.splitUnitOptionsValues
                        if self.table.cellWidget(row,col):
                            currentIndex=self.table.cellWidget(row,col).currentIndex()
                            currentCount=self.table.cellWidget(row,col).count()
                        if self.table.cellWidget(row,5).currentText() == "Count":
                            currentTable=self.splitUnitOptionsCounts
                        if currentIndex < currentCount:
                            value=currentTable[currentIndex]
                        self.table.item(row, col).setText(value)
                        self.enableComboBox(row,col,currentTable,value)
                    if col == 8:
                        self.enableComboBox(row,col,self.splitCmdOptions,self.table.item(row, col).text())                    

      
    def onCheckBoxChange(self, row):
        cb = self.table.item(row, 0)
        parmItem = self.table.item(row, 1)
        parm = parmItem.text()
        #also need to enable disable any associated cell widgets
        if cb.checkState() == QtCore.Qt.Checked:
            self.setEnable(parmItem, True)
            for col in range(2, 5):
                item = self.table.item(row, col)
                self.setEnableSelect(item, True)
            for col in range(5, 11):
                if parm  in self.iterateSettings["splittableAttrs"]:
                    item = self.table.item(row, col)
                    self.setEnableSelect(item, True)
                    if col == 5:
                        self.enableComboBox(row,col,self.splitTypeOptions,self.table.item(row, col).text())
                    if self.table.item(row, 5).text() and self.table.item(row, 5).text() != "None":
                        if col == 6:   
                            self.setEnableSelect(item, True)
                        if col == 7:
                            self.enableComboBox(row,col,self.splitUnitOptions,self.table.item(row, col).text())
                        if col == 8:
                            self.enableComboBox(row,col,self.splitCmdOptions,self.table.item(row, col).text())
                    elif col > 5 and col < 9:
                        self.setEnableSelect(item, False)
                    else:
                        if col == 9:
                            self.enableComboBox(row,col,self.mergeCmdOptions,self.table.item(row, col).text())
                        if col == 10:
                            self.enableComboBox(row,col,self.mergeOutputOptions,self.table.item(row, col).text())                    
                    
            if (
                "iteratedAttrs" in self.settingsCopy
                and parm not in self.settingsCopy["iteratedAttrs"]
            ):
                self.settingsCopy["iteratedAttrs"].append(parm)
        else:
            self.setEnable(parmItem, False)
            for col in range(2, 5):
                item = self.table.item(row, col)
                self.setEnableSelect(item, False)
            for col in range(5, 11):
                item = self.table.item(row, col)
                self.setEnableSelect(item, False)
                self.disableBox(row,col)                
            if (
                "iteratedAttrs" in self.settingsCopy
                and parm not in self.settingsCopy["iteratedAttrs"]
            ):
                self.settingsCopy["iteratedAttrs"].remove(parm)

    def enableComboBox(self,row,col,comboOptions,value=None):
        comboBox = QtGui.QComboBox()
        comboBox.setEditable(False)
        for comboOption in comboOptions:
            comboBox.addItem(comboOption)
        if value:
            if value not in comboOptions:
                comboBox.addItem(value);
                comboBox.setEditable(True)
            else:
                comboBox.addItem("Custom")
            comboBox.setCurrentText(value)
        else:
            comboBox.addItem("Custom")
        if value == "Custom":
            comboBox.setEditable(True)
        
        comboBox.currentIndexChanged.connect(lambda index: self.onCbIndexChange(index,cb=comboBox,r=row,c=col))
        self.table.setCellWidget(row,col,comboBox)
        if self.table.item(row,col).text() != comboBox.currentText():
            self.table.item(row,col).setText(comboBox.currentText())
    def onCbIndexChange (self,index,cb=None,r=None,c=None):
        if cb and cb.count()-1 == index:
            cb.setEditable(True)
        else:
            cb.setEditable(False)
        if cb and r is not None and c is not None:
            if self.table.item(r,c).text() != cb.currentText():
                self.table.item(r,c).setText(cb.currentText())   
         
    def disableBox(self,row,col):
        if self.table.cellWidget(row,col):
            self.table.removeCellWidget(row,col);

    def closeEvent(self, event):
        if self.iterateSettings == self.settingsCopy:
            event.accept()
            return
        qm = QMessageBox(self)
        qm.setWindowTitle("Save input")
        qm.setInformativeText("Save input to current settings?")
        qm.setStandardButtons(QMessageBox.Yes | QMessageBox.No | QMessageBox.Cancel)
        qm.setDefaultButton(QMessageBox.Cancel)
        reply = qm.exec_()
        if reply == QMessageBox.Yes:
            self.save()
        elif reply == QMessageBox.No:
            event.accept()
        else:
            event.ignore()
        return

    def updateSettings(self):
        newData = OrderedDict()
        newIteratedAttrs = []

        for i in range(self.table.rowCount()):
            #update any comboBoxes
            self.onCellChange(i,6)
            self.onCellChange(i,7)
            parm = itemToText(self.table.item(i, 1))
            newData[parm] = {
                "groupSize": itemToText(self.table.item(i, 2)),
                "threads": itemToText(self.table.item(i, 3)),
                "ram": itemToText(self.table.item(i, 4)),
                "splitSize": itemToText(self.table.item(i, 5)),
                "splitUnit" :itemToText(self.table.item(i, 6)),
                "splitCmd": itemToText(self.table.item(i, 7)),
                "mergeCmd": itemToText(self.table.item(i, 8)),
                "mergeOutput": itemToText(self.table.item(i, 9)),
            }
            if self.table.item(i, 0).checkState() == QtCore.Qt.Checked:
                newIteratedAttrs.append(parm)
        self.settingsCopy["data"] = newData
        self.settingsCopy["iteratedAttrs"] = newIteratedAttrs

    def save(self):
        self.updateSettings()
        try:
            self.iterateSettings = self.settingsCopy
            title = "Save settings"
            message = "Settings successfully saved"
            ret = QtGui.QMessageBox.information(
                self, title, message, QtGui.QMessageBox.Ok
            )
        except Exception as e:
            warning = QtGui.QMessageBox.warning(
                None, "", "Settings not saved - error: {}\n".format(str(e))
            )


class ServerDialog(QDialog):
    def __init__(self, serverSettings):
        super().__init__()
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
        self.setMinimumSize(720, 384)
        self.serverSettings = serverSettings
        self.setWindowTitle("Edit server settings")
        self.settingsFile = "/biodepot/serverSettings.json"
        self.addIcon = QtGui.QIcon("/icons/add.png")
        self.removeIcon = QtGui.QIcon("/icons/remove.png")
        self.table = TableWidgetDragRows()
        self.table.setColumnCount(3)
        self.table.horizontalHeader().setResizeMode(
            0, QtGui.QHeaderView.ResizeToContents
        )
        self.table.horizontalHeader().setResizeMode(
            1, QtGui.QHeaderView.ResizeToContents
        )
        self.table.horizontalHeader().setResizeMode(2, QtGui.QHeaderView.Stretch)
        self.table.horizontalHeader().setStretchLastSection(True)
        self.table.setHorizontalHeaderLabels(["Server IP", "Threads", "Volume map"])
        # start with blank slate - and keep temp copy of serverSettings
        # only update when dialog is done and if the
        self.serverSettingsCopy = {}
        try:
            self.importServers(filename=self.settingsFile)
        except Exception as e:
            self.serverSettingsCopy = {}

        if "data" not in self.serverSettingsCopy:
            self.serverSettingsCopy["data"] = OrderedDict()
        nRows = len(self.serverSettingsCopy["data"].keys())
        if nRows:
            self.table.setRowCount(nRows)
            rowNum = 0
            for addr in self.serverSettingsCopy["data"]:
                # sys.stderr.write('addr {} data {}\n'.format(addr,serverSettingsCopy['data'][addr]))
                self.table.setItem(rowNum, 0, QTableWidgetItem(addr))
                self.table.setItem(
                    rowNum,
                    1,
                    QTableWidgetItem(
                        self.serverSettingsCopy["data"][addr]["maxThreads"]
                    ),
                )
                self.table.setItem(
                    rowNum,
                    2,
                    QTableWidgetItem(
                        self.serverSettingsCopy["data"][addr]["volumeMapping"]
                    ),
                )
                rowNum = rowNum + 1

        # buttons for save and load

        loadBtn = gui.button(None, self, "Load", callback=self.importServers)
        loadBtn.setStyleSheet(self.css)
        loadBtn.setFixedSize(70, 20)

        saveBtn = gui.button(None, self, "Save", callback=self.save)
        saveBtn.setStyleSheet(self.css)
        saveBtn.setFixedSize(70, 20)

        saveAsBtn = gui.button(None, self, "Save As", callback=self.exportServers)
        saveAsBtn.setStyleSheet(self.css)
        saveAsBtn.setFixedSize(70, 20)

        addBtn = gui.button(None, self, "", callback=self.addRow)
        removeBtn = gui.button(None, self, "", callback=self.removeRow)

        addBtn.setIcon(self.addIcon)
        addBtn.setStyleSheet(self.addRemoveCSS)
        removeBtn.setIcon(self.removeIcon)
        removeBtn.setStyleSheet(self.addRemoveCSS)

        # table
        tableBox = QGroupBox()
        scroll_area = QScrollArea(verticalScrollBarPolicy=Qt.ScrollBarAlwaysOn)
        scroll_area.setWidget(tableBox)
        scroll_area.setWidgetResizable(True)
        tableLayout = QHBoxLayout()
        tableLayout.addWidget(self.table)
        tableBox.setLayout(tableLayout)

        # buttons
        buttonLayout = QHBoxLayout()
        buttonLayout.addWidget(loadBtn)
        buttonLayout.addWidget(saveBtn)
        buttonLayout.addWidget(saveAsBtn)
        buttonLayout.addStretch(1)
        buttonLayout.addWidget(addBtn)
        buttonLayout.addWidget(removeBtn)

        serverLayout = QVBoxLayout()
        serverLayout.addWidget(tableBox)
        serverLayout.addLayout(buttonLayout)
        self.setLayout(serverLayout)

        removeBtn.setEnabled(bool(len(self.table.selectionModel().selectedRows())))
        self.table.itemSelectionChanged.connect(
            lambda: removeBtn.setEnabled(
                bool(len(self.table.selectionModel().selectedRows()))
            )
        )

    def redrawTable(self):
        if "data" not in self.serverSettingsCopy or not self.serverSettingsCopy["data"]:
            return
        nRows = len(self.serverSettingsCopy["data"].keys())
        if nRows:
            self.table.setRowCount(nRows)
            rowNum = 0
            for addr in self.serverSettingsCopy["data"]:
                # sys.stderr.write('addr {} data {}\n'.format(addr,self.serverSettingsCopy['data'][addr]))
                self.table.setItem(rowNum, 0, QTableWidgetItem(addr))
                self.table.setItem(
                    rowNum,
                    1,
                    QTableWidgetItem(
                        self.serverSettingsCopy["data"][addr]["maxThreads"]
                    ),
                )
                self.table.setItem(
                    rowNum,
                    2,
                    QTableWidgetItem(
                        self.serverSettingsCopy["data"][addr]["volumeMapping"]
                    ),
                )
                rowNum = rowNum + 1

    def closeEvent(self, event):
        if self.serverSettings == self.serverSettingsCopy:
            event.accept()
            return
        qm = QMessageBox(self)
        qm.setWindowTitle("Save input")
        qm.setInformativeText("Save input to current settings?")
        qm.setStandardButtons(QMessageBox.Yes | QMessageBox.No | QMessageBox.Cancel)
        qm.setDefaultButton(QMessageBox.Cancel)
        reply = qm.exec_()
        if reply == QMessageBox.Yes:
            self.save()
        elif reply == QMessageBox.No:
            event.accept()
        else:
            event.ignore()
        return

    def updateSettings(self):
        newData = OrderedDict()
        for i in range(self.table.rowCount()):
            addr = itemToText(self.table.item(i, 0))
            newData[addr] = {
                "maxThreads": itemToText(self.table.item(i, 1)),
                "volumeMapping": itemToText(self.table.item(i, 2)),
            }
        self.serverSettingsCopy["data"] = newData

    def parseServers(self, filename):
        self.serverSettingsCopy["data"] = OrderedDict()
        with open(filename, "r") as f:
            data = f.read()
        if data.lstrip()[0] == "{":
            # json file
            try:
                temp = jsonpickle.decode(data)
                if temp:
                    self.serverSettingsCopy = temp
                return True
            except Exception as e:
                return False
        elif checkIP(data.lstrip().strip("\t")):
            # tsv file
            lines = data.splitlines()
            for line in lines:
                addr, maxThreads, volumeMapping = line.split("\t")
                if checkIP(addr) and addr not in self.serverSettingsCopy["data"]:
                    self.serverSettingsCopy["data"][addr] = {
                        "maxThreads": maxThreads,
                        "volumeMapping": volumeMapping,
                    }
            return True
        else:
            ret = QtGui.QMessageBox.warning(
                None, "", "Unable to detect valid json or tsv servers file"
            )
            return False

    def save(self):
        self.updateSettings()
        try:
            self.serverSettings = self.serverSettingsCopy
            with open(self.settingsFile, "w") as f:
                f.write(jsonpickle.encode(self.serverSettings))
            title = "Save settings"
            message = "Settings successfully saved"
            ret = QtGui.QMessageBox.information(
                self, title, message, QtGui.QMessageBox.Ok
            )
        except Exception as e:
            warning = QtGui.QMessageBox.warning(
                None, "", "Settings not saved - error: {}\n".format(str(e))
            )

    def importServers(self, filename=None):
        startDir = "/"
        if (
            "saveIPDir" in self.serverSettingsCopy
            and self.serverSettingsCopy["saveIPDir"]
        ):
            startDir = self.serverSettingsCopy["saveIPDir"]
        elif os.path.isdir("/data"):
            startDir = "/data"
        if not filename:
            filename = QFileDialog.getOpenFileName(
                self, "Load IPs from File", startDir
            )[0]
            if not filename:
                return QFileDialog.Rejected
        self.serverSettingsCopy["saveIPDir"] = os.path.dirname(filename)
        if self.parseServers(filename):
            self.redrawTable()

    def exportServers(self):
        startDir = "/"
        if (
            "saveIPDir" in self.serverSettingsCopy
            and self.serverSettingsCopy["saveIPDir"]
        ):
            startDir = self.serverSettingsCopy["saveIPDir"]
        elif os.path.isdir("/data"):
            startDir = "/data"
        filename, filter = QFileDialog.getSaveFileName(
            self, "Save servers", startDir, "tab delimited (*.tsv);; json (*.json)"
        )
        if not filename:
            return QFileDialog.Rejected
        self.serverSettingsCopy["saveIPDir"] = os.path.dirname(filename)
        saveJson = False
        fileStub, extension = os.path.splitext(filename)
        sys.stderr.write("extension of filename {} is {}\n".format(filename, extension))
        if (
            extension == ".json"
            or extension == ".JSON"
            or extension == ".jsn"
            or extension == ".JSN"
        ):
            saveJson = True
        elif extension != ".tsv":
            qm = QtGui.QMessageBox
            title = "Which save format"
            ret = qm.question(
                self, title, "Can save as tsv or json. Save as json?", qm.Yes | qm.No
            )
            if ret == qm.Yes:
                saveJson = True
        with open(filename, "w") as f:
            self.updateSettings()
            if saveJson:
                f.write(jsonpickle.encode(self.serverSettingsCopy))
            else:
                for addr in self.serverSettingsCopy["data"]:
                    if checkIP(addr):
                        output = [addr]
                        if "maxThreads" in self.serverSettingsCopy["data"][addr]:
                            output.append(
                                self.serverSettingsCopy["data"][addr]["maxThreads"]
                            )
                        else:
                            output.append("")
                        if "volumeMapping" in self.serverSettingsCopy["data"][addr]:
                            output.append(
                                self.serverSettingsCopy["data"][addr]["volumeMapping"]
                            )
                        else:
                            output.append("")

    def addRow(self):
        self.table.insertRow(self.table.rowCount())

    def removeRow(self):
        selections = self.table.selectionModel().selectedRows()
        indices = []
        for selection in selections:
            row = selection.row()
            if row not in indices:
                indices.append(row)
        for i in sorted(indices, reverse=True):
            self.table.removeRow(i)


# From https://stackoverflow.com/questions/26227885/drag-and-drop-rows-within-qtablewidget
class TableWidgetDragRows(QTableWidget):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        self.setDragEnabled(True)
        self.setAcceptDrops(True)
        self.viewport().setAcceptDrops(True)
        self.setDragDropOverwriteMode(False)
        self.setDropIndicatorShown(True)

        self.setSelectionMode(QAbstractItemView.ExtendedSelection)
        self.setSelectionBehavior(QAbstractItemView.SelectRows)
        self.setDragDropMode(QAbstractItemView.InternalMove)

    def dropEvent(self, event):
        if not event.isAccepted() and event.source() == self:
            drop_row = self.drop_on(event)

            rows = sorted(set(item.row() for item in self.selectedItems()))
            rows_to_move = [
                [
                    QTableWidgetItem(self.item(row_index, column_index))
                    for column_index in range(self.columnCount())
                ]
                for row_index in rows
            ]
            for row_index in reversed(rows):
                self.removeRow(row_index)
                if row_index < drop_row:
                    drop_row -= 1

            for row_index, data in enumerate(rows_to_move):
                row_index += drop_row
                self.insertRow(row_index)
                for column_index, column_data in enumerate(data):
                    self.setItem(row_index, column_index, column_data)
            event.accept()
            for row_index in range(len(rows_to_move)):
                self.item(drop_row + row_index, 0).setSelected(True)
                self.item(drop_row + row_index, 1).setSelected(True)
        super().dropEvent(event)

    def drop_on(self, event):
        index = self.indexAt(event.pos())
        if not index.isValid():
            return self.rowCount()

        return index.row() + 1 if self.is_below(event.pos(), index) else index.row()

    def is_below(self, pos, index):
        rect = self.visualRect(index)
        margin = 2
        if pos.y() - rect.top() < margin:
            return False
        elif rect.bottom() - pos.y() < margin:
            return True
        # noinspection PyTypeChecker
        return (
            rect.contains(pos, True)
            and not (int(self.model().flags(index)) & Qt.ItemIsDropEnabled)
            and pos.y() >= rect.center().y()
        )


def checkIP(addr):
    if not addr:
        return False
    try:
        socket.inet_aton(addr)
        return True
    except socket.error:
        return False


def itemToText(item):
    if item:
        return item.text()
    return None


def registerDirectory(baseToolPath):
    os.system("cd {} && pip install -e .".format(baseToolPath))


def editIPs(serverSettings):
    serverDialog = ServerDialog(serverSettings)
    serverDialog.exec_()
    serverSettings = serverDialog.serverSettings
