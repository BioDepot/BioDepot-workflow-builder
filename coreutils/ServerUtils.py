import os
import re
import sys
import json
import jsonpickle
import pickle
import csv
import tempfile,shutil
import OWImageBuilder
import workflowTools
import socket
from xml.dom import minidom
from glob import glob
from pathlib import Path
from shutil import copyfile
from createWidget import mergeWidget, createWidget,findIconFile
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
    QWidget, QButtonGroup, QGroupBox, QRadioButton, QSlider,
    QDoubleSpinBox, QComboBox, QSpinBox, QListView, QLabel,
    QScrollArea, QVBoxLayout, QHBoxLayout, QFormLayout,
    QSizePolicy, QApplication, QCheckBox
)

defaultIconFile='/icons/default.png'

class ServerDialog(QDialog):
    def __init__(self, serverSettings):
        super().__init__()
        self.css = '''
        QPushButton {background-color: #1588c5; color: white; height: 20px; border: 1px solid black; border-radius: 2px;}
        QPushButton:hover {background-color: #1555f5; }
        QPushButton:hover:pressed { background-color: #1588c5; color: black; border-style: inset; border: 1px solid white} 
        QPushButton:disabled { background-color: lightGray; border: 1px solid gray; } 
        '''
        self.addRemoveCSS='''
        QPushButton {background-color: lightBlue; color: white; height: 20px; border: 1px solid black; border-radius: 2px;}
        QPushButton:hover {background-color: blue; }
        QPushButton:hover:pressed { background-color: lightBlue; color: black; border-style: inset; border: 1px solid white} 
        QPushButton:disabled { background-color: white; border: 1px solid gray; } 
        '''   
        self.addIcon=QtGui.QIcon('/icons/add.png')
        self.removeIcon=QtGui.QIcon('/icons/remove.png')
        
        self.serverSettings=serverSettings
        if 'data' not in self.serverSettings:
            self.serverSettings['data']=OrderedDict()
        self.table=TableWidgetDragRows()
        self.table.setColumnCount(3)
        self.table.horizontalHeader().setStretchLastSection(True)
        self.table.setHorizontalHeaderLabels(['Server IP', 'Threads', 'Volume map'])
        nRows=len(serverSettings['data'].keys())
        if nRows:
            self.table.setRowCount(nRows)
            rowNum=0
            for addr in serverSettings['data']:
            #sys.stderr.write('addr {} data {}\n'.format(addr,serverSettings['data'][addr])) 
                self.table.setItem(rowNum,0,QTableWidgetItem(addr))
                self.table.setItem(rowNum,1,QTableWidgetItem(serverSettings['data'][addr]['maxThreads']))
                self.table.setItem(rowNum,2,QTableWidgetItem(serverSettings['data'][addr]['volumeMapping']))
                rowNum=rowNum+1

        #buttons for save and load
        
        loadBtn = gui.button(None, self, "Load", callback=self.importServers)
        loadBtn.setStyleSheet(self.css)
        loadBtn.setFixedSize(70,20)
        
        saveBtn = gui.button(None, self, "Save", callback=self.exportServers)
        saveBtn.setStyleSheet(self.css)
        saveBtn.setFixedSize(70,20)
        
        addBtn=gui.button(None, self, "", callback=self.addRow)
        removeBtn=gui.button(None, self, "", callback=self.removeRow)
        
        addBtn.setIcon(self.addIcon)
        addBtn.setStyleSheet(self.addRemoveCSS)
        removeBtn.setIcon(self.removeIcon)
        removeBtn.setStyleSheet(self.addRemoveCSS)
                
        #table
        tableBox=QGroupBox()
        scroll_area = QScrollArea(verticalScrollBarPolicy=Qt.ScrollBarAlwaysOn)
        scroll_area.setWidget(tableBox)
        scroll_area.setWidgetResizable(True)
        tableLayout=QHBoxLayout()
        tableLayout.addWidget(self.table)
        tableBox.setLayout(tableLayout)
        
        #buttons
        buttonLayout=QHBoxLayout()
        buttonLayout.addWidget(loadBtn)
        buttonLayout.addWidget(saveBtn)
        buttonLayout.addStretch(1)
        buttonLayout.addWidget(addBtn)
        buttonLayout.addWidget(removeBtn)

        serverLayout=QVBoxLayout()
        serverLayout.addWidget(tableBox)
        serverLayout.addLayout(buttonLayout)
        self.setLayout(serverLayout)
        

        removeBtn.setEnabled(bool(len(self.table.selectionModel().selectedRows())))
        self.table.itemSelectionChanged.connect(lambda: removeBtn.setEnabled(bool(len(self.table.selectionModel().selectedRows()))))
        # boxEdit.itemSelectionChanged.connect(lambda: self.onListWidgetSelect(boxEdit,addBtn,removeBtn,lineItem))
        # boxEdit.itemMoved.connect(lambda oldRow,newRow : self.onItemMoved(oldRow,newRow,boxEdit))
        
    def redrawTable(self):
        if not self.serverSettings['data']:
            return
        nRows=len(self.serverSettings['data'].keys())
        if nRows:
            self.table.setRowCount(nRows)
            rowNum=0
            for addr in self.serverSettings['data']:
            #sys.stderr.write('addr {} data {}\n'.format(addr,self.serverSettings['data'][addr])) 
                self.table.setItem(rowNum,0,QTableWidgetItem(addr))
                self.table.setItem(rowNum,1,QTableWidgetItem(self.serverSettings['data'][addr]['maxThreads']))
                self.table.setItem(rowNum,2,QTableWidgetItem(self.serverSettings['data'][addr]['volumeMapping']))
                rowNum=rowNum+1
                        
    def closeEvent(self, event):
        self.updateTable()
        event.accept()
        
    def updateTable(self):
        newData=OrderedDict()
        for i in range(self.table.rowCount()):
            addr=itemToText(self.table.item(i,0))
            newData[addr]={'maxThreads':itemToText(self.table.item(i,1)),'volumeMapping':itemToText(self.table.item(i,2))}
            sys.stderr.write('row {}: {} {} {}\n'.format(i,itemToText(self.table.item(i,0)),itemToText(self.table.item(i,1)),itemToText(self.table.item(i,2))))
        self.serverSettings['data']=newData
        
    def parseServers(self,filename):
        self.serverSettings['data']=OrderedDict()
        with open (filename,'r') as f:
            for line in f:
                line=line.rstrip("\n\r")
                addr,maxThreads,volumeMapping=line.split('\t')
                if checkIP(addr) and addr not in self.serverSettings['data']:
                    self.serverSettings['data'][addr]={'maxThreads':maxThreads,'volumeMapping':volumeMapping}
    def importServers(self):
        startDir='/'
        if 'saveIPDir' in self.serverSettings and self.serverSettings['saveIPDir']:
            startDir=self.serverSettings['saveIPDir']
        elif os.path.isdir('/data'):
            startDir='/data'
        filename = QFileDialog.getOpenFileName(self,"Load IPs from File",startDir)[0]
        if not filename:
            return QFileDialog.Rejected
        self.serverSettings['saveIPDir']=os.path.dirname(filename)
        self.parseServers(filename)
        self.redrawTable()  
              
    def exportServers(self):
        startDir='/'
        if 'saveIPDir' in self.serverSettings and self.serverSettings['saveIPDir']:
            startDir=self.serverSettings['saveIPDir']
        elif os.path.isdir('/data'):
            startDir='/data'
        filename, filter = QFileDialog.getSaveFileName(self, 'Save servers', startDir, "tab delimited (*.tsv)")
        if not filename:
            return QFileDialog.Rejected
        self.serverSettings['saveIPDir']=os.path.dirname(filename)
        with open (filename,'w') as f:
            for i in range(self.table.rowCount()):
                addr=itemToText(self.table.item(i,0))
                if checkIP(addr):
                    output=[addr,itemToText(self.table.item(i,1)),itemToText(self.table.item(i,2))]
                f.write('\t'.join(output)+'\n')
    def addRow(self):
        self.table.insertRow(self.table.rowCount())
    def removeRow(self):
        selections = self.table.selectionModel().selectedRows()
        indices=[]
        for selection in selections:
            row=selection.row()
            if row not in indices:
                indices.append(row)
        for i in sorted(indices,reverse=True):
            self.table.removeRow(i)


#From https://stackoverflow.com/questions/26227885/drag-and-drop-rows-within-qtablewidget
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
            rows_to_move = [[QTableWidgetItem(self.item(row_index, column_index)) for column_index in range(self.columnCount())]
                            for row_index in rows]
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
        return rect.contains(pos, True) and not (int(self.model().flags(index)) & Qt.ItemIsDropEnabled) and pos.y() >= rect.center().y()

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
    os.system('cd {} && pip install -e .'.format(baseToolPath))


    
def editIPs(parent,serverSettings):
    serverDialog=ServerDialog(serverSettings)
    serverDialog.exec_()

