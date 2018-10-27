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
from PyQt5.QtWidgets import QInputDialog, QLineEdit, QMainWindow, QApplication, QPushButton, QWidget, QAction, QTabWidget, QVBoxLayout, QMessageBox, QFileDialog


from AnyQt.QtWidgets import (
    QWidget, QButtonGroup, QGroupBox, QRadioButton, QSlider,
    QDoubleSpinBox, QComboBox, QSpinBox, QListView, QLabel,
    QScrollArea, QVBoxLayout, QHBoxLayout, QFormLayout,
    QSizePolicy, QApplication, QCheckBox
)

defaultIconFile='/icons/default.png'


def registerDirectory(baseToolPath):
    os.system('cd {} && pip install -e .'.format(baseToolPath))

def parseIPFile(serverSettings,ipFile):
    if 'addrs' not in serverSettings:
        serverSettings['addrs']={}
    with open (ipFile,'r') as f:
        for line in f:
            parts=line.split()
            if parts and parts[0]:
                addr=parts.pop(0)
                try:
                    #store dict with ip and volumeMappings as keys
                    if addr not in serverSettings['addrs']:
                            serverSettings['addrs'][addr]={}
                    if parts:
                        #form attr=value
                        settings=serverSettings['addrs'][addr]
                        for attrValue in parts:
                            attr,value=attrValue.split('=')
                            if not attr in settings:
                                settings[attr]=[value]
                            elif not value in settings[attr]:
                                settings[attr].append(value)
                except socket.error:
                    QtGui.QMessageBox.warning('IP error','','Rejected {} - Not a valid IP'.format(addr))
    return
    
def importIPs(parent,serverSettings):
    startDir='/'
    if 'saveIPDir' in serverSettings and serverSettings['saveIPDir']:
        startDir=serverSettings['saveIPDir']
    elif os.path.isdir('/data'):
        startDir='/data'

    ipFile = QFileDialog.getOpenFileName(parent,"Load IPs from File",startDir)[0]
    if not ipFile:
        return QFileDialog.Rejected
    serverSettings['saveIPDir']=os.path.dirname(ipFile)
    parseIPFile(serverSettings,ipFile)
    return QFileDialog.Accepted
