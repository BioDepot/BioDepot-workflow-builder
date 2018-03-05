import os
import glob
import sys
import functools
import jsonpickle
from collections import OrderedDict
from Orange.widgets import widget, gui, settings
from orangebiodepot.util.DockerClient import DockerClient
from orangebiodepot.util.BwBase import OWBwBWidget, ConnectionDict, ContainerPaths, BwbGuiElements
from PyQt5 import QtWidgets, QtGui

class OWdownloadURL(OWBwBWidget):
    name = "downloadURL"
    description = "Downloads files from URL"
    category = "Utility"
    priority = 10
    icon = "/biodepot/orangebiodepot/icons/download.png"
    want_main_area = False
    docker_image_name = "biodepot/alpine-download"
    docker_image_tag = "1.0"
    inputs = [("directory",str,"handleInputsdirectory"),("URL",str,"handleInputsURL"),("trigger",str,"handleInputstrigger")]
    outputs = [("directory",str)]
    pset=functools.partial(settings.Setting,schema_only=True)
    runMode=pset(0)
    runTriggers=pset([])
    inputConnectionsStore=pset({})
    optionsChecked=pset({})
    directory=pset(None)
    URL=pset(None)
    decompress=pset(True)
    def __init__(self):
        super().__init__(self.docker_image_name, self.docker_image_tag)
        with open("/biodepot/orangebiodepot/json/downloadURL.json") as f:
            self.data=jsonpickle.decode(f.read())
            f.close()
        self.initVolumes()
        self.inputConnections = ConnectionDict(self.inputConnectionsStore)
        self.drawGUI()
    def handleInputsdirectory(self, value, sourceId=None):
        self.handleInputs(value, "directory", sourceId=None)
    def handleInputsURL(self, value, sourceId=None):
        self.handleInputs(value, "URL", sourceId=None)
    def handleInputstrigger(self, value, sourceId=None):
        self.handleInputs(value, "trigger", sourceId=None)
    def handleOutputs(self):
        outputValue=None
        if hasattr(self,"directory"):
            outputValue=getattr(self,"directory")
        self.send("directory", outputValue)
