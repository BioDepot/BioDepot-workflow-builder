import os
import glob
import sys
import functools
import jsonpickle
from collections import OrderedDict
from Orange.widgets import widget, gui, settings
import Orange.data
from Orange.data.io import FileFormat
from DockerClient import DockerClient
from BwBase import OWBwBWidget, ConnectionDict, BwbGuiElements, getIconName, getJsonName
from PyQt5 import QtWidgets, QtGui

class OWDtoxSAnalysis(OWBwBWidget):
    name = "DtoxSAnalysis"
    description = "Step 2 of Dtoxs SOP. Uses edgeR for differential expression analysis"
    category = "RNA-seq"
    priority = 2
    icon = getIconName(__file__,"dtoxs-analysis2.svg")
    want_main_area = False
    docker_image_name = "biodepot/dtoxs_analysis"
    docker_image_tag = "1.0__ubuntu-16.04__bioc-3.6__r-3.4.3__072818"
    inputs = [("RepositoryDirectory",str,"handleInputsRepositoryDirectory"),("ConfigurationFile",str,"handleInputsConfigurationFile"),("Trigger",str,"handleInputsTrigger")]
    outputs = [("ResultsDirectory",str),("topGenesFile",str)]
    pset=functools.partial(settings.Setting,schema_only=True)
    runMode=pset(0)
    exportGraphics=pset(False)
    runTriggers=pset([])
    triggerReady=pset({})
    inputConnectionsStore=pset({})
    optionsChecked=pset({})
    RepositoryDirectory=pset(None)
    ConfigurationFile=pset(None)
    def __init__(self):
        super().__init__(self.docker_image_name, self.docker_image_tag)
        with open(getJsonName(__file__,"DtoxSAnalysis")) as f:
            self.data=jsonpickle.decode(f.read())
            f.close()
        self.initVolumes()
        self.inputConnections = ConnectionDict(self.inputConnectionsStore)
        self.drawGUI()
    def handleInputsRepositoryDirectory(self, value, *args):
        if args and len(args) > 0: 
            self.handleInputs("RepositoryDirectory", value, args[0][0])
        else:
            self.handleInputs("inputFile", value, None)
    def handleInputsConfigurationFile(self, value, *args):
        if args and len(args) > 0: 
            self.handleInputs("ConfigurationFile", value, args[0][0])
        else:
            self.handleInputs("inputFile", value, None)
    def handleInputsTrigger(self, value, *args):
        if args and len(args) > 0: 
            self.handleInputs("Trigger", value, args[0][0])
        else:
            self.handleInputs("inputFile", value, None)
    def handleOutputs(self):
        outputValue="/data"
        if hasattr(self,"ResultsDirectory"):
            outputValue=getattr(self,"ResultsDirectory")
        self.send("ResultsDirectory", outputValue)
        outputValue=None
        if hasattr(self,"topGenesFile"):
            outputValue=getattr(self,"topGenesFile")
        self.send("topGenesFile", outputValue)
