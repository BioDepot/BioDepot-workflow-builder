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
from BwBase import OWBwBWidget, ConnectionDict, BwbGuiElements
from PyQt5 import QtWidgets, QtGui

class OWfastqc(OWBwBWidget):
    name = "fastqc"
    description = "fastqc"
    category = "Utilities"
    priority = 5
    icon = "/widgets/fastqc/icon/fastqc_icon_100.png"
    want_main_area = False
    docker_image_name = "biodepot/fastqc"
    docker_image_tag = "latest"
    inputs = [("inputDir",str,"handleInputsinputDir")]
    outputs = [("outputDir",str)]
    pset=functools.partial(settings.Setting,schema_only=True)
    runMode=pset(0)
    exportGraphics=pset(False)
    runTriggers=pset([])
    triggerReady=pset({})
    inputConnectionsStore=pset({})
    optionsChecked=pset({})
    inputFiles=pset([])
    outputDir=pset(None)
    inputDir=pset(None)
    def __init__(self):
        super().__init__(self.docker_image_name, self.docker_image_tag)
        with open("/widgets/fastqc/fastqc.json") as f:
            self.data=jsonpickle.decode(f.read())
            f.close()
        self.initVolumes()
        self.inputConnections = ConnectionDict(self.inputConnectionsStore)
        self.drawGUI()
    def handleInputsinputDir(self, value, *args):
        if args and len(args) > 0: 
            self.handleInputs("inputDir", value, args[0][0])
        else:
            self.handleInputs("inputFile", value, None)
    def handleOutputs(self):
        outputValue=None
        if hasattr(self,"outputDir"):
            outputValue=getattr(self,"outputDir")
        self.send("outputDir", outputValue)
