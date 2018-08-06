import os
import glob
import sys
import functools
import jsonpickle
from collections import OrderedDict
from Orange.widgets import widget, gui, settings
import Orange.data
from Orange.data.io import FileFormat
from orangebiodepot.util.DockerClient import DockerClient
from orangebiodepot.util.BwBase import OWBwBWidget, ConnectionDict, BwbGuiElements
from PyQt5 import QtWidgets, QtGui

class OWkallistoindex(OWBwBWidget):
    name = "kallisto index"
    description = "Generates index files for kallisto"
    category = "RNA-seq"
    priority = 3
    icon = "/widgets/kallistoIndex/icon/kallistoindex.png"
    want_main_area = False
    docker_image_name = "biodepot/kallisto"
    docker_image_tag = "0.44.0__ubuntu-16.04__072818"
    inputs = [("trigger",str,"handleInputstrigger")]
    outputs = [("outputFilename",str)]
    pset=functools.partial(settings.Setting,schema_only=True)
    runMode=pset(0)
    exportGraphics=pset(False)
    runTriggers=pset([])
    triggerReady=pset({})
    inputConnectionsStore=pset({})
    optionsChecked=pset({})
    fastaFiles=pset([])
    outputFilename=pset(None)
    kmerSize=pset(31)
    makeUnique=pset(False)
    def __init__(self):
        super().__init__(self.docker_image_name, self.docker_image_tag)
        with open("/widgets/kallistoIndex/kallistoIndex.json") as f:
            self.data=jsonpickle.decode(f.read())
            f.close()
        self.initVolumes()
        self.inputConnections = ConnectionDict(self.inputConnectionsStore)
        self.drawGUI()
    def handleInputstrigger(self, value, *args):
        if args and len(args) > 0: 
            self.handleInputs("trigger", value, args[0][0])
        else:
            self.handleInputs("inputFile", value, None)
    def handleOutputs(self):
        outputValue=None
        if hasattr(self,"outputFilename"):
            outputValue=getattr(self,"outputFilename")
        self.send("outputFilename", outputValue)
