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
    priority = 10
    icon = "/biodepot/RNA_seq/kallistoIndex/kallistoindex.png"
    want_main_area = False
    docker_image_name = "biodepot/kallisto"
    docker_image_tag = "0.44.0__ubuntu-16.04__072818"
    inputs = [("trigger",str,"handleInputstrigger")]
    outputs = [("outputFilename",str)]
    pset=functools.partial(settings.Setting,schema_only=True)
    runMode=pset(0)
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
        with open("/biodepot/RNA_seq/json/kallistoIndex.json") as f:
            self.data=jsonpickle.decode(f.read())
            f.close()
        self.initVolumes()
        self.inputConnections = ConnectionDict(self.inputConnectionsStore)
        self.drawGUI()
    def handleInputstrigger(self, value, sourceId=None):
        self.handleInputs(value, "trigger", sourceId=None)
    def handleOutputs(self):
        outputValue=None
        if hasattr(self,"outputFilename"):
            outputValue=getattr(self,"outputFilename")
        self.send("outputFilename", outputValue)
