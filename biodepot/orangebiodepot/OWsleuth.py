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

class OWsleuth(OWBwBWidget):
    name = "sleuth"
    description = "Analysis of kallisto quantified RNA-seq files"
    category = "RNA-seq"
    priority = 10
    icon = "/biodepot/orangebiodepot/icons/sleuth2.png"
    want_main_area = False
    docker_image_name = "biodepot/ubuntu-sleuth"
    docker_image_tag = "1.0"
    inputs = [("kallistoDirectory",str,"handleInputskallistoDirectory"),("trigger",str,"handleInputstrigger")]
    outputs = [("OutputDirectory",str),("topGenes",Orange.data.Table)]
    pset=functools.partial(settings.Setting,schema_only=True)
    runMode=pset(0)
    runTriggers=pset([])
    triggerReady=pset({})
    inputConnectionsStore=pset({})
    optionsChecked=pset({})
    kallistoDirectory=pset(None)
    full_model=pset(None)
    def __init__(self):
        super().__init__(self.docker_image_name, self.docker_image_tag)
        with open("/biodepot/orangebiodepot/json/sleuth.json") as f:
            self.data=jsonpickle.decode(f.read())
            f.close()
        self.initVolumes()
        self.inputConnections = ConnectionDict(self.inputConnectionsStore)
        self.drawGUI()
    def handleInputskallistoDirectory(self, value, sourceId=None):
        self.handleInputs(value, "kallistoDirectory", sourceId=None)
    def handleInputstrigger(self, value, sourceId=None):
        self.handleInputs(value, "trigger", sourceId=None)
    def handleOutputs(self):
        outputValue=None
        if hasattr(self,"OutputDirectory"):
            outputValue=getattr(self,"OutputDirectory")
        self.send("OutputDirectory", outputValue)
        outputValue=None
        if hasattr(self,"topGenes"):
            outputValue=getattr(self,"topGenes")
        self.send("topGenes", outputValue)
