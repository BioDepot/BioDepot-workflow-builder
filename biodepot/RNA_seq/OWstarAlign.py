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

class OWStaralign(OWBwBWidget):
    name = "Star align"
    description = "STAR aligner alignment module"
    category = "RNA-seq"
    priority = 10
    icon = "/biodepot/RNA_seq/icons/staralign.png"
    want_main_area = False
    docker_image_name = "biodepot/alpine-star"
    docker_image_tag = "3.7-020201"
    inputs = [("Trigger",str,"handleInputsTrigger"),("GenomeDir",str,"handleInputsGenomeDir"),("fastqDir",str,"handleInputsfastqDir")]
    outputs = [("outputDir",str)]
    pset=functools.partial(settings.Setting,schema_only=True)
    runMode=pset(0)
    runTriggers=pset([])
    triggerReady=pset({})
    inputConnectionsStore=pset({})
    optionsChecked=pset({})
    fastqDir=pset(None)
    GenomeDir=pset(None)
    def __init__(self):
        super().__init__(self.docker_image_name, self.docker_image_tag)
        with open("/biodepot/RNA_seq/json/starAlign.json") as f:
            self.data=jsonpickle.decode(f.read())
            f.close()
        self.initVolumes()
        self.inputConnections = ConnectionDict(self.inputConnectionsStore)
        self.drawGUI()
    def handleInputsTrigger(self, value, sourceId=None):
        self.handleInputs(value, "Trigger", sourceId=None)
    def handleInputsGenomeDir(self, value, sourceId=None):
        self.handleInputs(value, "GenomeDir", sourceId=None)
    def handleInputsfastqDir(self, value, sourceId=None):
        self.handleInputs(value, "fastqDir", sourceId=None)
    def handleOutputs(self):
        outputValue=None
        if hasattr(self,"outputDir"):
            outputValue=getattr(self,"outputDir")
        self.send("outputDir", outputValue)
