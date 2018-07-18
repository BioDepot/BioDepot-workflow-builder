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

class OWDESeq2(OWBwBWidget):
    name = "DESeq2"
    description = "Differential expression using DESeq2 package"
    category = "RNA-seq"
    priority = 10
    icon = "/biodepot/RNA_seq/icons/deseq2.png"
    want_main_area = False
    docker_image_name = "biodepot/deseq2"
    docker_image_tag = "1.20-ubuntu-16.04-bioc-3.7-r-3.5.1"
    inputs = [("Trigger",str,"handleInputsTrigger")]
    outputs = [("topGenes",Orange.data.Table)]
    pset=functools.partial(settings.Setting,schema_only=True)
    runMode=pset(0)
    runTriggers=pset([])
    triggerReady=pset({})
    inputConnectionsStore=pset({})
    optionsChecked=pset({})
    countsFile=pset(None)
    descriptionFile=pset(None)
    control=pset(None)
    treatment=pset(None)
    outputDir=pset(None)
    ngenes=pset(50)
    def __init__(self):
        super().__init__(self.docker_image_name, self.docker_image_tag)
        with open("/biodepot/RNA_seq/json/deseq2.json") as f:
            self.data=jsonpickle.decode(f.read())
            f.close()
        self.initVolumes()
        self.inputConnections = ConnectionDict(self.inputConnectionsStore)
        self.drawGUI()
    def handleInputsTrigger(self, value, sourceId=None):
        self.handleInputs(value, "Trigger", sourceId=None)
    def handleOutputs(self):
        outputValue=None
        if hasattr(self,"topGenes"):
            outputValue=getattr(self,"topGenes")
        self.send("topGenes", outputValue)
