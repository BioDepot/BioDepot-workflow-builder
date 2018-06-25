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

class OWStartoDESeq2(OWBwBWidget):
    name = "Star to DESeq2"
    description = "Convert Star quantMode to DESeq2 style counts"
    category = "RNA-seq"
    priority = 10
    icon = "/biodepot/RNA_seq/icons/startodeseq2.png"
    want_main_area = False
    docker_image_name = "alpine-star-deseq2"
    docker_image_tag = "3.7-1.0"
    inputs = [("outputFile",str,"handleInputsoutputFile"),("inputFile",str,"handleInputsinputFile")]
    outputs = [("outputFile",str)]
    pset=functools.partial(settings.Setting,schema_only=True)
    runMode=pset(0)
    runTriggers=pset([])
    triggerReady=pset({})
    inputConnectionsStore=pset({})
    optionsChecked=pset({})
    outputFile=pset(None)
    Column=pset(4)
    def __init__(self):
        super().__init__(self.docker_image_name, self.docker_image_tag)
        with open("/biodepot/RNA_seq/json/startodeseq2.json") as f:
            self.data=jsonpickle.decode(f.read())
            f.close()
        self.initVolumes()
        self.inputConnections = ConnectionDict(self.inputConnectionsStore)
        self.drawGUI()
    def handleInputsoutputFile(self, value, sourceId=None):
        self.handleInputs(value, "outputFile", sourceId=None)
    def handleInputsinputFile(self, value, sourceId=None):
        self.handleInputs(value, "inputFile", sourceId=None)
    def handleOutputs(self):
        outputValue=None
        if hasattr(self,"outputFile"):
            outputValue=getattr(self,"outputFile")
        self.send("outputFile", outputValue)
