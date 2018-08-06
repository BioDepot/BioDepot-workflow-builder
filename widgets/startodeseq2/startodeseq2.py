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
    description = "Convert Star quantMode counts file to DESeq2 style counts file"
    category = "RNA-seq"
    priority = 10
    icon = "/widgets/startodeseq2/icon/startodeseq2.png"
    want_main_area = False
    docker_image_name = "biodepot/star2deseq"
    docker_image_tag = "1.0__alpine-3.7__07-29-18"
    inputs = [("inputDirs",str,"handleInputsinputDirs"),("Trigger",str,"handleInputsTrigger")]
    outputs = [("outputFile",str)]
    pset=functools.partial(settings.Setting,schema_only=True)
    runMode=pset(0)
    exportGraphics=pset(False)
    runTriggers=pset([])
    triggerReady=pset({})
    inputConnectionsStore=pset({})
    optionsChecked=pset({})
    outputFile=pset(None)
    inputFile=pset("ReadsPerGene.out.tab")
    column=pset(4)
    inputDirs=pset([])
    def __init__(self):
        super().__init__(self.docker_image_name, self.docker_image_tag)
        with open("/widgets/startodeseq2/startodeseq2.json") as f:
            self.data=jsonpickle.decode(f.read())
            f.close()
        self.initVolumes()
        self.inputConnections = ConnectionDict(self.inputConnectionsStore)
        self.drawGUI()
    def handleInputsinputDirs(self, value, *args):
        if args and len(args) > 0: 
            self.handleInputs("inputDirs", value, args[0][0])
        else:
            self.handleInputs("inputFile", value, None)
    def handleInputsTrigger(self, value, *args):
        if args and len(args) > 0: 
            self.handleInputs("Trigger", value, args[0][0])
        else:
            self.handleInputs("inputFile", value, None)
    def handleOutputs(self):
        outputValue=None
        if hasattr(self,"outputFile"):
            outputValue=getattr(self,"outputFile")
        self.send("outputFile", outputValue)
