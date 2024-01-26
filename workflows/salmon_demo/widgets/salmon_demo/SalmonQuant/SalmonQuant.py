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

class OWSalmonQuant(OWBwBWidget):
    name = "SalmonQuant"
    description = "Alignment and quantification of reads from fastq files"
    priority = 4
    icon = getIconName(__file__,"salmon-feat.png")
    want_main_area = False
    docker_image_name = "biodepot/salmon"
    docker_image_tag = "latest"
    inputs = [("index",str,"handleInputsindex"),("trigger",str,"handleInputstrigger"),("mates1",str,"handleInputsmates1"),("mates2",str,"handleInputsmates2"),("unmatedReads",str,"handleInputsunmatedReads"),("outputDirs",str,"handleInputsoutputDirs")]
    outputs = [("trigger",str)]
    pset=functools.partial(settings.Setting,schema_only=True)
    runMode=pset(0)
    exportGraphics=pset(False)
    runTriggers=pset([])
    triggerReady=pset({})
    inputConnectionsStore=pset({})
    optionsChecked=pset({})
    index=pset(None)
    libType=pset("A")
    mates1=pset([])
    mates2=pset([])
    unmatedReads=pset([])
    outputDirs=pset([])
    nThreads=pset(2)
    geneMap=pset(None)
    def __init__(self):
        super().__init__(self.docker_image_name, self.docker_image_tag)
        with open(getJsonName(__file__,"SalmonQuant")) as f:
            self.data=jsonpickle.decode(f.read())
            f.close()
        self.initVolumes()
        self.inputConnections = ConnectionDict(self.inputConnectionsStore)
        self.drawGUI()
    def handleInputsindex(self, value, *args):
        if args and len(args) > 0: 
            self.handleInputs("index", value, args[0][0], test=args[0][3])
        else:
            self.handleInputs("inputFile", value, None, False)
    def handleInputstrigger(self, value, *args):
        if args and len(args) > 0: 
            self.handleInputs("trigger", value, args[0][0], test=args[0][3])
        else:
            self.handleInputs("inputFile", value, None, False)
    def handleInputsmates1(self, value, *args):
        if args and len(args) > 0: 
            self.handleInputs("mates1", value, args[0][0], test=args[0][3])
        else:
            self.handleInputs("inputFile", value, None, False)
    def handleInputsmates2(self, value, *args):
        if args and len(args) > 0: 
            self.handleInputs("mates2", value, args[0][0], test=args[0][3])
        else:
            self.handleInputs("inputFile", value, None, False)
    def handleInputsunmatedReads(self, value, *args):
        if args and len(args) > 0: 
            self.handleInputs("unmatedReads", value, args[0][0], test=args[0][3])
        else:
            self.handleInputs("inputFile", value, None, False)
    def handleInputsoutputDirs(self, value, *args):
        if args and len(args) > 0: 
            self.handleInputs("outputDirs", value, args[0][0], test=args[0][3])
        else:
            self.handleInputs("inputFile", value, None, False)
    def handleOutputs(self):
        outputValue=None
        if hasattr(self,"trigger"):
            outputValue=getattr(self,"trigger")
        self.send("trigger", outputValue)
