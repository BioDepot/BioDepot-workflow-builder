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

class OWkallistoQuant(OWBwBWidget):
    name = "kallistoQuant"
    description = "Alignment and quantification of reads from fastq files"
    priority = 4
    icon = getIconName(__file__,"kallistoquant.png")
    want_main_area = False
    docker_image_name = "biodepot/kallisto"
    docker_image_tag = "0.44.0__ubuntu-16.04__072818"
    inputs = [("indexFile",str,"handleInputsindexFile"),("fastqFiles",str,"handleInputsfastqFiles"),("trigger",str,"handleInputstrigger")]
    outputs = [("outputDir",str)]
    pset=functools.partial(settings.Setting,schema_only=True)
    runMode=pset(0)
    exportGraphics=pset(False)
    runTriggers=pset([])
    triggerReady=pset({})
    inputConnectionsStore=pset({})
    optionsChecked=pset({})
    outputDir=pset(None)
    indexFile=pset(None)
    fastqFiles=pset([])
    bootstrap=pset(30)
    seed=pset(42)
    plaintext=pset(False)
    fusion=pset(False)
    single=pset(False)
    single_overhang=pset(False)
    fr_stranded=pset(False)
    rf_stranded=pset(False)
    fragment_length=pset(None)
    stdev=pset(None)
    nThreads=pset(1)
    pseudoBam=pset(False)
    genomeBam=pset(False)
    gtf=pset(None)
    multiSample=pset(1)
    def __init__(self):
        super().__init__(self.docker_image_name, self.docker_image_tag)
        with open(getJsonName(__file__,"kallistoQuant")) as f:
            self.data=jsonpickle.decode(f.read())
            f.close()
        self.initVolumes()
        self.inputConnections = ConnectionDict(self.inputConnectionsStore)
        self.drawGUI()
    def handleInputsindexFile(self, value, *args):
        if args and len(args) > 0: 
            self.handleInputs("indexFile", value, args[0][0]), test=args[0][3]))
        else:
            self.handleInputs("inputFile", value, None)
    def handleInputsfastqFiles(self, value, *args):
        if args and len(args) > 0: 
            self.handleInputs("fastqFiles", value, args[0][0]), test=args[0][3]))
        else:
            self.handleInputs("inputFile", value, None)
    def handleInputstrigger(self, value, *args):
        if args and len(args) > 0: 
            self.handleInputs("trigger", value, args[0][0]), test=args[0][3]))
        else:
            self.handleInputs("inputFile", value, None)
    def handleOutputs(self):
        outputValue=None
        if hasattr(self,"outputDir"):
            outputValue=getattr(self,"outputDir")
        self.send("outputDir", outputValue)
