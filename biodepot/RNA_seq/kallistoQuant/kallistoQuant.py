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

class OWkallistoquant(OWBwBWidget):
    name = "kallisto quant"
    description = "Alignment and quantification of reads from fastq files"
    category = "RNA-seq"
    priority = 10
    icon = "/biodepot/RNA_seq/kallistoQuant/icon/kallistoquant.png"
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
        with open("/biodepot/RNA_seq/json/kallistoQuant.json") as f:
            self.data=jsonpickle.decode(f.read())
            f.close()
        self.initVolumes()
        self.inputConnections = ConnectionDict(self.inputConnectionsStore)
        self.drawGUI()
    def handleInputsindexFile(self, value, sourceId=None):
        self.handleInputs(value, "indexFile", sourceId=None)
    def handleInputsfastqFiles(self, value, sourceId=None):
        self.handleInputs(value, "fastqFiles", sourceId=None)
    def handleInputstrigger(self, value, sourceId=None):
        self.handleInputs(value, "trigger", sourceId=None)
    def handleOutputs(self):
        outputValue=None
        if hasattr(self,"outputDir"):
            outputValue=getattr(self,"outputDir")
        self.send("outputDir", outputValue)
