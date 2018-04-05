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

class OWkallisto_quant(OWBwBWidget):
    name = "kallisto_quant"
    description = "Alignment and quantification of reads from fastq files"
    category = "RNA-seq"
    priority = 10
    icon = "/biodepot/orangebiodepot/icons/kallisto-analysis.svg"
    want_main_area = False
    docker_image_name = "biodepot/kallisto"
    docker_image_tag = "0.44"
    inputs = [("indexFile",str,"handleInputsindexFile"),("fastqFiles",str,"handleInputsfastqFiles")]
    outputs = [("outputDir",str)]
    pset=functools.partial(settings.Setting,schema_only=True)
    runMode=pset(0)
    runTriggers=pset([])
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
    def __init__(self):
        super().__init__(self.docker_image_name, self.docker_image_tag)
        with open("/biodepot/orangebiodepot/json/kallisto.json") as f:
            self.data=jsonpickle.decode(f.read())
            f.close()
        self.initVolumes()
        self.inputConnections = ConnectionDict(self.inputConnectionsStore)
        self.drawGUI()
    def handleInputsindexFile(self, value, sourceId=None):
        self.handleInputs(value, "indexFile", sourceId=None)
    def handleInputsfastqFiles(self, value, sourceId=None):
        self.handleInputs(value, "fastqFiles", sourceId=None)
    def handleOutputs(self):
        outputValue=None
        if hasattr(self,"outputDir"):
            outputValue=getattr(self,"outputDir")
        self.send("outputDir", outputValue)
