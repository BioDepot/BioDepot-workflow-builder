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
    icon = "/widgets/deseq2/icon/deseq2.png"
    want_main_area = False
    docker_image_name = "biodepot/deseq2"
    docker_image_tag = "1.20__ubuntu-16.04__bioc-3.7__r-3.5.1__072918"
    inputs = [("Trigger",str,"handleInputsTrigger"),("countsFile",str,"handleInputscountsFile")]
    outputs = [("topGenesFile",str)]
    pset=functools.partial(settings.Setting,schema_only=True)
    runMode=pset(0)
    exportGraphics=pset(False)
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
        with open("/widgets/deseq2/deseq2.json") as f:
            self.data=jsonpickle.decode(f.read())
            f.close()
        self.initVolumes()
        self.inputConnections = ConnectionDict(self.inputConnectionsStore)
        self.drawGUI()
    def handleInputsTrigger(self, value, *args):
        if args and len(args) > 0: 
            self.handleInputs("Trigger", value, args[0][0])
        else:
            self.handleInputs("inputFile", value, None)
    def handleInputscountsFile(self, value, *args):
        if args and len(args) > 0: 
            self.handleInputs("countsFile", value, args[0][0])
        else:
            self.handleInputs("inputFile", value, None)
    def handleOutputs(self):
        topGenesFile="top{}Genes.tsv".format(self.ngenes)
        sys.stderr.write("Top file is {}\n".format(topGenesFile))
        if self.outputDir is not None:
            topGenesFile=self.outputDir+"/"+topGenesFile
            self.send('topGenesFile',topGenesFile)
