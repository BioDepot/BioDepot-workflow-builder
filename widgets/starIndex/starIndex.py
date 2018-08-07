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

class OWStarindex(OWBwBWidget):
    name = "Star index"
    description = "Construct indices for STAR aligner "
    category = "RNA-seq"
    priority = 11
    icon = "/widgets/starIndex/icon/starIndex.png"
    want_main_area = False
    docker_image_name = "biodepot/star"
    docker_image_tag = "2.6.0c__debian-8.11-slim__072918"
    inputs = [("Trigger",str,"handleInputsTrigger")]
    outputs = [("genomeDir",str)]
    pset=functools.partial(settings.Setting,schema_only=True)
    runMode=pset(0)
    exportGraphics=pset(False)
    runTriggers=pset([])
    triggerReady=pset({})
    inputConnectionsStore=pset({})
    optionsChecked=pset({})
    rmode=pset("genomeGenerate")
    genomeDir=pset(None)
    genomeFastaFiles=pset([])
    genomeChrBinNbits=pset("18")
    genomeSAindexNbases=pset(14)
    genomeSAsparseD=pset(1)
    genomeSuffixLengthMax=pset(-1)
    runThreadN=pset(1)
    sjdbGTFfile=pset(None)
    sjdbFileChrStartEnd =pset([])
    sjdbGTFchrPrefix =pset("chr")
    sjdbGTFfeatureExon=pset("exon")
    sjdbGTFtagExonParentTranscript=pset("transcript_id")
    sjdbGTFtagExonParentGene=pset("gene_id")
    sjdbOverhang=pset(100)
    sjdbScore=pset(2)
    sjdbInsertSave =pset("Basic")
    def __init__(self):
        super().__init__(self.docker_image_name, self.docker_image_tag)
        with open("/widgets/starIndex/starIndex.json") as f:
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
    def handleOutputs(self):
        outputValue=None
        if hasattr(self,"genomeDir"):
            outputValue=getattr(self,"genomeDir")
        self.send("genomeDir", outputValue)
