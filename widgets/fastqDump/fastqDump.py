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
from BwBase import OWBwBWidget, ConnectionDict, BwbGuiElements
from PyQt5 import QtWidgets, QtGui

class OWfastqDump(OWBwBWidget):
    name = "fastqDump"
    description = "Download fastq files from GEO"
    category = "Utilities"
    priority = 2
    icon = "/widgets/fastqDump/icon/pfqDump.png"
    want_main_area = False
    docker_image_name = "biodepot/sratoolkit"
    docker_image_tag = "2.8.2-1__minideb-jessie__072818"
    inputs = [("Trigger",str,"handleInputsTrigger")]
    outputs = [("OutputDir",str)]
    pset=functools.partial(settings.Setting,schema_only=True)
    runMode=pset(0)
    exportGraphics=pset(False)
    runTriggers=pset([])
    triggerReady=pset({})
    inputConnectionsStore=pset({})
    optionsChecked=pset({})
    accession=pset(None)
    OutputDir=pset("/data")
    fqdversion=pset(False)
    help=pset(False)
    sraids=pset([])
    tableName=pset(None)
    splitSpot=pset(False)
    minSpotID=pset(None)
    maxSpotID=pset(None)
    spotGroups=pset([])
    clip=pset(False)
    minReadLen=pset(None)
    readFilter=pset(None)
    qualFilter=pset(None)
    qualFilter1=pset(None)
    aligned=pset(False)
    unaligned=pset(False)
    alignedRegion=pset(None)
    matePairDistance=pset(None)
    useStdout=pset(False)
    gzip=pset(False)
    bzip2=pset(False)
    splitFiles=pset(False)
    split3=pset(False)
    spotGroup=pset(False)
    groupinDirs=pset(False)
    keepEmpty=pset(False)
    dumpcs=pset(None)
    qOffset=pset(33)
    fasta=pset(False)
    suppressQual=pset(False)
    origfmt=pset(False)
    readids=pset(False)
    helicos=pset(False)
    deflineSeq=pset(None)
    deflineQual=pset(None)
    disablemt=pset(False)
    logLevel=pset(None)
    verbose=pset(False)
    ncbiError=pset(False)
    legacyError=pset(False)
    def __init__(self):
        super().__init__(self.docker_image_name, self.docker_image_tag)
        with open("/widgets/fastqDump/fastqDump.json") as f:
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
        outputValue="/data"
        if hasattr(self,"OutputDir"):
            outputValue=getattr(self,"OutputDir")
        self.send("OutputDir", outputValue)
