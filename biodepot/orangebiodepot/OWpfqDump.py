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

class OWpfastqDump(OWBwBWidget):
    name = "pfastqDump"
    description = "Download fastq files from GEO"
    category = "Utilities"
    priority = 10
    icon = "/biodepot/orangebiodepot/icons/pfqDump.png"
    want_main_area = False
    docker_image_name = "biodepot/pfastq-dump"
    docker_image_tag = "1.0"
    inputs = [("Trigger",str,"handleInputsTrigger")]
    outputs = [("OutputDir",str)]
    pset=functools.partial(settings.Setting,schema_only=True)
    runMode=pset(0)
    runTriggers=pset([])
    inputConnectionsStore=pset({})
    optionsChecked=pset({})
    IDs=pset([])
    nthreads=pset(1)
    OutputDir=pset("/data")
    tempdir=pset("/data")
    version=pset(False)
    help=pset(False)
    def __init__(self):
        super().__init__(self.docker_image_name, self.docker_image_tag)
        with open("/biodepot/orangebiodepot/json/pfqDump.json") as f:
            self.data=jsonpickle.decode(f.read())
            f.close()
        self.initVolumes()
        self.inputConnections = ConnectionDict(self.inputConnectionsStore)
        self.drawGUI()
    def handleInputsTrigger(self, value, sourceId=None):
        self.handleInputs(value, "Trigger", sourceId=None)
    def handleOutputs(self):
        outputValue="/data"
        if hasattr(self,"OutputDir"):
            outputValue=getattr(self,"OutputDir")
        self.send("OutputDir", outputValue)
