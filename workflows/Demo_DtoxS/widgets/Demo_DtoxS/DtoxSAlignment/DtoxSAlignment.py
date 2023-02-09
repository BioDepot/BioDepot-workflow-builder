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

class OWDtoxSAlignment(OWBwBWidget):
    name = "DtoxSAlignment"
    description = "Alignment part of DtoxS standard operating procedure (SOP)"
    priority = 1
    icon = getIconName(__file__,"dtoxs-alignment2.svg")
    want_main_area = False
    docker_image_name = "biodepot/dtoxsalignment"
    docker_image_tag = "bwa-0.7.17__buster-slim__441c5219__9b21142b__a20f8157"
    inputs = [("barcodesFile",str,"handleInputsbarcodesFile"),("topDir",str,"handleInputstopDir"),("trigger",str,"handleInputstrigger")]
    outputs = [("topDir",str)]
    pset=functools.partial(settings.Setting,schema_only=True)
    runMode=pset(0)
    exportGraphics=pset(False)
    runTriggers=pset([])
    triggerReady=pset({})
    inputConnectionsStore=pset({})
    optionsChecked=pset({})
    topDir=pset(None)
    barcodesFile=pset("barcodes_trugrade_96_set4.dat")
    lanes=pset(6)
    seriesName=pset("20150409")
    def __init__(self):
        super().__init__(self.docker_image_name, self.docker_image_tag)
        with open(getJsonName(__file__,"DtoxSAlignment")) as f:
            self.data=jsonpickle.decode(f.read())
            f.close()
        self.initVolumes()
        self.inputConnections = ConnectionDict(self.inputConnectionsStore)
        self.drawGUI()
    def handleInputsbarcodesFile(self, value, *args):
        if args and len(args) > 0: 
            self.handleInputs("barcodesFile", value, args[0][0], test=args[0][3])
        else:
            self.handleInputs("inputFile", value, None, False)
    def handleInputstopDir(self, value, *args):
        if args and len(args) > 0: 
            self.handleInputs("topDir", value, args[0][0], test=args[0][3])
        else:
            self.handleInputs("inputFile", value, None, False)
    def handleInputstrigger(self, value, *args):
        if args and len(args) > 0: 
            self.handleInputs("trigger", value, args[0][0], test=args[0][3])
        else:
            self.handleInputs("inputFile", value, None, False)
    def handleOutputs(self):
        outputValue=None
        if hasattr(self,"topDir"):
            outputValue=getattr(self,"topDir")
        self.send("topDir", outputValue)
