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


class OWdeseq2(OWBwBWidget):
    name = "deseq2"
    description = "Differential expression using DESeq2 package"
    priority = 14
    icon = getIconName(__file__, "deseq2.png")
    want_main_area = False
    docker_image_name = "biodepot/deseq2"
    docker_image_tag = "bioc-r_3.16-ubuntu-22.04-r-4.2.2__375472a0__186da45b__add50ffe"
    inputs = [
        ("Trigger", str, "handleInputsTrigger"),
        ("countsFile", str, "handleInputscountsFile"),
    ]
    outputs = [("topGenesFile", str)]
    pset = functools.partial(settings.Setting, schema_only=True)
    runMode = pset(0)
    exportGraphics = pset(False)
    runTriggers = pset([])
    triggerReady = pset({})
    inputConnectionsStore = pset({})
    optionsChecked = pset({})
    countsFile = pset(None)
    descriptionFile = pset(None)
    control = pset(None)
    treatment = pset(None)
    outputDir = pset(None)
    ngenes = pset(50)
    topGenesFile = pset(None)

    def __init__(self):
        super().__init__(self.docker_image_name, self.docker_image_tag)
        with open(getJsonName(__file__, "deseq2")) as f:
            self.data = jsonpickle.decode(f.read())
            f.close()
        self.initVolumes()
        self.inputConnections = ConnectionDict(self.inputConnectionsStore)
        self.drawGUI()

    def handleInputsTrigger(self, value, *args):
        if args and len(args) > 0:
            self.handleInputs("Trigger", value, args[0][0], test=args[0][3])
        else:
            self.handleInputs("inputFile", value, None)

    def handleInputscountsFile(self, value, *args):
        if args and len(args) > 0:
            self.handleInputs("countsFile", value, args[0][0], test=args[0][3])
        else:
            self.handleInputs("inputFile", value, None)

    def handleOutputs(self):
        outputValue = None
        if hasattr(self, "topGenesFile"):
            outputValue = getattr(self, "topGenesFile")
        self.send("topGenesFile", outputValue)
