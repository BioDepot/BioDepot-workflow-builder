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


class OWPython3(OWBwBWidget):
    name = "Python3"
    description = "Minimum Python3 container with pip"
    priority = 10
    icon = getIconName(__file__, "python3.png")
    want_main_area = False
    docker_image_name = "biodepot/python"
    docker_image_tag = "3.10.9__alpine-3.17.1__8132a508__2132c733__29cfe2a0"
    inputs = [
        ("inputFile", str, "handleInputsinputFile"),
        ("Trigger", str, "handleInputsTrigger"),
    ]
    outputs = [("OutputDir", str)]
    pset = functools.partial(settings.Setting, schema_only=True)
    runMode = pset(0)
    exportGraphics = pset(False)
    runTriggers = pset([])
    triggerReady = pset({})
    inputConnectionsStore = pset({})
    optionsChecked = pset({})
    InputFile = pset(None)

    def __init__(self):
        super().__init__(self.docker_image_name, self.docker_image_tag)
        with open(getJsonName(__file__, "Python3")) as f:
            self.data = jsonpickle.decode(f.read())
            f.close()
        self.initVolumes()
        self.inputConnections = ConnectionDict(self.inputConnectionsStore)
        self.drawGUI()

    def handleInputsinputFile(self, value, *args):
        if args and len(args) > 0:
            self.handleInputs("inputFile", value, args[0][0], test=args[0][3])
        else:
            self.handleInputs("inputFile", value, None)

    def handleInputsTrigger(self, value, *args):
        if args and len(args) > 0:
            self.handleInputs("Trigger", value, args[0][0], test=args[0][3])
        else:
            self.handleInputs("inputFile", value, None)

    def handleOutputs(self):
        outputValue = None
        if hasattr(self, "OutputDir"):
            outputValue = getattr(self, "OutputDir")
        self.send("OutputDir", outputValue)
