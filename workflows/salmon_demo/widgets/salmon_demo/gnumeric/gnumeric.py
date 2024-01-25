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


class OWgnumeric(OWBwBWidget):
    name = "gnumeric"
    description = "Open source spreadsheet"
    priority = 2
    icon = getIconName(__file__, "gnumeric.png")
    want_main_area = False
    docker_image_name = "biodepot/gnumeric"
    docker_image_tag = "1.12.44__buster-slim__e4105988__f8d184ab__7b51306f"
    inputs = [
        ("inputFile", str, "handleInputsinputFile"),
        ("Trigger", str, "handleInputsTrigger"),
    ]
    pset = functools.partial(settings.Setting, schema_only=True)
    runMode = pset(0)
    exportGraphics = pset(False)
    runTriggers = pset([])
    triggerReady = pset({})
    inputConnectionsStore = pset({})
    optionsChecked = pset({})
    inputFile = pset("")

    def __init__(self):
        super().__init__(self.docker_image_name, self.docker_image_tag)
        with open(getJsonName(__file__, "gnumeric")) as f:
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
