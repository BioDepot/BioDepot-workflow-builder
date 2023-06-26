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

class OWfiji(OWBwBWidget):
    name = "fiji"
    description = "Run a FIJI script or macro."
    priority = 10
    icon = getIconName(__file__,"fiji.png")
    want_main_area = False
    docker_image_name = "biodepot/fiji"
    docker_image_tag = "20201104-1356__update20211210__ubuntu_20.04__1fa37497"
    inputs = [("fijidir",str,"handleInputsfijidir"),("installfiji",str,"handleInputsinstallfiji"),("trigger",str,"handleInputstrigger"),("imagefile",str,"handleInputsimagefile"),("pluginsdir",str,"handleInputspluginsdir"),("macrotrigger",str,"handleInputsmacrotrigger"),("datatrigger",str,"handleInputsdatatrigger")]
    outputs = [("fijidir",str),("installfiji",str)]
    pset=functools.partial(settings.Setting,schema_only=True)
    runMode=pset(0)
    exportGraphics=pset(False)
    runTriggers=pset([])
    triggerReady=pset({})
    inputConnectionsStore=pset({})
    optionsChecked=pset({})
    fijidir=pset(None)
    installfiji=pset(None)
    overwrite=pset(False)
    headless=pset(False)
    pluginsdir=pset(None)
    macrofile=pset(None)
    script=pset(False)
    updatefiji=pset(False)
    param=pset(None)
    quitnow=pset(False)
    addupdatesite=pset([])
    additional=pset(None)
    arguments=pset([])
    def __init__(self):
        super().__init__(self.docker_image_name, self.docker_image_tag)
        with open(getJsonName(__file__,"fiji")) as f:
            self.data=jsonpickle.decode(f.read())
            f.close()
        self.initVolumes()
        self.inputConnections = ConnectionDict(self.inputConnectionsStore)
        self.drawGUI()
    def handleInputsfijidir(self, value, *args):
        if args and len(args) > 0: 
            self.handleInputs("fijidir", value, args[0][0], test=args[0][3])
        else:
            self.handleInputs("inputFile", value, None, False)
    def handleInputsinstallfiji(self, value, *args):
        if args and len(args) > 0: 
            self.handleInputs("installfiji", value, args[0][0], test=args[0][3])
        else:
            self.handleInputs("inputFile", value, None, False)
    def handleInputstrigger(self, value, *args):
        if args and len(args) > 0: 
            self.handleInputs("trigger", value, args[0][0], test=args[0][3])
        else:
            self.handleInputs("inputFile", value, None, False)
    def handleInputsimagefile(self, value, *args):
        if args and len(args) > 0: 
            self.handleInputs("imagefile", value, args[0][0], test=args[0][3])
        else:
            self.handleInputs("inputFile", value, None, False)
    def handleInputspluginsdir(self, value, *args):
        if args and len(args) > 0: 
            self.handleInputs("pluginsdir", value, args[0][0], test=args[0][3])
        else:
            self.handleInputs("inputFile", value, None, False)
    def handleInputsmacrotrigger(self, value, *args):
        if args and len(args) > 0: 
            self.handleInputs("macrotrigger", value, args[0][0], test=args[0][3])
        else:
            self.handleInputs("inputFile", value, None, False)
    def handleInputsdatatrigger(self, value, *args):
        if args and len(args) > 0: 
            self.handleInputs("datatrigger", value, args[0][0], test=args[0][3])
        else:
            self.handleInputs("inputFile", value, None, False)
    def handleOutputs(self):
        outputValue=None
        if hasattr(self,"fijidir"):
            outputValue=getattr(self,"fijidir")
        self.send("fijidir", outputValue)
        outputValue=None
        if hasattr(self,"installfiji"):
            outputValue=getattr(self,"installfiji")
        self.send("installfiji", outputValue)
