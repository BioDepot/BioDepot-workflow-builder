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

class OWgen3_download(OWBwBWidget):
    name = "gen3_download"
    description = "Download using GEN3 Client or GDC API"
    priority = 10
    icon = getIconName(__file__,"gen3-download.png")
    want_main_area = False
    docker_image_name = "biodepot/gen3-download"
    docker_image_tag = "2021.03__alpine_3.12__833f38e5"
    inputs = [("manifest",str,"handleInputsmanifest"),("guids",str,"handleInputsguids"),("downloadDir",str,"handleInputsdownloadDir"),("Trigger",str,"handleInputsTrigger"),("credentials",str,"handleInputscredentials"),("gdctoken",str,"handleInputsgdctoken")]
    outputs = [("manifest",str),("guids",str),("downloadDir",str)]
    pset=functools.partial(settings.Setting,schema_only=True)
    runMode=pset(0)
    exportGraphics=pset(False)
    runTriggers=pset([])
    triggerReady=pset({})
    inputConnectionsStore=pset({})
    optionsChecked=pset({})
    profile=pset("user")
    credentials=pset("/data/credentials.json")
    downloadDir=pset("/data")
    guids=pset([])
    manifest=pset(None)
    filenameformat=pset("original")
    rename=pset(False)
    skipcompleted=pset(True)
    numparallel=pset(1)
    protocol=pset(None)
    decompress=pset(True)
    gdctoken=pset(None)
    datacommons_url=pset(None)
    def __init__(self):
        super().__init__(self.docker_image_name, self.docker_image_tag)
        with open(getJsonName(__file__,"gen3_download")) as f:
            self.data=jsonpickle.decode(f.read())
            f.close()
        self.initVolumes()
        self.inputConnections = ConnectionDict(self.inputConnectionsStore)
        self.drawGUI()
    def handleInputsmanifest(self, value, *args):
        if args and len(args) > 0: 
            self.handleInputs("manifest", value, args[0][0], test=args[0][3])
        else:
            self.handleInputs("inputFile", value, None, False)
    def handleInputsguids(self, value, *args):
        if args and len(args) > 0: 
            self.handleInputs("guids", value, args[0][0], test=args[0][3])
        else:
            self.handleInputs("inputFile", value, None, False)
    def handleInputsdownloadDir(self, value, *args):
        if args and len(args) > 0: 
            self.handleInputs("downloadDir", value, args[0][0], test=args[0][3])
        else:
            self.handleInputs("inputFile", value, None, False)
    def handleInputsTrigger(self, value, *args):
        if args and len(args) > 0: 
            self.handleInputs("Trigger", value, args[0][0], test=args[0][3])
        else:
            self.handleInputs("inputFile", value, None, False)
    def handleInputscredentials(self, value, *args):
        if args and len(args) > 0: 
            self.handleInputs("credentials", value, args[0][0], test=args[0][3])
        else:
            self.handleInputs("inputFile", value, None, False)
    def handleInputsgdctoken(self, value, *args):
        if args and len(args) > 0: 
            self.handleInputs("gdctoken", value, args[0][0], test=args[0][3])
        else:
            self.handleInputs("inputFile", value, None, False)
    def handleOutputs(self):
        outputValue=None
        if hasattr(self,"manifest"):
            outputValue=getattr(self,"manifest")
        self.send("manifest", outputValue)
        outputValue=None
        if hasattr(self,"guids"):
            outputValue=getattr(self,"guids")
        self.send("guids", outputValue)
        outputValue=None
        if hasattr(self,"downloadDir"):
            outputValue=getattr(self,"downloadDir")
        self.send("downloadDir", outputValue)
