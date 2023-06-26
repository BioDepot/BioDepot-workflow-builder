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

class OWgcloud_download(OWBwBWidget):
    name = "gcloud_download"
    description = "Download google bucket files to a local directory"
    priority = 10
    icon = getIconName(__file__,"gcloud_download.png")
    want_main_area = False
    docker_image_name = "biodepot/gcpdownload"
    docker_image_tag = "277.0.0-alpine__1c1fdc4b"
    inputs = [("Trigger",str,"handleInputsTrigger"),("credentials_file",str,"handleInputscredentials_file"),("bucket",str,"handleInputsbucket"),("downloadDir",str,"handleInputsdownloadDir")]
    outputs = [("downloadDir",str)]
    pset=functools.partial(settings.Setting,schema_only=True)
    runMode=pset(0)
    exportGraphics=pset(False)
    runTriggers=pset([])
    triggerReady=pset({})
    inputConnectionsStore=pset({})
    optionsChecked=pset({})
    credentials_file=pset("/data/credentials.json")
    bucket=pset("myBucket")
    downloadDir=pset("/data")
    dirs=pset([])
    files=pset([])
    def __init__(self):
        super().__init__(self.docker_image_name, self.docker_image_tag)
        with open(getJsonName(__file__,"gcloud_download")) as f:
            self.data=jsonpickle.decode(f.read())
            f.close()
        self.initVolumes()
        self.inputConnections = ConnectionDict(self.inputConnectionsStore)
        self.drawGUI()
    def handleInputsTrigger(self, value, *args):
        if args and len(args) > 0: 
            self.handleInputs("Trigger", value, args[0][0], test=args[0][3])
        else:
            self.handleInputs("inputFile", value, None, False)
    def handleInputscredentials_file(self, value, *args):
        if args and len(args) > 0: 
            self.handleInputs("credentials_file", value, args[0][0], test=args[0][3])
        else:
            self.handleInputs("inputFile", value, None, False)
    def handleInputsbucket(self, value, *args):
        if args and len(args) > 0: 
            self.handleInputs("bucket", value, args[0][0], test=args[0][3])
        else:
            self.handleInputs("inputFile", value, None, False)
    def handleInputsdownloadDir(self, value, *args):
        if args and len(args) > 0: 
            self.handleInputs("downloadDir", value, args[0][0], test=args[0][3])
        else:
            self.handleInputs("inputFile", value, None, False)
    def handleOutputs(self):
        outputValue=None
        if hasattr(self,"downloadDir"):
            outputValue=getattr(self,"downloadDir")
        self.send("downloadDir", outputValue)
