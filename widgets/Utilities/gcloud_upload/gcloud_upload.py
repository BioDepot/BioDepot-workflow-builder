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

class OWgcloud_upload(OWBwBWidget):
    name = "gcloud_upload"
    description = "Upload a directory to gcp bucket"
    priority = 10
    icon = getIconName(__file__,"gcloud_upload.png")
    want_main_area = False
    docker_image_name = "biodepot/gcpupload"
    docker_image_tag = "277.0.0-alpine__173ca532"
    inputs = [("Trigger",str,"handleInputsTrigger"),("credentials_file",str,"handleInputscredentials_file"),("bucket",str,"handleInputsbucket"),("bucketDir",str,"handleInputsbucketDir"),("uploadDir",str,"handleInputsuploadDir")]
    outputs = [("credentials_file",str)]
    pset=functools.partial(settings.Setting,schema_only=True)
    runMode=pset(0)
    exportGraphics=pset(False)
    runTriggers=pset([])
    triggerReady=pset({})
    inputConnectionsStore=pset({})
    optionsChecked=pset({})
    credentials_file=pset("/data/credentials.json")
    bucket=pset("myBucket")
    uploadDir=pset("/data/source")
    bucketDir=pset("destination")
    def __init__(self):
        super().__init__(self.docker_image_name, self.docker_image_tag)
        with open(getJsonName(__file__,"gcloud_upload")) as f:
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
    def handleInputsbucketDir(self, value, *args):
        if args and len(args) > 0: 
            self.handleInputs("bucketDir", value, args[0][0], test=args[0][3])
        else:
            self.handleInputs("inputFile", value, None, False)
    def handleInputsuploadDir(self, value, *args):
        if args and len(args) > 0: 
            self.handleInputs("uploadDir", value, args[0][0], test=args[0][3])
        else:
            self.handleInputs("inputFile", value, None, False)
    def handleOutputs(self):
        outputValue=None
        if hasattr(self,"credentials_file"):
            outputValue=getattr(self,"credentials_file")
        self.send("credentials_file", outputValue)
