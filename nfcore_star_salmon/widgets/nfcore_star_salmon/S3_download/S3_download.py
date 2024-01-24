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

class OWS3_download(OWBwBWidget):
    name = "S3_download"
    description = "Enter and output a file"
    priority = 10
    icon = getIconName(__file__,"downloadS3.png")
    want_main_area = False
    docker_image_name = "biodepot/s3download"
    docker_image_tag = "1.16.272__python_3.8.0__alpine-3.10__e2962454"
    inputs = [("Trigger",str,"handleInputsTrigger"),("awsdir",str,"handleInputsawsdir"),("bucket",str,"handleInputsbucket"),("downloadDir",str,"handleInputsdownloadDir")]
    outputs = [("downloadDir",str)]
    pset=functools.partial(settings.Setting,schema_only=True)
    runMode=pset(0)
    exportGraphics=pset(False)
    runTriggers=pset([])
    triggerReady=pset({})
    inputConnectionsStore=pset({})
    optionsChecked=pset({})
    awsdir=pset("/data/.aws")
    bucket=pset("myBucket")
    downloadDir=pset("/data")
    dirs=pset([])
    nthreads=pset(None)
    key=pset(None)
    secret=pset(None)
    region=pset("us-east-2")
    def __init__(self):
        super().__init__(self.docker_image_name, self.docker_image_tag)
        with open(getJsonName(__file__,"S3_download")) as f:
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
    def handleInputsawsdir(self, value, *args):
        if args and len(args) > 0: 
            self.handleInputs("awsdir", value, args[0][0], test=args[0][3])
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
