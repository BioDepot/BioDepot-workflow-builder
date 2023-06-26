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

class OWS3_upload(OWBwBWidget):
    name = "S3_upload"
    description = "Enter and output a file"
    priority = 10
    icon = getIconName(__file__,"cloud_upload.png")
    want_main_area = False
    docker_image_name = "biodepot/s3upload"
    docker_image_tag = "1.16.272__python_3.8.0__alpine-3.10__381c2b66"
    inputs = [("Trigger",str,"handleInputsTrigger"),("credentials_dir",str,"handleInputscredentials_dir"),("uploadDir",str,"handleInputsuploadDir"),("bucket",str,"handleInputsbucket"),("s3Dir",str,"handleInputss3Dir")]
    outputs = [("credentials_dir",str)]
    pset=functools.partial(settings.Setting,schema_only=True)
    runMode=pset(0)
    exportGraphics=pset(False)
    runTriggers=pset([])
    triggerReady=pset({})
    inputConnectionsStore=pset({})
    optionsChecked=pset({})
    credentials_dir=pset("/data/.aws")
    uploadDir=pset("/data/source")
    bucket=pset("myBucket")
    s3Dir=pset("destination")
    nthreads=pset(None)
    key=pset(None)
    secret=pset(None)
    region=pset("us-east-2")
    def __init__(self):
        super().__init__(self.docker_image_name, self.docker_image_tag)
        with open(getJsonName(__file__,"S3_upload")) as f:
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
    def handleInputscredentials_dir(self, value, *args):
        if args and len(args) > 0: 
            self.handleInputs("credentials_dir", value, args[0][0], test=args[0][3])
        else:
            self.handleInputs("inputFile", value, None, False)
    def handleInputsuploadDir(self, value, *args):
        if args and len(args) > 0: 
            self.handleInputs("uploadDir", value, args[0][0], test=args[0][3])
        else:
            self.handleInputs("inputFile", value, None, False)
    def handleInputsbucket(self, value, *args):
        if args and len(args) > 0: 
            self.handleInputs("bucket", value, args[0][0], test=args[0][3])
        else:
            self.handleInputs("inputFile", value, None, False)
    def handleInputss3Dir(self, value, *args):
        if args and len(args) > 0: 
            self.handleInputs("s3Dir", value, args[0][0], test=args[0][3])
        else:
            self.handleInputs("inputFile", value, None, False)
    def handleOutputs(self):
        outputValue=None
        if hasattr(self,"credentials_dir"):
            outputValue=getattr(self,"credentials_dir")
        self.send("credentials_dir", outputValue)
