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

class OWnextflow(OWBwBWidget):
    name = "nextflow"
    description = "Minimum perl container"
    priority = 20
    icon = getIconName(__file__,"nfcore1.png")
    want_main_area = False
    docker_image_name = "nextflow/nextflow"
    docker_image_tag = "23.04.2"
    inputs = [("inputFile",str,"handleInputsinputFile"),("Trigger",str,"handleInputsTrigger")]
    outputs = [("OutputDir",str)]
    pset=functools.partial(settings.Setting,schema_only=True)
    runMode=pset(0)
    exportGraphics=pset(False)
    runTriggers=pset([])
    triggerReady=pset({})
    inputConnectionsStore=pset({})
    optionsChecked=pset({})
    nxfusrmap=pset("1000")
    nxfdockeropts=pset("'-u 1000'")
    workingdir=pset("  /data")
    nxfassets=pset("/data/.nextflow/assets")
    configfiles=pset([])
    jvmprops=pset(None)
    background1=pset(False)
    config=pset(None)
    logfile=pset(None)
    quiet=pset(False)
    command=pset("run")
    repo=pset("nf-core/rnaseq ")
    profile=pset("test,docker")
    parametersfile=pset(None)
    script=pset(None)
    withdocker=pset(" -with-docker")
    pipeline=pset(None)
    outdir=pset("")
    def __init__(self):
        super().__init__(self.docker_image_name, self.docker_image_tag)
        with open(getJsonName(__file__,"nextflow")) as f:
            self.data=jsonpickle.decode(f.read())
            f.close()
        self.initVolumes()
        self.inputConnections = ConnectionDict(self.inputConnectionsStore)
        self.drawGUI()
    def handleInputsinputFile(self, value, *args):
        if args and len(args) > 0: 
            self.handleInputs("inputFile", value, args[0][0], test=args[0][3])
        else:
            self.handleInputs("inputFile", value, None, False)
    def handleInputsTrigger(self, value, *args):
        if args and len(args) > 0: 
            self.handleInputs("Trigger", value, args[0][0], test=args[0][3])
        else:
            self.handleInputs("inputFile", value, None, False)
    def handleOutputs(self):
        outputValue=None
        if hasattr(self,"OutputDir"):
            outputValue=getattr(self,"OutputDir")
        self.send("OutputDir", outputValue)
