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

class OWrsem(OWBwBWidget):
    name = "rsem"
    description = "Rsem"
    priority = 1
    icon = getIconName(__file__,"rsem.png")
    want_main_area = False
    docker_image_name = "quay.io/ucsc_cgl/rsem"
    docker_image_tag = "1.2.25--d4275175cc8df36967db460b06337a14f40d2f21"
    inputs = [("RsemIndexDir",str,"handleInputsRsemIndexDir"),("BamFile",str,"handleInputsBamFile"),("Trigger",str,"handleInputsTrigger")]
    outputs = [("OutputDir",str)]
    pset=functools.partial(settings.Setting,schema_only=True)
    runMode=pset(0)
    exportGraphics=pset(False)
    runTriggers=pset([])
    triggerReady=pset({})
    inputConnectionsStore=pset({})
    optionsChecked=pset({})
    PairedEnd=pset(None)
    NoQualities=pset(None)
    Threads=pset(None)
    ForwardProb=pset(None)
    SeedLength=pset(None)
    FragmentLengthMean=pset(None)
    Bam=pset(None)
    ReferencePrefix=pset(None)
    OutputPrefix=pset(None)
    def __init__(self):
        super().__init__(self.docker_image_name, self.docker_image_tag)
        with open(getJsonName(__file__,"rsem")) as f:
            self.data=jsonpickle.decode(f.read())
            f.close()
        self.initVolumes()
        self.inputConnections = ConnectionDict(self.inputConnectionsStore)
        self.drawGUI()
    def handleInputsRsemIndexDir(self, value, *args):
        if args and len(args) > 0: 
            self.handleInputs("RsemIndexDir", value, args[0][0], test=args[0][3])
        else:
            self.handleInputs("inputFile", value, None, False)
    def handleInputsBamFile(self, value, *args):
        if args and len(args) > 0: 
            self.handleInputs("BamFile", value, args[0][0], test=args[0][3])
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
