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

class OWStart(OWBwBWidget):
    name = "Start"
    description = "Enter workflow parameters and start"
    priority = 10
    icon = getIconName(__file__,"start.png")
    want_main_area = False
    docker_image_name = "biodepot/salmon-start-simple"
    docker_image_tag = "alpine_3.18"
    outputs = [("work_dir",str),("index",str),("mate_1",str),("mate_2",str),("unpaired",str),("sample_dl_links",str),("output",str)]
    pset=functools.partial(settings.Setting,schema_only=True)
    runMode=pset(0)
    exportGraphics=pset(False)
    runTriggers=pset([])
    triggerReady=pset({})
    inputConnectionsStore=pset({})
    optionsChecked=pset({})
    work_dir=pset(None)
    samples=pset([])
    index=pset(None)
    mate_1=pset([])
    mate_2=pset([])
    unpaired=pset([])
    sample_dl_links=pset([])
    output=pset([])
    def __init__(self):
        super().__init__(self.docker_image_name, self.docker_image_tag)
        with open(getJsonName(__file__,"Start")) as f:
            self.data=jsonpickle.decode(f.read())
            f.close()
        self.initVolumes()
        self.inputConnections = ConnectionDict(self.inputConnectionsStore)
        self.drawGUI()
    def handleOutputs(self):
        outputValue=None
        if hasattr(self,"work_dir"):
            outputValue=getattr(self,"work_dir")
        self.send("work_dir", outputValue)
        outputValue=None
        if hasattr(self,"index"):
            outputValue=getattr(self,"index")
        self.send("index", outputValue)
        outputValue=None
        if hasattr(self,"mate_1"):
            outputValue=getattr(self,"mate_1")
        self.send("mate_1", outputValue)
        outputValue=None
        if hasattr(self,"mate_2"):
            outputValue=getattr(self,"mate_2")
        self.send("mate_2", outputValue)
        outputValue=None
        if hasattr(self,"unpaired"):
            outputValue=getattr(self,"unpaired")
        self.send("unpaired", outputValue)
        outputValue=None
        if hasattr(self,"sample_dl_links"):
            outputValue=getattr(self,"sample_dl_links")
        self.send("sample_dl_links", outputValue)
        outputValue=None
        if hasattr(self,"output"):
            outputValue=getattr(self,"output")
        self.send("output", outputValue)
