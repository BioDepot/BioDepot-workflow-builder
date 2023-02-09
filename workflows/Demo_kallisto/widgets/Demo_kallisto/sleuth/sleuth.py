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


class OWsleuth(OWBwBWidget):
    name = "sleuth"
    description = "Analysis of kallisto quantified RNA-seq files"
    priority = 5
    icon = getIconName(__file__, "sleuth2.png")
    want_main_area = False
    docker_image_name = "biodepot/sleuth"
    docker_image_tag = "0.30.1__bioc-r_3.16-ubuntu-22.04-r-4.2.2__face979c__4bc3b0b7__394dc8fc"
    inputs = [("trigger", str, "handleInputstrigger")]
    outputs = [("output_file", str)]
    pset = functools.partial(settings.Setting, schema_only=True)
    runMode = pset(0)
    exportGraphics = pset(False)
    runTriggers = pset([])
    triggerReady = pset({})
    inputConnectionsStore = pset({})
    optionsChecked = pset({})
    experimental_description = pset(None)
    full_model = pset(None)
    filter_fun = pset(None)
    norm_fun_counts = pset(None)
    norm_fun_tpm = pset(None)
    aggregation_column = pset("ext_gene")
    transformation_function = pset(None)
    num_cores = pset(1)
    read_bootstrap_tpm = pset(False)
    extra_bootstrap_summary = pset(False)
    column = pset(None)
    wald = pset(False)
    nTopGenes = pset(40)
    output_file = pset(None)
    qvalue = pset(0.05)
    geneNamesFile = pset(None)

    def __init__(self):
        super().__init__(self.docker_image_name, self.docker_image_tag)
        with open(getJsonName(__file__, "sleuth")) as f:
            self.data = jsonpickle.decode(f.read())
            f.close()
        self.initVolumes()
        self.inputConnections = ConnectionDict(self.inputConnectionsStore)
        self.drawGUI()

    def handleInputstrigger(self, value, *args):
        if args and len(args) > 0:
            self.handleInputs("trigger", value, args[0][0], test=args[0][3])
        else:
            self.handleInputs("inputFile", value, None)

    def handleOutputs(self):
        outputValue = None
        if hasattr(self, "output_file"):
            outputValue = getattr(self, "output_file")
        self.send("output_file", outputValue)
