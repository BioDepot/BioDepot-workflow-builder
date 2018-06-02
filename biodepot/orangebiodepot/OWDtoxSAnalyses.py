import os
import glob
import sys
import functools
import jsonpickle
from collections import OrderedDict
from Orange.widgets import widget, gui, settings
import Orange.data
from Orange.data.io import FileFormat
from orangebiodepot.util.DockerClient import DockerClient
from orangebiodepot.util.BwBase import OWBwBWidget, ConnectionDict, BwbGuiElements
from PyQt5 import QtWidgets, QtGui

class OWDtoxsAnalysis(OWBwBWidget):
    name = "Dtoxs Analysis"
    description = "Step 2 of Dtoxs SOP. Uses edgeR for differential expression analysis"
    category = "RNA-seq"
    priority = 10
    icon = "/biodepot/orangebiodepot/icons/dtoxs-analysis2.svg"
    want_main_area = False
    docker_image_name = "biodepot/dtoxs_analysis"
    docker_image_tag = "1.0"
    inputs = [("RepositoryDirectory",str,"handleInputsRepositoryDirectory"),("ConfigurationFile",str,"handleInputsConfigurationFile"),("Trigger",str,"handleInputsTrigger")]
    outputs = [("ResultsDirectory",str),("Top40",Orange.data.Table)]
    pset=functools.partial(settings.Setting,schema_only=True)
    runMode=pset(0)
    runTriggers=pset([])
    triggerReady=pset({})
    inputConnectionsStore=pset({})
    optionsChecked=pset({})
    RepositoryDirectory=pset(None)
    ConfigurationFile=pset(None)
    def __init__(self):
        super().__init__(self.docker_image_name, self.docker_image_tag)
        with open("/biodepot/orangebiodepot/json/DtoxSAnalyses.json") as f:
            self.data=jsonpickle.decode(f.read())
            f.close()
        self.initVolumes()
        self.inputConnections = ConnectionDict(self.inputConnectionsStore)
        self.drawGUI()
    def handleInputsRepositoryDirectory(self, value, sourceId=None):
        self.handleInputs(value, "RepositoryDirectory", sourceId=None)
    def handleInputsConfigurationFile(self, value, sourceId=None):
        self.handleInputs(value, "ConfigurationFile", sourceId=None)
    def handleInputsTrigger(self, value, sourceId=None):
        self.handleInputs(value, "Trigger", sourceId=None)
    def handleOutputs(self):
        resultsDir=os.path.join(getattr(self,"RepositoryDirectory"), 'Results');
        self.send("ResultsDirectory",resultsDir)
        tsvFile = os.path.join(resultsDir, 'FDR-0.1/TOP-40.tsv');
        tsvReader = FileFormat.get_reader(tsvFile)
        if os.path.isfile(tsvFile):
            try:
                data = tsvReader.read()
            except Exception as e:
                sys.stderr.write('unable to read output data exception {}\n'.format(e.message))
            return
            self.send("Top40", data)
