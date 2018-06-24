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

class OWsleuth(OWBwBWidget):
    name = "sleuth"
    description = "Analysis of kallisto quantified RNA-seq files"
    category = "RNA-seq"
    priority = 10
    icon = "/biodepot/RNA_seq/icons/sleuth2.png"
    want_main_area = False
    docker_image_name = "biodepot/ubuntu-sleuth"
    docker_image_tag = "18.04-4.5.1-py27-0.29"
    inputs = [("trigger",str,"handleInputstrigger")]
    outputs = [("output_file",str),("topGenes",Orange.data.Table)]
    pset=functools.partial(settings.Setting,schema_only=True)
    runMode=pset(0)
    runTriggers=pset([])
    triggerReady=pset({})
    inputConnectionsStore=pset({})
    optionsChecked=pset({})
    experimental_description=pset(None)
    full_model=pset(None)
    filter_fun=pset(None)
    norm_fun_counts=pset(None)
    norm_fun_tpm=pset(None)
    aggregation_column=pset("ext_gene")
    transformation_function=pset(None)
    num_cores=pset(1)
    read_bootstrap_tpm=pset(False)
    extra_bootstrap_summary=pset(False)
    column=pset(None)
    wald=pset(False)
    nTopGenes=pset(40)
    output_file=pset(None)
    qvalue=pset(0.05)
    geneNamesFile=pset(None)
    def __init__(self):
        super().__init__(self.docker_image_name, self.docker_image_tag)
        with open("/biodepot/RNA_seq/json/sleuth.json") as f:
            self.data=jsonpickle.decode(f.read())
            f.close()
        self.initVolumes()
        self.inputConnections = ConnectionDict(self.inputConnectionsStore)
        self.drawGUI()
    def handleInputstrigger(self, value, sourceId=None):
        self.handleInputs(value, "trigger", sourceId=None)
    def handleOutputs(self):
        outputValue=None
        if hasattr(self,"output_file"):
            outputValue=getattr(self,"output_file")
        self.send("output_file", outputValue)
        sys.stderr.write('output_file is {}\n'.format(outputValue))
        if outputValue:
            self.sendTable(outputValue)
        
    def sendTable(self,filename):
        outname=filename+'.tsv'
        
        with open (filename, 'r') as fin, open (outname,'w') as fout:
            line=fin.readline()
            fout.write(line)
            for line in fin:
                fout.write(line.split("\t",1)[1])
            fin.close()
            fout.close()
        tsvReader = FileFormat.get_reader(outname)
        dataTable= tsvReader.read()
        self.send("topGenes", dataTable)
        os.unlink(outname)


        
    
