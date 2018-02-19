import os
import glob
import sys
from collections import OrderedDict
from Orange.widgets import widget, gui, settings
from orangebiodepot.util.DockerClient import DockerClient
from orangebiodepot.util.BwBase import OWBwBWidget, ConnectionDict, ContainerPaths, BwbGuiElements
from PyQt5 import QtWidgets, QtGui


class OWKallistoQuant(OWBwBWidget):    
    name = "Kallisto Quant"
    description = 'Alignment and quantification of reads from fastq files'
    category = "RNASeq"
    icon = "icons/kallisto-analysis.svg"
    priority = 10
    want_main_area = False
    docker_image_name = "biodepot/kallisto"
    docker_image_tag = "latest"
    setting=settings.Setting
    inputs = [("fastqFiles", str, '_set_fastqFiles', widget.Multiple),
              ("indexFile", str, '_set_indexFile')]
    outputs = [("outputDir", str)]
    
    #persistent settings
    outputDir = settings.Setting('/data', schema_only=True)
    indexFile = settings.Setting('/data/kallisto_hg38.idx', schema_only=True)
    fastqFiles = settings.Setting('/data/fastq', schema_only=True)
    nThreads = settings.Setting(1, schema_only=True)
    nThreadsChecked = settings.Setting(False, schema_only=True)
    nBootstraps=settings.Setting(0, schema_only=True)
    nBootstrapsChecked = settings.Setting(False, schema_only=True)
    runMode=settings.Setting(0, schema_only=True)
    runTriggers=settings.Setting([], schema_only=True)
    #this is a dict of lists of form {inputSlot:(widgetId1,widgetId2..widgetIdN) where None is a valid entry for the case when there are no Multiple connections }
    inputConnectionsStore = settings.Setting({},schema_only=True);
    
    def __init__(self):
        super().__init__(self.docker_image_name, self.docker_image_tag)
        self.data={
            'name':'Kallisto Quant',
            'description':'Alignment and quantification of reads from fastq files',
            'category' : "RNASeq",
            'icon' : "icons/kallisto-analysis.svg",
            'priority' : 10,
            'want_main_area' : False,
            'docker_image_name' : 'biodepot/kallisto',
            'docker_image_tag' : 'latest',
            'persistentSettings' : 'all',
            'command' : 'kallisto quant',
            'inputs'   : OrderedDict([ ('indexFile', {'type':str, 'callback' : None}),
                                       ('fastqFiles',  {'type':str, 'callback' : None})
                                    ]),
            'outputs'  : OrderedDict([ ('outputDir', {'type':str, 'output' : {'value': 'outputDir', 'type': str, 'valueType' :'attribute'}, 'callback' : None})
                                    ]),
            'volumeMappings' : [{'conVolume':'/root/output','attr':'outputDir','default':'/data/output'},
                                {'conVolume':'/root/fastq' , 'attr':'fastqFiles'},
                                {'conVolume':'/root/reference', 'attr':'indexFile','default':'/data/reference'}
                               ],
            'requiredParameters' : ['outputDir','indexFile','fastqFiles'],
            'parameters': OrderedDict([('outputDir',{'flags':['-o','--output-dir='], 'label':'Output directory', 'type': "directory"}),
                                       ('indexFile',{'flags':['-i','--index='], 'label':'Index file', 'type':'file'}),
                                       ('fastqFiles' ,{'flags':[],'label':'fastq files', 'type': "files",'filePattern':["*.fastq.*"]}),
                                       ('bias',     {'flags':['--bias'],'label':'Perform sequence based bias correction','type': 'binary'}),
                                       ('bootstrap',{'flags':['-b','--bootstrap-samples='],'label':'Number of bootstrap samples','type': 'int','default': 1}),
                                       ('seed',     {'flags':['--seed='],'label':'Seed for bootstrap sampling','type': 'int','gui':'Ledit', 'default' : 42}),
                                       ('plaintext',{'flags':['--plaintext'],'label':'Output plaintext instead of HDF5','type': 'bool', 'default' : False}),                                         
                                       ('fusion',   {'flags':['--fusion'],'label':'Search for fusion genes','type': 'bool', 'default' : False}),
                                       ('single',   {'flags':['--single'],'label':'Quantify single-end reads','type': 'bool', 'default' : False}),    
                                       ('single-overhang',   {'flags':['--single-overhang'],'label':'Include reads that go beyond transcript start','type': 'bool', 'default' : False}),                                       
                                       ('fr-stranded',   {'flags':['--fr-stranded'],'label':'strand specific read - first read forward','type': 'bool', 'default' : False}),
                                       ('rf-stranded',   {'flags':['--rf-stranded'],'label':'strand specific read - first read reverse','type': 'bool', 'default' : False}),
                                       ('fragment-length',   {'flags':['-l','--fragment-length'],'label':'Estimated fragment length','type': 'double', 'default' : None}),
                                       ('stdev',    {'flags':['-s','--sd'],'label':'Standard deviation of fragment length','type': 'double', 'default' : None}),
                                       ('nThreads' ,{'flags':['-t','--threads='],'label':'Number of threads','type': 'int','default':1}),
                                       ('pseudoBam',{'flags':['--pseudobam'],'label':'Save alignments to BAM file','type': 'bool','default':False}),                                                                            
                                       ('genomeBam',{'flags':['--genomebam'],'label':'Project alignments to sorted BAM file','type': 'bool','default':False}),
                                       ('gtf' ,{'flags':['--gtf'],'label':'GTF file', 'type': "file"}),                                                                                
                                       ('chromosomes' ,{'flags':['-c','--chromosomes'],'label':'Chromosome file', 'type': "file"})
                           ]
                          )
        }
        self.initVolumes()
        self.inputConnections=ConnectionDict(self.inputConnectionsStore)
        self.drawGUI()
        
    #input callbacks
    def _set_fastqFiles(self, path, sourceId=None):
        self.handleInputs(path,'fastqFiles',sourceId=sourceId)

    def _set_indexFile(self, path, sourceId=None):
        self.handleInputs(path,'indexFile',sourceId=sourceId)

