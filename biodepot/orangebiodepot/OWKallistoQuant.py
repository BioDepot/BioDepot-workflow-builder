import os
import glob
import sys
from collections import OrderedDict
from Orange.widgets import widget, gui, settings
from orangebiodepot.util.DockerClient import DockerClient
from orangebiodepot.util.BwBase import OWBwBWidget, ConnectionDict, ContainerPaths, BwbGuiElement, BwbGuiValue
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
    inputs = [("fastqDir", str, 'setfastqDir', widget.Multiple),
              ("indexFile", str, 'setindexFile')]
    outputs = [("outputDir", str)]
    
    #persistent settings
    outputDir = settings.Setting('/data', schema_only=True)
    indexFile = settings.Setting('/data/kallisto_hg38.idx', schema_only=True)
    fastqDir = settings.Setting('/data/fastq', schema_only=True)
#    pseudoBam = settings.Setting(False, schema_only=True)
    nThreads = settings.Setting(1, schema_only=True)
    nThreadsChecked = settings.Setting(False, schema_only=True)
    nBootstraps=settings.Setting(0, schema_only=True)
    nBootstrapsChecked = settings.Setting(False, schema_only=True)
    #this is a dict of lists of form {inputSlot:(widgetId1,widgetId2..widgetIdN) where None is a valid entry for the case when there are no Multiple connections }
    inputConnectionsStore = settings.Setting({},schema_only=True);
    
    def __init__(self):
        super().__init__(self.docker_image_name, self.docker_image_tag)
        self.con=ContainerPaths()
        self.con.outputDir="/root/output"
        self.con.fastqDir="/root/fastq"
        self.con.indexDir="/root/index"
        self.inputConnections=ConnectionDict(self.inputConnectionsStore)
        self.setDirectories(self.con.outputDir,self.outputDir)
        #init gui elements
        self.ledit=BwbGuiElement()
        self.button=BwbGuiElement()
        self.checkbox=BwbGuiElement()
        self.checkboxValue=BwbGuiValue()
        self.spin=BwbGuiElement()
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
            'inputs'   : OrderedDict([ ('indexFile', str),
                                       ('fastqDir', str)
                                    ]),
            'outputs'  : OrderedDict([ ('outputDir', str)
                                    ]),
            'requiredParameters' : ['outputDir','indexFile','fastqDir'],
            'parameters': OrderedDict([('outputDir',{'flags':['-i','--output-dir='], 'label':'Output directory', 'type': "directory"}),
                                       ('indexFile',{'flags':['-i','--index='], 'label':'Index file', 'type':'file'}),
                                       ('fastqDir' ,{'flags':[],'label':'fastq directory', 'type': "directory",'filePattern':["*.fastq.*"]}),
                                       ('bias',     {'flags':['--bias'],'label':'Perform sequence based bias correction','type': 'binary'}),
                                       ('bootstrap',{'flags':['-b','--bootstrap-samples='],'label':'Number of bootstrap samples','type': 'int','default':0}),
                                       ('seed',     {'flags':['--seed='],'label':'Seed for bootstrap sampling','type': 'int','gui':'Ledit', 'default' : 42}),
                                       ('plaintext',{'flags':['--plaintext'],'label':'Output plaintext instead of HDF5','type': 'bool', 'default' : False}),                                         
                                       ('fusion',   {'flags':['--fusion'],'label':'Search for fusion genes','type': 'bool', 'default' : False}),
                                       ('single',   {'flags':['--single'],'label':'Quantify single-end reads','type': 'bool', 'default' : False}),    
                                       ('single-overhang',   {'flags':['--single-overhang'],'label':'Include reads that go beyond transcript start','type': 'bool', 'default' : False}),                                       
                                       ('fr-stranded',   {'flags':['--fr-stranded'],'label':'strand specific read - first read forward','type': 'bool', 'default' : False}),
                                       ('rf-stranded',   {'flags':['--rf-stranded'],'label':'strand specific read - first read reverse','type': 'bool', 'default' : False}),
                                       ('fragment-length',   {'flags':['-l','--fragment-length'],'label':'estimated fragment length','type': 'double', 'default' : None}),
                                       ('stdev',    {'flags':['-s','--sd'],'label':'standard deviation of fragment length','type': 'double', 'default' : None}),
                                       ('nThreads' ,{'flags':['-t','--threads='],'label':'Number of threads','type': 'int','default':1}),
                                       ('pseudoBam',{'flags':['--pseudobam'],'label':'Save alignments to BAM file','type': 'bool','default':False}),                                                                            
                                       ('genomeBam',{'flags':['--genomebam'],'label':'Project alignments to sorted BAM file','type': 'bool','default':False}),
                                       ('gtf' ,{'flags':['--gtf'],'label':'GTF file', 'type': "file"}),                                                                                
                                       ('chromosomes' ,{'flags':['-c','--chromosomes'],'label':'Chromosome file', 'type': "file"})
                           ]
                          )
        }
        self.drawGUI()

    def OnRunClicked(self):
        self.startJob()

    def OnSettingsChanged(self):
        self.drawGui()

    def startJob(self):
        volumes = {}
        for conDir in (self.con.outputDir,self.con.fastqDir,self.con.indexDir):
            if not self.getDirectory(conDir):
                self.infoLabel.setText("Incorrect mapping of directory to " + conDir)
                return
            volumes[conDir]=self.getDirectory(conDir)
        
        flags={"-i":self.con.indexFile,"-o":self.con.outputDir}
        if self.nThreadsChecked :
            flags['-t']=str(self.nThreads)
        if self.nBootstrapsChecked :
            flags['-b']=str(self.nBootstraps)
        logFile=os.path.join(self.con.outputDir,'kallistoLog')
        bareCmd=self.generateCmd(('kallisto', 'quant'),flags,args=[])
        cmd="bash -c 'find "+ self.con.fastqDir+"/*.fastq.* -type f | sort | xargs " + bareCmd + " >& {}'".format(logFile)
#        self.infoLabel.setText(cmd)
        self.dockerRun(volumes,cmd)

    def Event_OnRunFinished(self):
        self.infoLabel.setText("Finished")
        self.send("outputDir", self.outputDir)
        #self.send("runTrigger", "done")
        
    def Event_OnRunMessage(self, message):
        self.infoLabel.setText(message)
    
    def setfastqDir(self, path, sourceId=None):
        self.setDir(path,'fastqDir',sourceId=sourceId)

    def setindexFile(self, path, sourceId=None):
        self.setFile(path,'indexFile',sourceId=sourceId)
    
