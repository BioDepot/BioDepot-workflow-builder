import os
import glob
from Orange.widgets import widget, gui, settings
from orangebiodepot.util.DockerClient import DockerClient
from orangebiodepot.util.BwBase import OWBwBWidget
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

    inputs = [("hostFastqDir", str, "setHostFastqDir", widget.Multiple),
              ("indexFile", str, "setIndexFile")]
    outputs = [("hostOutputDir", str)]
    
    #persistent settings
    outputPath = settings.Setting('/data', schema_only=True)
    indexFile = settings.Setting('kallisto_hg38.idx', schema_only=True)
    hostFastqDir = settings.Setting('/data/fastq', schema_only=True)
    pseudoBam = settings.Setting(False, schema_only=True)
    nThreads = settings.Setting(1, schema_only=True)
    nThreadsChecked = settings.Setting(False, schema_only=True)
    indexFileConnected = settings.Setting(False, schema_only=True)
    hostFastqDirConnected = settings.Setting(False, schema_only=True)
    
    def __init__(self):
        super().__init__(self.docker_image_name, self.docker_image_tag)
        self.conOutputDir="/root/output"
        self.conFastqDir="/root/fastq"
        self.conIndexDir="/root/Index"
        self.setDirectories(self.conOutputDir,self.outputPath)
        self.drawGUI()

    def drawGUI(self):
        controlBox = gui.widgetBox(self.controlArea, "Status")
        self.infoLabel = gui.widgetLabel(controlBox, 'Waiting...')
        self.infoLabel.setWordWrap(True)
        optionsBox = gui.widgetBox(self.controlArea, "Options")
        outputBox = gui.vBox(optionsBox, "Choose output directory")
        gui.lineEdit(outputBox, self, 'outputPath')
        browseOutputBtn=gui.button(outputBox, self, "☰ Browse", callback=self.browseOutputDir, autoDefault=True, width=80)
        
        self.indexFileBox = gui.vBox(optionsBox, "Choose index file",disabled=self.indexFileConnected)
        self.indexFileLedit=gui.lineEdit(self.indexFileBox, self, 'indexFile',disabled=self.indexFileConnected)
        self.browseIndexBtn=gui.button(self.indexFileBox, self, "☰ Browse", callback=self.browseIndexFile, autoDefault=True, width=80, disabled=self.indexFileConnected)
        
        self.hostFastqDirBox = gui.vBox(optionsBox, "Choose fastq directory",disabled=self.hostFastqDirConnected)
        self.hostFastqDirLedit=gui.lineEdit(self.hostFastqDirBox, self, 'hostFastqDir',disabled=self.hostFastqDirConnected)
        self.browseHostFastqDirBtn=gui.button(self.hostFastqDirBox, self, "☰ Browse", callback=self.browseHostFastqDir, autoDefault=True, width=80, disabled=self.hostFastqDirConnected)
        
        parmBox = gui.vBox(optionsBox, "Optional flags")
        gui.checkBox(parmBox, self, "pseudoBam", "Output pseudoBam files")
        gui.spin(parmBox, self, "nThreads", minv=1, maxv=128,
        label="Nthreads:", checked='nThreadsChecked')
        btnRun = gui.button(optionsBox, self, "Run", callback=self.OnRunClicked)

    def OnRunClicked(self):
        self.startJob()

    def OnSettingsChanged(self):
        self.drawGui()

    def startJob(self):
        volumes = {}
        vstr=""
        message=""
        for conDir in (self.conFastqDir,self.conOutputDir,self.conFastqDir):
            if not self.getDirectory(conDir):
                self.infoLabel.setText("Incorrect mapping of directory to " + conDir)
                return
            volumes[self.getDirectory(conDir)]=conDir
        volumes['/data/logs']='/root/logs'
        
        for h,c in volumes.items():
            message += " -v " + h+":"+c
        bareCmd=self.generateCmd(('kallisto', 'quant'),{"-i":self.indexFile,"-o":self.conOutputDir},args=[])
        cmd="bash -c 'find "+ self.conFastqDir+"/*.fastq.* -type f | sort xargs " + bareCmd + " >& /root/logs/log'"
        self.infoLabel.setText(message+" "+cmd)
        self.dockerRun(volumes,cmd)

    def Event_OnRunFinished(self):
        self.infoLabel.setText("Finished")
        self.send("hostOutputDir", self.outputPath)
        #self.send("runTrigger", "done")
        
    def Event_OnRunMessage(self, message):
        self.infoLabel.setText(message)

    def setHostFastqDir(self, path, sourceId=None):
        #This only works for single connections
        #change to Dict to keep track of connections
        if path is None:
            self.hostFastqDirConnected = False
            self.hostFastqDir = None
            self.setDirectories(self.conFastqDir, None)
        else:
            self.hostFastqDir=path
            self.hostFastqDirConnected = True
            self.setDirectories(self.conFastqDir,path)
        self.browseHostFastqDirBtn.setEnabled(not self.hostFastqDirConnected)
        self.hostFastqDirLedit.setEnabled(not self.hostFastqDirConnected)
        self.hostFastqDirBox.setEnabled(not self.hostFastqDirConnected)
        self.setDirectories(self.conFastqDir,path)

    def setIndexFile(self, path, sourceId=None):
        if path is None:
            self.indexFileConnected = False
            self.indexFile = None
            self.hostIndexDir = None
            self.setDirectories(self.conIndexDir, None)
        else:
            self.indexFile=path
            self.hostIndexDir=os.path.dirname(self.indexFile)
            self.indexFileConnected = True
            self.setDirectories(self.conIndexDir,path)
        self.browseIndexBtn.setEnabled(not self.indexFileConnected)
        self.indexFileLedit.setEnabled(not self.indexFileConnected)
        self.indexFileBox.setEnabled(not self.indexFileConnected)
        
    def browseOutputDir(self):
        defaultDir = '/root'
        if os.path.exists('/data'):
            defaultDir = '/data'
        self.outputPath = QtWidgets.QFileDialog.getExistingDirectory(self, caption="Locate Output Directory", directory=defaultDir)
        self.setDirectories(self.conOutputDir,self.outputPath)

    def browseIndexFile(self):
        defaultDir = '/root'
        if os.path.exists('/data'):
            defaultDir = '/data'
        self.indexFile = QtWidgets.QFileDialog.getOpenFileName(self, "Locate Index File", defaultDir)[0]
        self.setDirectories(self.conIndexDir,self.indexFile)
        
    def browseHostFastqDir(self):
        defaultDir = '/root'
        if os.path.exists('/data'):
            defaultDir = '/data'
        self.hostFastqDir = QtWidgets.QFileDialog.getExistingDirectory(self, caption="Locate Output Directory", directory=defaultDir)
        self.setDirectories(self.conFastqDir,self.hostFastqDir)

    def setRunTrigger(self):
        if self.readyToGo():
            self.startJob

    
