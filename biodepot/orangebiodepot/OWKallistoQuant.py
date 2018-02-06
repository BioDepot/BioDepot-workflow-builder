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
              ("hostIndexFile", str, "setHostIndexFile")]
    outputs = [("hostOutputDir", str)]
    
    #persistent settings
    hostOutputDir = settings.Setting('/data', schema_only=True)
    hostIndexFile = settings.Setting('/data/kallisto_hg38.idx', schema_only=True)
    hostFastqDir = settings.Setting('/data/fastq', schema_only=True)
#    pseudoBam = settings.Setting(False, schema_only=True)
    nThreads = settings.Setting(1, schema_only=True)
    nThreadsChecked = settings.Setting(False, schema_only=True)
    nBootstraps=settings.Setting(0, schema_only=True)
    nBootstrapsChecked = settings.Setting(False, schema_only=True)
    hostIndexFileConnected = settings.Setting(False, schema_only=True)
    hostFastqDirConnected = settings.Setting(False, schema_only=True)
    
    def __init__(self):
        super().__init__(self.docker_image_name, self.docker_image_tag)
        self.conOutputDir="/root/output"
        self.conFastqDir="/root/fastq"
        self.conIndexDir="/root/index"
        #for attr in ('hostFastqDir','hostIndexFile','hostOutputDir'):
        #    path=getattr(self,attr)
        #    setattr(self,attr,self.dockerClient.to_host_directory(path))
        self.setDirectories(self.conOutputDir,self.hostOutputDir)
        self.drawGUI()

    def drawGUI(self):
        controlBox = gui.widgetBox(self.controlArea, "Status")
        self.infoLabel = gui.widgetLabel(controlBox, 'Waiting...')
        self.infoLabel.setWordWrap(True)
        optionsBox = gui.widgetBox(self.controlArea, "Options")
        outputBox = gui.vBox(optionsBox, "Choose output directory")
        gui.lineEdit(outputBox, self, 'hostOutputDir')
        browseOutputBtn=gui.button(outputBox, self, "☰ Browse", callback=self.browseOutputDir, autoDefault=True, width=80)
        
        self.hostIndexFileBox = gui.vBox(optionsBox, "Choose index file",disabled=self.hostIndexFileConnected)
        self.hostIndexFileLedit=gui.lineEdit(self.hostIndexFileBox, self, 'hostIndexFile',disabled=self.hostIndexFileConnected)
        self.browseIndexBtn=gui.button(self.hostIndexFileBox, self, "☰ Browse", callback=self.browseIndexFile, autoDefault=True, width=80, disabled=self.hostIndexFileConnected)
        
        self.hostFastqDirBox = gui.vBox(optionsBox, "Choose fastq directory",disabled=self.hostFastqDirConnected)
        self.hostFastqDirLedit=gui.lineEdit(self.hostFastqDirBox, self, 'hostFastqDir',disabled=self.hostFastqDirConnected)
        self.browseHostFastqDirBtn=gui.button(self.hostFastqDirBox, self, "☰ Browse", callback=self.browseHostFastqDir, autoDefault=True, width=80, disabled=self.hostFastqDirConnected)
        
        parmBox = gui.vBox(optionsBox, "Optional flags")
#        gui.checkBox(parmBox, self, "pseudoBam", "Output pseudoBam files")
        gui.spin(parmBox, self, "nThreads", minv=1, maxv=128,
        label="Nthreads:", checked='nThreadsChecked')
        gui.spin(parmBox, self, "nBootstraps", minv=1, maxv=128,
        label="Bootstrap samples:", checked='nBootstrapsChecked')
        btnRun = gui.button(optionsBox, self, "Run", callback=self.OnRunClicked)

    def OnRunClicked(self):
        self.startJob()

    def OnSettingsChanged(self):
        self.drawGui()

    def startJob(self):
        volumes = {}
        for conDir in (self.conOutputDir,self.conFastqDir,self.conIndexDir):
            if not self.getDirectory(conDir):
                self.infoLabel.setText("Incorrect mapping of directory to " + conDir)
                return
            volumes[conDir]=self.getDirectory(conDir)
        
        flags={"-i":self.conIndexFile,"-o":self.conOutputDir}
        #if self.pseudoBam :
            #flags['--pseudobam']=None
        if self.nThreadsChecked :
            flags['-t']=str(self.nThreads)
        if self.nBootstrapsChecked :
            flags['-b']=str(self.nBootstraps)
        logFile=os.path.join(self.conOutputDir,'kallistoLog')
        bareCmd=self.generateCmd(('kallisto', 'quant'),flags,args=[])
        cmd="bash -c 'find "+ self.conFastqDir+"/*.fastq.* -type f | sort | xargs " + bareCmd + " >& {}'".format(logFile)
#        self.infoLabel.setText(cmd)
        self.dockerRun(volumes,cmd)

    def Event_OnRunFinished(self):
        self.infoLabel.setText("Finished")
        self.send("hostOutputDir", self.hostOutputDir)
        #self.send("runTrigger", "done")
        
    def Event_OnRunMessage(self, message):
        self.infoLabel.setText(message)

    def setHostFastqDir(self, path, sourceId=None):
        #This only works for single connections
        #change to Dict to keep track of connections
        if path is None:
            self.hostFastqDirConnected = False
            self.hostFastqDir = None
        else:
            #self.hostFastqDir=self.dockerClient.to_host_directory(path)
            self.hostFastqDir=path
            self.hostFastqDirConnected = True
        self.setDirectories(self.conFastqDir,self.hostFastqDir)
        self.browseHostFastqDirBtn.setEnabled(not self.hostFastqDirConnected)
        self.hostFastqDirLedit.setEnabled(not self.hostFastqDirConnected)
        self.hostFastqDirBox.setEnabled(not self.hostFastqDirConnected)


    def setHostIndexFile(self, path, sourceId=None):
        if path is None:
            self.hostIndexFileConnected = False
            self.hostIndexFile = None
            self.hostIndexDir = None
        else:
            #self.hostIndexFile=self.dockerClient.to_host_directory(path)
            self.hostIndexFile=path
            self.hostIndexDir=os.path.dirname(self.hostIndexFile)
            self.hostIndexFileConnected = True
        self.setDirectories(self.conIndexDir,self.hostIndexDir)
        self.browseIndexBtn.setEnabled(not self.hostIndexFileConnected)
        self.hostIndexFileLedit.setEnabled(not self.hostIndexFileConnected)
        self.hostIndexFileBox.setEnabled(not self.hostIndexFileConnected)
        self.conIndexFile=os.path.normpath(str.join(os.sep,(self.conIndexDir,os.path.basename(self.hostIndexFile))))
        
    def browseOutputDir(self):
        defaultDir = '/root'
        if os.path.exists('/data'):
            defaultDir = '/data'
        self.hostOutputDir = QtWidgets.QFileDialog.getExistingDirectory(self, caption="Locate Output Directory", directory=defaultDir)
        #self.hostOutputDir=self.dockerClient.to_host_directory(self.hostOutputDir)
        self.setDirectories(self.conOutputDir,self.hostOutputDir)

    def browseIndexFile(self):
        defaultDir = '/root'
        if os.path.exists('/data'):
            defaultDir = '/data'
        self.hostIndexFile = QtWidgets.QFileDialog.getOpenFileName(self, "Locate Index File", defaultDir)[0]
        #self.hostIndexFile =self.dockerClient.to_host_directory(self.hostIndexFile)
        self.hostIndexDir=os.path.dirname(self.hostIndexFile)
        self.setDirectories(self.conIndexDir,self.hostIndexDir)
        self.conIndexFile=os.path.normpath(str.join(os.sep,(self.conIndexDir,os.path.basename(self.hostIndexFile))))
        
    def browseHostFastqDir(self):
        defaultDir = '/root'
        if os.path.exists('/data'):
            defaultDir = '/data'
        self.hostFastqPath = QtWidgets.QFileDialog.getExistingDirectory(self, caption="Locate Fastq Directory", directory=defaultDir)
        #self.hostFastqPath=self.dockerClient.to_host_directory(self.hostFastqPath)
        self.setDirectories(self.conFastqDir,self.hostFastqPath)

    def setRunTrigger(self):
        if self.readyToGo():
            self.startJob

    
