import os
import glob
from Orange.widgets import widget, gui, settings
from orangebiodepot.util.DockerClient import DockerClient
from orangebiodepot.util.BwBase import OWBwBWidget, ConnectionDict
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
    #this is a dict of lists of form {inputSlot:(widgetId1,widgetId2..widgetIdN) where None is a valid entry for the case when there are no Multiple connections }
    inputConnectionsStore = settings.Setting({},schema_only=True);
    
    def __init__(self):
        super().__init__(self.docker_image_name, self.docker_image_tag)
        self.conOutputDir="/root/output"
        self.conFastqDir="/root/fastq"
        self.conIndexDir="/root/index"
        #for attr in ('hostFastqDir','hostIndexFile','hostOutputDir'):
        #    path=getattr(self,attr)
        #    setattr(self,attr,self.dockerClient.to_host_directory(path))
        self.inputConnections=ConnectionDict(self.inputConnectionsStore)
        self.setDirectories(self.conOutputDir,self.hostOutputDir)
        self.defaultFileIcon=QtGui.QIcon('icons/file.png')
        self.drawGUI()

    def drawGUI(self):
        controlBox = gui.widgetBox(self.controlArea, "Status")
        self.infoLabel = gui.widgetLabel(controlBox, 'Waiting...')
        self.infoLabel.setWordWrap(True)
        optionsBox = gui.widgetBox(self.controlArea, "Options")
        self.outputLedit=gui.lineEdit(None, self, 'hostOutputDir')
        self.browseOutputBtn=gui.button(None, self, "", callback=self.browseOutputDir, autoDefault=True, width=19, height=19)
        self.bwbFileEntry(optionsBox,self.browseOutputBtn,self.outputLedit, label='Output directory:', entryType='directory')
        
        self.hostIndexFileLedit=gui.lineEdit(None, self, 'hostIndexFile',disabled=self.inputConnections.isConnected('hostIndexFile'))
        self.browseIndexBtn=gui.button(None, self, "", callback=self.browseIndexFile, autoDefault=True,  width=19, height=19, disabled=self.inputConnections.isConnected('hostIndexFile'))
        self.bwbFileEntry(optionsBox,self.browseIndexBtn, self.hostIndexFileLedit, label='Index file:', entryType='file')
        
        self.hostFastqDirLedit=gui.lineEdit(None, self, 'hostFastqDir',disabled=self.inputConnections.isConnected('hostFastqDir'))
        self.browseHostFastqDirBtn=gui.button(None, self, "", callback=self.browseHostFastqDir, autoDefault=True, width=19, height=19, disabled=self.inputConnections.isConnected('hostFastqDir'))
        self.bwbFileEntry(optionsBox,self.browseHostFastqDirBtn, self.hostFastqDirLedit, label='Fastq directory: ', entryType='directory')      
                           
        parmBox = gui.vBox(optionsBox, " ")
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
        if path is None:
            self.inputConnections.remove('hostFastqDir',sourceId)
            self.hostFastqDir = None
            self.hostFastqDirLedit.clear()
        else:
            self.hostFastqDir=path
            self.inputConnections.add('hostFastqDir',sourceId)
            
        self.setDirectories(self.conFastqDir,self.hostFastqDir)
        
        self.browseHostFastqDirBtn.setEnabled(not self.inputConnections.isConnected('hostFastqDir'))
        self.hostFastqDirLedit.setEnabled(not self.inputConnections.isConnected('hostFastqDir'))

    def setHostIndexFile(self, path, sourceId=None):
        if path is None:
            self.inputConnections.remove('hostIndexFile',sourceId)
            self.hostIndexFile = None
            self.hostIndexDir = None
            self.hostIndexFileLedit.clear()
        else:
            self.hostIndexFile=path
            self.hostIndexDir=os.path.dirname(self.hostIndexFile)
            self.inputConnections.add('hostIndexFile',sourceId)
        self.setDirectories(self.conIndexDir,self.hostIndexDir)
        
        self.browseIndexBtn.setEnabled(not self.inputConnections.isConnected('hostIndexFile'))
        self.hostIndexFileLedit.setEnabled(not self.inputConnections.isConnected('hostIndexFile'))
        self.conIndexFile=os.path.normpath(str.join(os.sep,(self.conIndexDir,os.path.basename(self.hostIndexFile))))
        
    def browseOutputDir(self):
        defaultDir = '/root'
        if os.path.exists('/data'):
            defaultDir = '/data'
        self.hostOutputDir = QtWidgets.QFileDialog.getExistingDirectory(self, caption="Locate Output Directory", directory=defaultDir)
        self.setDirectories(self.conOutputDir,self.hostOutputDir)

    def browseIndexFile(self):
        defaultDir = '/root'
        if os.path.exists('/data'):
            defaultDir = '/data'
        self.hostIndexFile = QtWidgets.QFileDialog.getOpenFileName(self, "Locate Index File", defaultDir)[0]
        self.hostIndexDir=os.path.dirname(self.hostIndexFile)
        self.setDirectories(self.conIndexDir,self.hostIndexDir)
        self.conIndexFile=os.path.normpath(str.join(os.sep,(self.conIndexDir,os.path.basename(self.hostIndexFile))))
        
    def browseHostFastqDir(self):
        defaultDir = '/root'
        if os.path.exists('/data'):
            defaultDir = '/data'
        self.hostFastqDir = QtWidgets.QFileDialog.getExistingDirectory(self, caption="Locate Fastq Directory", directory=defaultDir)
        self.setDirectories(self.conFastqDir,self.hostFastqDir)

    def setRunTrigger(self):
        if self.readyToGo():
            self.startJob

            

    
