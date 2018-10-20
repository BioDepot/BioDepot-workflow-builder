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

class OWstarAlign(OWBwBWidget):
    name = "starAlign"
    description = "Star aligner alignment module"
    priority = 12
    icon = getIconName(__file__,"staralign.png")
    want_main_area = False
    docker_image_name = "biodepot/star"
    docker_image_tag = "2.6.0c__debian-8.11-slim__072918"
    inputs = [("trigger",str,"handleInputstrigger"),("outputDir",str,"handleInputsoutputDir"),("genomeDir",str,"handleInputsgenomeDir")]
    outputs = [("outputDir",str),("genomeDir",str)]
    pset=functools.partial(settings.Setting,schema_only=True)
    runMode=pset(0)
    exportGraphics=pset(False)
    runTriggers=pset([])
    triggerReady=pset({})
    inputConnectionsStore=pset({})
    optionsChecked=pset({})
    readFilesIn=pset([])
    apelist=pset(False)
    spelist=pset(False)
    parametersFiles=pset(None)
    sysShell=pset("/bin/sh")
    runThreadN=pset(1)
    runRNGseed=pset(777)
    readFilesCommand=pset(None)
    readMapNumber=pset(-1)
    readMatesLengthsIn=pset("NotEqual")
    readNameSeparator=pset("/")
    clip3pNbases=pset(['0'])
    clip5pNbases=pset(['0'])
    clip3pAdapterSeq=pset([])
    clip3pAdapterMMp=pset(['0.1'])
    clip3pAfterAdapterNbases=pset(['0'])
    limitIObufferSize=pset(150000000)
    limitOutSAMoneReadBytes=pset(100000)
    outputDir=pset("")
    outTmpDir=pset("")
    outStd=pset("Log")
    outReadsUnmapped=pset("")
    outMultimapperOrder=pset("Old_2.4")
    outSAMtype=pset("SAM")
    outSAMmode=pset("Full")
    outSAMstrandField=pset(None)
    outSAMattributes=pset("Standard")
    outSAMattrIHstart=pset(1)
    outSAMunmapped=pset(None)
    outSAMorder=pset("Paired")
    outSAMprimaryFlag=pset("OneBestScore")
    outSAMreadID=pset("Standard")
    outSAMmapqUnique=pset(255)
    outSAMflagOR=pset(0)
    outSAMflagAND=pset(65535)
    outSAMattrRGline=pset(None)
    outSAMheaderHD=pset(None)
    outSAMheaderPG=pset(None)
    outSAMheaderCommentFile=pset(None)
    outSAMfilter=pset(None)
    outSAMmultNmax=pset(-1)
    outBAMcompression=pset(1)
    outBAMsortingThreadN=pset(0)
    bamRemoveDuplicatesType=pset(None)
    bamRemoveDuplicatesMate2basesN=pset(0)
    outWigType=pset("")
    outWigStrand=pset("Stranded")
    outWigReferencesPrefix=pset(None)
    outWigNorm=pset("RPM")
    outFilterType=pset("Normal")
    outFilterMultimapScoreRange=pset(1)
    outFilterMultimapNmax=pset(10)
    outFilterMismatchNmax=pset(10)
    outFilterMismatchNoverLmax =pset(0.3)
    outFilterMismatchNoverReadLmax=pset(1.0)
    outFilterScoreMin=pset(0)
    outFilterScoreMinOverLread=pset(0.66)
    outFilterMatchNmin=pset(0)
    outFilterMatchNminOverLread=pset(0.66)
    outFilterIntronMotifs=pset(None)
    scoreGap=pset(0)
    scoreGapNoncan=pset(-8)
    scoreGapGCAG=pset(-4)
    scoreGapATAC=pset(-8)
    scoreGenomicLengthLog2scale=pset(-0.25)
    scoreDelOpen=pset(-2)
    scoreDelBase=pset(-2)
    scoreInsOpen=pset(-2)
    scoreInsBase=pset(-2)
    scoreStitchSJshift=pset(1)
    seedSearchStartLmax=pset(50)
    seedSearchStartLmaxOverLread=pset(1.0)
    seedSearchLmax=pset(0)
    seedMultimapNmax=pset(10000)
    seedPerReadNmax=pset(1000)
    seedPerWindowNmax=pset(50)
    seedNoneLociPerWindow=pset(10)
    alignIntronMin=pset(21)
    alignIntronMax=pset(0)
    alignMatesGapMax=pset(0)
    alignSJoverhangMin=pset(5)
    alignSJstitchMismatchNmax=pset(None)
    alignSJDBoverhangMin=pset(3)
    alignSplicedMateMapLmin=pset(0)
    alignSplicedMateMapLminOverLmate=pset(0.66)
    alignWindowsPerReadNmax=pset(10000)
    alignTranscriptsPerWindowNmax=pset(100)
    alignTranscriptsPerReadNmax=pset(10000)
    alignEndsType=pset("Local")
    alignEndsProtrude=pset(['0', 'ConcordantPair'])
    alignSoftClipAtReferenceEnds=pset("Yes")
    winAnchorMultimapNmax=pset(50)
    winBinNbits=pset(16)
    winAnchorDistNbins=pset(9)
    winFlankNbins=pset(4)
    winReadCoverageRelativeMin=pset(0.5)
    winReadCoverageBasesMin=pset(0)
    chimOutType=pset("SeparateSAMold")
    chimSegmentMin=pset(0)
    chimScoreMin=pset(0)
    chimScoreDropMax=pset(20)
    chimScoreSeparation=pset(10)
    chimScoreJunctionNonGTAG=pset(-1)
    chimJunctionOverhangMin=pset(20)
    chimSegmentReadGapMax=pset(0)
    chimFilter=pset("banGenomicN")
    quantMode=pset(None)
    quantTranscriptomeBAMcompression=pset(1)
    quantTranscriptomeBan=pset("IndelSoftclipSingleend")
    twopassMode=pset("")
    twopass1readsN =pset(-1)
    genomeDir=pset("/data/GenomeDir")
    genomeLoad=pset("NoSharedMemory")
    outputFilePrefix=pset(None)
    multipleSample=pset(False)
    def __init__(self):
        super().__init__(self.docker_image_name, self.docker_image_tag)
        with open(getJsonName(__file__,"starAlign")) as f:
            self.data=jsonpickle.decode(f.read())
            f.close()
        self.initVolumes()
        self.inputConnections = ConnectionDict(self.inputConnectionsStore)
        self.drawGUI()
    def handleInputstrigger(self, value, *args):
        if args and len(args) > 0: 
            self.handleInputs("trigger", value, args[0][0]), test=args[0][3]))
        else:
            self.handleInputs("inputFile", value, None)
    def handleInputsoutputDir(self, value, *args):
        if args and len(args) > 0: 
            self.handleInputs("outputDir", value, args[0][0]), test=args[0][3]))
        else:
            self.handleInputs("inputFile", value, None)
    def handleInputsgenomeDir(self, value, *args):
        if args and len(args) > 0: 
            self.handleInputs("genomeDir", value, args[0][0]), test=args[0][3]))
        else:
            self.handleInputs("inputFile", value, None)
    def handleOutputs(self):
        outputValue="/data"
        if hasattr(self,"outputDir"):
            outputValue=getattr(self,"outputDir")
        self.send("outputDir", outputValue)
        outputValue="/data/GenomeDir"
        if hasattr(self,"genomeDir"):
            outputValue=getattr(self,"genomeDir")
        self.send("genomeDir", outputValue)
