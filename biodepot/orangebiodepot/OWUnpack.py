import sys
import os
import platform
import subprocess
import re

import Orange.data
from Orange.widgets import widget, gui
from orangebiodepot.util.DockerClient import DockerClient

class OWUnpack(widget.OWWidget):
    name = "UnPack"
    description = "Unpack compressed fastqc files"
    icon = "icons/unpack.png"

    priority = 10

    inputs = [("Unpack", str, "setDirectory", widget.Default)]
    outputs = [("Directory", str)]

    want_main_area = False

    dockerClient = DockerClient('unix:///var/run/docker.sock', 'local')

    def __init__(self):
        super().__init__()

        self.inputDirectory = None
        self.bDirectorySet = False;

        # GUI
        box = gui.widgetBox(self.controlArea, "Info")
        self.infoa = gui.widgetLabel(box, 'Please specify a directory.')
        self.infob = gui.widgetLabel(box, '')
        self.infoa.setWordWrap(True)
        self.infob.setWordWrap(True)

        self.AutoUnpack = False
        gui.checkBox(self.controlArea, self, 'AutoUnpack', 'Run automatically when directory was set.')

        self.mergeFile = True
        gui.checkBox(self.controlArea, self, "mergeFile", "Merge files into R1.fastq and R2.fastq")

        self.btnUnpack = gui.button(self.controlArea, self, "Unpack", callback=self.StartUnpack)
        self.btnUnpack.setEnabled(False)

    """
       Set Unpack Directory
    """
    def setDirectory(self, path):
        # When a user removes a connected Directory widget,
        # it sends a signal with path=None
        if path is None:
            self.bDirectorySet = False
        else:
            if not os.path.exists(path):
                self.bDirectorySet = False
                print('References set to invalid directory')
            else:
                self.inputDirectory = path
                self.bDirectorySet = True
                self.infoa.setText("Directory: {0}".format(self.inputDirectory))
        self.btnUnpack.setEnabled(self.bDirectorySet)

        if self.bDirectorySet and self.AutoUnpack:
            self.StartUnpack()
    """
        Do Unpack
        tar  file:///path/to/sample.tar
        #   fq  paired  file:///path/to/R1.fq.gz,file:///path/to/R2.fq.gz
        #   tar single  http://sample-depot.com/single-end-sample.tar
        #   tar paired  s3://my-bucket-name/directory/paired-sample.tar.gz
        #   fq  single  s3://my-bucket-name/directory/single-end-file.fq

        now, only consider [fq paired]
    """
    def StartUnpack(self):
        if not self.bDirectorySet:
            return

        self.progressBarInit()
        command = ['tar', '-zxf']
        if platform.system() == "Darwin":
            #command = ['gunzip', '-k']
            command = ['gunzip']

        fastq_files = []
        for root, subdir, files in os.walk(self.inputDirectory):
            fastq_files.extend([os.path.join(root, x) for x in files])

        if self.mergeFile:
            r1, r2 = [], []
            # Pattern convention: Look for "R1" / "R2" in the filename, or "_1" / "_2" before the extension
            pattern = re.compile('(?:^|[._-])(R[12]|[12]\.f)')
            for fastq in sorted(fastq_files):
                match = pattern.search(os.path.basename(fastq))
                if not match:
                    print ('Invalid FASTQ file: {0}, ignore it.'.format(fastq))
                    continue
                elif '1' in match.group():
                    r1.append(fastq)
                elif '2' in match.group():
                    r2.append(fastq)

            if len(r1) != len(r2):
                self.error('Check fastq names, uneven number of pairs found.')
                print ('Check fastq names, uneven number of pairs found.\nr1: {}\nr2: {}'.format(r1, r2))
                return

            if len(r1) == 0:
                return

            if platform.system() == "Darwin":
                command = ['gunzip', '-c'] if r1[0].endswith('.gz') and r2[0].endswith('.gz') else ['cat']
            else:
                command = ['zcat'] if r1[0].endswith('.gz') and r2[0].endswith('.gz') else ['cat']

            R1_file = open(os.path.join(self.inputDirectory, 'R1.fastq'), 'a')
            R2_file = open(os.path.join(self.inputDirectory, 'R2.fastq'), 'a')

            for i in range(len(r1)):
                p1 = subprocess.Popen(command + [str(r1[i])], stdout = R1_file)
                p2 = subprocess.Popen(command + [str(r2[i])], stdout = R2_file)
                p1.wait()
                p2.wait()
                self.progressBarSet((i + 1) / len(r1) * 100)

            R1_file.flush()
            R2_file.flush()

        else:
            for index, fastq in enumerate(sorted(fastq_files)):
                print("Unpacking....{0}".format(fastq))
                if not fastq.endswith('.gz'):
                    print ("Not supported file format {0}".format(fastq))
                else:
                    runcommand = command + [fastq]
                    #print (runcommand)
                    processor = subprocess.Popen(runcommand)
                    processor.wait()

                self.progressBarSet((index+1)/len(fastq_files)*100)
                #time.sleep(0.1)

        self.progressBarFinished()
        self.send("Directory", self.inputDirectory)

def main(argv=sys.argv):
    from AnyQt.QtWidgets import QApplication
    app = QApplication(list(argv))
    filename = "iris"

    ow = OWUnpack()
    ow.setDirectory("~/Seqs/")
    ow.show()
    ow.raise_()

    dataset = Orange.data.Table(filename)
    ow.handleNewSignals()
    app.exec_()
    ow.handleNewSignals()
    return 0

if __name__=="__main__":
    sys.exit(main())