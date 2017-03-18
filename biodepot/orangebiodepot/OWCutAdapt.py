import sys
import os

from Orange.widgets import widget, gui
from orangebiodepot.util.DockerClient import DockerClient, PullImageThread
from PyQt5.QtCore import QThread, pyqtSignal

class OWCutAdapt(widget.OWWidget):
    name = "CutAdapt"
    description = "CutAdapt for RNA-seq"
    icon = "icons/cutadapt.png"

    priority = 10

    inputs = [("CutAdapt", str, "setDirectory")]
    outputs = [("Directory", str)]

    want_main_area = False

    dockerClient = DockerClient('unix:///var/run/docker.sock', 'local')

    image_name = "quay.io/ucsc_cgl/cutadapt"
    image_version = "1.9--6bd44edd2b8f8f17e25c5a268fedaab65fa851d2"

    def __init__(self):
        super().__init__()

        # GUI
        box = gui.widgetBox(self.controlArea, "Info")
        self.infoa = gui.widgetLabel(box, 'Please specify a directory.')
        self.infob = gui.widgetLabel(box, '')
        self.infoa.setWordWrap(True)
        self.infob.setWordWrap(True)

        self.AutoCutAdapt = False
        gui.checkBox(self.controlArea, self, 'AutoCutAdapt', 'Run automatically when directory was set.')

        self.btnCutAdapt = gui.button(self.controlArea, self, "CutAdapt", callback=self.StartCutAdapt)
        self.btnCutAdapt.setEnabled(False)

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
        self.btnCutAdapt.setEnabled(self.bDirectorySet)

        if self.bDirectorySet and self.AutoCutAdapt:
            self.StartCutAdapt()

    def StartCutAdapt(self):
        if not self.bDirectorySet:
            return

        if not self.dockerClient.has_image(self.image_name, self.image_version):
            self.pull_image()
        elif self.bDirectorySet:
            self.run_cutadpat()


    """
        Pull image
    """
    def pull_image(self):
        self.infoa.setText('Pulling \'' + self.image_name + ":" + self.image_version)
        self.setStatusMessage("Downloading...")
        self.progressBarInit()
        self.btnCutAdapt.setEnabled(False)
        # Pull the image in a new thread
        self.pull_image_thread = PullImageThread(self.dockerClient, self.image_name, self.image_version)
        self.pull_image_thread.pull_progress.connect(self.pull_image_progress)
        self.pull_image_thread.finished.connect(self.pull_image_done)
        self.pull_image_thread.start()

    def pull_image_progress(self, val=0):
        self.progressBarSet(val)

    def pull_image_done(self):
        self.infoa.setText('Finished pulling \'' + self.image_name + ":" + self.image_version + '\'')
        self.progressBarFinished()
        self.run_cutadpat()

    def run_cutadpat(self):
        self.btnCutAdapt.setEnabled(False)
        self.infoa.setText('Running cutadapt...')
        self.setStatusMessage('Running...')
        self.progressBarInit()
        # Run the container in a new thread
        self.run_cutadpat_thread = CutAdaptThread(self.dockerClient,
                                                     self.image_name,
                                                     self.image_version,
                                                     self.inputDirectory,
                                                     self.inputDirectory)

        self.run_cutadpat_thread.analysis_progress.connect(self.run_cutadpat_progress)
        self.run_cutadpat_thread.finished.connect(self.run_cutadpat_done)
        self.run_cutadpat_thread.start()

    def run_cutadpat_progress(self, val):
        self.progressBarSet(val)

    def run_cutadpat_done(self):
        self.infoa.setText("Finished running analysis!")
        self.btnCutAdapt.setEnabled(True)
        self.btnCutAdapt.setText('Run again')
        self.setStatusMessage('Finished!')
        self.progressBarFinished()
        self.send("Directory", self.inputDirectory)
        self.btnCutAdapt.setEnabled(True)


class CutAdaptThread(QThread):
    analysis_progress = pyqtSignal(int)

    container_cutadpat_dir = '/data'
    container_results_dir = '/data'

    # quay.io/ucsc_cgl/cutadapt:1.9--6bd44edd2b8f8f17e25c5a268fedaab65fa851d2
    # Forward adapter: AGATCGGAAGAG
    # Reverse adapater: AGATCGGAAGAG

    parameters = ['-a', "AGATCGGAAGAG",
                  '-m', '25']

    parameters.extend(['-A', "AGATCGGAAGAG",
                       '-o', '/data/R1_cutadapt.fastq',
                       '-p', '/data/R2_cutadapt.fastq',
                       '/data/R1.fastq', '/data/R2.fastq'])

    def __init__(self, cli, image_name, image_version, host_cutadpat_dir, host_results_dir):
        QThread.__init__(self)
        self.docker = cli
        self.image_name = image_name
        self.image_version = image_version
        self.host_cutadpat_dir = host_cutadpat_dir
        self.host_results_dir = host_results_dir

    def __del__(self):
        self.wait()

    """
    Run should first create a container and then start it
    """

    def run(self):
        commands = ' '.join((str(w) for w in self.parameters))
        print(commands)

        volumes = {self.host_cutadpat_dir: self.container_cutadpat_dir}
        response = self.docker.create_container(self.image_name+":"+self.image_version,
                                                volumes=volumes,
                                                commands=commands)
        if response['Warnings'] == None:
            self.containerId = response['Id']
            self.docker.start_container(self.containerId)
        else:
            print(response['Warnings'])
        i = 1
        # Keep running until container is exited
        while self.docker.container_running(self.containerId):
            self.analysis_progress.emit(i)
            self.sleep(0.0001)
            i += 1
        # Remove the container now that it is finished
        self.docker.remove_container(self.containerId)

def main(argv=sys.argv):
    from AnyQt.QtWidgets import QApplication
    app = QApplication(list(argv))

    ow = OWCutAdapt()
    ow.setDirectory("/Users/Jimmy/Downloads/Seqs/")
    ow.show()
    ow.raise_()

    ow.handleNewSignals()
    app.exec_()
    ow.set_data(None)
    ow.handleNewSignals()
    return 0

if __name__=="__main__":
    sys.exit(main())