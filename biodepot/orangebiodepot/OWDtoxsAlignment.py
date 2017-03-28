import os, time
from PyQt5.QtWidgets import *

from Orange.widgets import widget, gui
from Orange.widgets.settings import Setting
from AnyQt.QtCore import QObject, QThread, pyqtSignal, pyqtSlot

from orangebiodepot.util.DockerClient import DockerClient, PullImageThread

class OWDtoxsAlignment(widget.OWWidget):
    name = "Dtoxs Alignment"
    description = "Step 1 of the Dtoxs SOP. Uses a Burrows-Wheeler Aligner (BWA)."
    category = "RNASeq"
    icon = "icons/dtoxs-alignment2.svg"
    priority = 10

    image_name = "biodepot/dtoxs_alignment"
    image_version = "latest"

    inputs = [("References", str, "set_refs"),
              ("Seqs", str, "set_seqs")]
    outputs = [("Counts", str)]

    auto_run = Setting(True)

    want_main_area = False

    def __init__(self):
        super().__init__()

        # Docker Client
        # This client talks to your local docker
        self.docker = DockerClient('unix:///var/run/docker.sock', 'local')

        # The directory of the seq data needed to run
        # the container will be set by the user before the
        # container can be run
        self.host_ref_dir = None
        self.host_seq_dir = None
        self.ref_dir_set = False
        self.seq_dir_set = False

        # The default write location of all widgets should be
        # ~/BioDepot/WidgetName/
        # TODO is this an issue if multiple containers write to the same place?
        # TODO add timestamp to directory name to differentiate runs
        counts_dir = '~/BioDepot/Dtoxs_Alignment/Counts'
        if not os.path.exists(counts_dir):
            os.makedirs(counts_dir)
        # Jimmy, Mar-7-2-17, I don't know what does "toHostDir" do. Since our DockerClient missed that code, I comment it temporary
        self.host_counts_dir = counts_dir #self.docker.toHostDir(counts_dir)

        # GUI
        box = gui.widgetBox(self.controlArea, "Info")
        self.infoLabel = gui.widgetLabel(box, 'Connect to a Directory widget '
                                              'to specify location of reference '
                                              'and seq fastq files')
        self.infoLabel.setWordWrap(True)

        gui.checkBox(self.controlArea, self, 'auto_run', 'Run automatically when input set')
        self.btn_run = gui.button(self.controlArea, self, "Run", callback=self.btn_run_pushed)
        self.is_running = False


    """
    Called when the user pushes 'Run'
    """
    def btn_run_pushed(self):
        self.start_container(run_btn_pushed=True)

    """
    Set references
    """
    def set_refs(self, path):
        # When a user removes a connected Directory widget,
        # it sends a signal with path=None
        if path is None:
            self.ref_dir_set = False
        else:
            # Jimmy, Mar-7-2-17, I don't know what does "toHostDir" do. Since our DockerClient missed that code, I comment it temporary
            self.host_ref_dir = path # self.docker.toHostDir(path)
            if self.host_ref_dir is None:
                # TODO emit error
                self.ref_dir_set = False
                print('References set to invalid directory')
            else:
                self.ref_dir_set = True
        self.start_container(run_btn_pushed=False)

    """
    Set seqs
    """
    def set_seqs(self, path):
        # When a user removes a connected Directory widget,
        # it sends a signal with path=None
        if path is None:
            self.seq_dir_set = False
        else:
            # Jimmy, Mar-7-2-17, I don't know what does "toHostDir" do. Since our DockerClient missed that code, I comment it temporary
            self.host_seq_dir = path # self.docker.toHostDir(path)
            if self.host_seq_dir is None:
                # TODO emit error
                self.seq_dir_set = False
                print('Seq set to invalid directory')
            else:
                self.seq_dir_set = True
        self.start_container(run_btn_pushed=False)


    """
    Pull image
    """
    def pull_image(self):
        self.infoLabel.setText('Pulling \'' + self.image_name + ":" + self.image_version + '\' from Dockerhub...')
        self.setStatusMessage("Downloading...")
        self.progressBarInit()
        self.is_running = True
        self.btn_run.setEnabled(False)
        # Pull the image in a new thread
        self.pull_image_thread = QThread()
        self.pull_image_worker = PullImageThread(self.docker, self.image_name, self.image_version)
        self.pull_image_worker.pull_progress[int].connect(self.pull_image_progress)
        self.pull_image_worker.finished.connect(self.pull_image_finished)
        self.pull_image_worker.moveToThread(self.pull_image_thread)
        self.pull_image_thread.started.connect(self.pull_image_worker.work)
        self.pull_image_thread.start()

    @pyqtSlot(int, name="pullImageProgress")
    def pull_image_progress(self, val):
        self.progressBarSet(val)

    @pyqtSlot(name="pullImageFinished")
    def pull_image_finished(self):
        self.pull_image_thread.terminate()
        self.pull_image_thread.wait()
        self.infoLabel.setText('Finished pulling \'' + self.image_name + ":" + self.image_version + '\'')
        self.progressBarFinished()
        self.start_container()

    """
    Alignment
    """
    def start_container(self, run_btn_pushed=True):
        # Make sure both inputs are set
        if self.ref_dir_set and self.seq_dir_set:
            if not self.is_running and (self.auto_run or run_btn_pushed):
                # Make sure the docker image is downloaded
                if not self.docker.has_image(self.image_name, self.image_version):
                    self.pull_image()
                else:
                    self.run_container()
            else:
                self.infoLabel.setText('References and Seqs directories set.\nWaiting to run...')
        elif self.ref_dir_set:
            self.infoLabel.setText("Waiting for user to set Seqs directory.")
        elif self.seq_dir_set:
            self.infoLabel.setText("Waiting for user to set References directory.")


    def run_container(self):
        self.is_running = True
        self.infoLabel.setText('Running alignment...')
        self.setStatusMessage('Running...')
        self.progressBarInit()
        # Run the container in a new thread
        self.run_container_thread = QThread()
        self.run_container_worker = RunAlignmentWorker(self.docker,
                                                       self.image_name,
                                                       self.host_ref_dir,
                                                       self.host_seq_dir,
                                                       self.host_counts_dir)
        self.run_container_worker.progress[int].connect(self.run_container_progress)
        self.run_container_worker.finished.connect(self.run_container_finished)
        self.run_container_worker.moveToThread(self.run_container_thread)
        self.run_container_thread.started.connect(self.run_container_worker.work)
        self.run_container_thread.start()

    @pyqtSlot(int, name="runContainerProgress")
    def run_container_progress(self, val):
        self.progressBarSet(val)

    @pyqtSlot(name="runContainerFinished")
    def run_container_finished(self):
        self.run_container_thread.terminate()
        self.run_container_thread.wait()
        self.infoLabel.setText("Finished running alignment!")
        self.btn_run.setEnabled(True)
        self.is_running = False
        self.btn_run.setText('Run again')
        self.setStatusMessage('Finished!')
        self.progressBarFinished()
        self.send("Counts", self.docker.toContainerDir(self.host_counts_dir))
        # TODO create Orange.data.Table from TOP-40.tsv


"""
Run Alignment Worker
"""
class RunAlignmentWorker(QObject):
    progress = pyqtSignal(int, name="runContainerProgress")
    finished = pyqtSignal(name="runContainerFinished")

    container_ref_dir = "/root/LINCS/References"
    container_seq_dir = "/root/LINCS/Seqs"
    container_counts_dir = "/root/LINCS/Counts"

    commands = ["cd /root/LINCS/; "
                "Programs/Broad-DGE/run-alignment-analysis.sh >& run-alignment-analysis.log; "
                "exit"]

    def __init__(self, cli, image_name, host_ref_dir, host_seq_dir, host_counts_dir):
        super().__init__()
        self.docker = cli
        self.image_name = image_name
        self.host_ref_dir = host_ref_dir
        self.host_seq_dir = host_seq_dir
        self.host_counts_dir = host_counts_dir
        self.containerId = ""

    def work(self):
        volumes = {self.host_ref_dir: self.container_ref_dir,
                   self.host_seq_dir: self.container_seq_dir,
                   self.host_counts_dir: self.container_counts_dir}
        response = self.docker.create_container(self.image_name,
                                                volumes=volumes,
                                                commands=self.commands)
        if response['Warnings'] is None:
            self.containerId = response['Id']
            self.docker.start_container(self.containerId)
        else:
            # TODO emit warning to Widget
            print(response['Warnings'])

        i = 1
        # Keep running until container is exited
        while self.docker.isRunning(self.containerId):
            # self.docker.printStats(self.containerId)
            self.progress.emit(i)
            time.sleep(2)
            i += 2
        # Remove the container now that it is finished
        self.docker.remove_container(self.containerId)
        self.finished.emit()


if __name__ == "__main__":
    import sys
    from AnyQt.QtWidgets import QApplication

    a = QApplication(sys.argv)
    ow = OWDtoxsAlignment()
    ow.show()
    a.exec_()
    ow.saveSettings()
