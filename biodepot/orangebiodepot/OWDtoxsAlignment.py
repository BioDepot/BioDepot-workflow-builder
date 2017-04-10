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

        self.host_ref_dir = None
        self.host_seq_dir = None
        self.ref_dir_set = False
        self.seq_dir_set = False

        self.result_folder_name = "/Results"
        self.aligns_folder_name = "/Aligns"

        # TODO is this an issue if multiple containers write to the same place?
        # TODO add timestamp to directory name to differentiate runs
        # Jimmy, Mar/7/2017, I don't know what does "toHostDir" do. Since our DockerClient missed that code, I comment it temporary
        #self.host_counts_dir = counts_dir #self.docker.toHostDir(counts_dir)

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
            self.host_ref_dir = path
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
            self.host_seq_dir = path

            # Jimmy March-28-2017, once the seq input was set, automatically create a result and aligns folder as a sibling of "Seqs"
            parent_path = os.path.abspath(os.path.join(self.host_seq_dir, '..'))
            self.host_counts_dir = os.path.join(parent_path + self.result_folder_name)
            self.host_aligns_dir = os.path.join(parent_path + self.aligns_folder_name)

            if not os.path.exists(self.host_counts_dir):
                os.makedirs(self.host_counts_dir)
            if not os.path.exists(self.host_aligns_dir):
                os.makedirs(self.host_aligns_dir)

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
        self.pull_image_worker = PullImageThread(self.docker, self.image_name, self.image_version)
        self.pull_image_worker.pull_progress.connect(self.pull_image_progress)
        self.pull_image_worker.finished.connect(self.pull_image_finished)
        self.pull_image_worker.start()

    def pull_image_progress(self, val):
        self.progressBarSet(val)

    def pull_image_finished(self):
        self.infoLabel.setText('Finished pulling \'' + self.image_name + ":" + self.image_version + '\'')
        self.progressBarFinished()
        self.run_container()

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
        self.run_container_thread = RunAlignmentThread(self.docker,
                                                       self.image_name,
                                                       self.host_ref_dir,
                                                       self.host_seq_dir,
                                                       self.host_counts_dir,
                                                       self.host_aligns_dir)
        self.run_container_thread.progress.connect(self.run_container_progress)
        self.run_container_thread.finished.connect(self.run_container_finished)
        self.run_container_thread.start()

    def run_container_progress(self, val):
        self.progressBarSet(val)

    def run_container_finished(self):
        self.infoLabel.setText("Finished running alignment!")
        self.btn_run.setEnabled(True)
        self.is_running = False
        self.btn_run.setText('Run again')
        self.setStatusMessage('Finished!')
        self.progressBarFinished()
        self.send("Counts", self.host_counts_dir)


"""
Run Alignment Worker
"""
class RunAlignmentThread(QThread):
    progress = pyqtSignal(int)


    container_ref_dir = "/root/LINCS/References"
    container_seq_dir = "/root/LINCS/Seqs"
    container_counts_dir = "/root/LINCS/Counts"
    container_aligns_dir = "/root/LINCS/Aligns"

    commands = ["/root/LINCS/Programs/Broad-DGE/run-alignment-analysis.sh >& /root/LINCS/Counts/run-alignment-analysis.log; "
                "exit"]

    def __init__(self, cli, image_name, host_ref_dir, host_seq_dir, host_counts_dir, host_aligns_dir):
        QThread.__init__(self)

        self.docker = cli
        self.image_name = image_name
        self.host_ref_dir = host_ref_dir
        self.host_seq_dir = host_seq_dir
        self.host_counts_dir = host_counts_dir
        self.host_aligns_dir = host_aligns_dir
        self.containerId = ""

    def __del__(self):
        self.wait()

    def run(self):
        volumes = {self.host_ref_dir: self.container_ref_dir,
                   self.host_seq_dir: self.container_seq_dir,
                   self.host_counts_dir: self.container_counts_dir,
                   self.host_aligns_dir: self.container_aligns_dir}

        response = self.docker.create_container(self.image_name,
                                                volumes=volumes,
                                                commands=self.commands)
        if response['Warnings'] == None:
            self.containerId = response['Id']
            self.docker.start_container(self.containerId)
        else:
            print(response['Warnings'])

        i = 1
        # Keep running until container is exited
        while self.docker.container_running(self.containerId):
            # self.docker.printStats(self.containerId)
            self.progress.emit(i)
            time.sleep(2)
            i += 2
        # Remove the container now that it is finished
        self.docker.remove_container(self.containerId)


if __name__ == "__main__":
    import sys
    from AnyQt.QtWidgets import QApplication

    a = QApplication(sys.argv)
    ow = OWDtoxsAlignment()
    ow.show()
    a.exec_()
    ow.saveSettings()
